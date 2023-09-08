#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/string.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/srgb.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/volumegrid.h>
#include <drjit/dynamic.h>
#include <drjit/texture.h>

NAMESPACE_BEGIN(mitsuba)


/**!
.. _volume-gridvolume:

Grid-based volume data source (:monosp:`gridvolume`)
----------------------------------------------------

.. pluginparameters::

 * - filename
   - |string|
   - Filename of the volume to be loaded

 * - grid
   - :monosp:`VolumeGrid object`
   - When creating a grid volume at runtime, e.g. from Python or C++,
     an existing ``VolumeGrid`` instance can be passed directly rather than
     loading it from the filesystem with :paramtype:`filename`.

 * - data
   - |tensor|
   - Tensor array containing the grid data. This parameter can only be specified
     when building this plugin at runtime from Python or C++ and cannot be
     specified in the XML scene description. The :paramtype:`raw` parameter must
     also be set to :monosp:`true` when using a tensor.
   - |exposed|, |differentiable|

 * - filter_type
   - |string|
   - Specifies how voxel values are interpolated. The following options are
     currently available:

     - ``trilinear`` (default): perform trilinear interpolation.

     - ``nearest``: disable interpolation. In this mode, the plugin
       performs nearest neighbor lookups of volume values.

 * - wrap_mode
   - |string|
   - Controls the behavior of volume evaluations that fall outside of the
     :math:`[0, 1]` range. The following options are currently available:

     - ``clamp`` (default): clamp coordinates to the edge of the volume.

     - ``repeat``: tile the volume infinitely.

     - ``mirror``: mirror the volume along its boundaries.

 * - raw
   - |bool|
   - Should the transformation to the stored color data (e.g. sRGB to linear,
     spectral upsampling) be disabled? You will want to enable this when working
     with non-color, 3-channel volume data. Currently, no plugin needs this option
     to be set to true (Default: false)

 * - to_world
   - |transform|
   - Specifies an optional 4x4 transformation matrix that will be applied to volume coordinates.

 * - accel
   - |bool|
   - Hardware acceleration features can be used in CUDA mode. These features can
     cause small differences as hardware interpolation methods typically have a
     loss of precision (not exactly 32-bit arithmetic). (Default: true)

This class implements access to volume data stored on a 3D grid using a
simple binary exchange format (compatible with Mitsuba 0.6). When appropriate,
spectral upsampling is applied at loading time to convert RGB values to
spectra that can be used in the renderer.
We provide a small `helper utility <https://github.com/mitsuba-renderer/mitsuba2-vdb-converter>`_
to convert OpenVDB files to this format. The format uses a
little endian encoding and is specified as follows:

.. list-table:: Volume file format
   :widths: 8 30
   :header-rows: 1

   * - Position
     - Content
   * - Bytes 1-3
     - ASCII Bytes ’V’, ’O’, and ’L’
   * - Byte 4
     - File format version number (currently 3)
   * - Bytes 5-8
     - Encoding identified (32-bit integer). Currently, only a value of 1 is
       supported (float32-based representation)
   * - Bytes 9-12
     - Number of cells along the X axis (32 bit integer)
   * - Bytes 13-16
     - Number of cells along the Y axis (32 bit integer)
   * - Bytes 17-20
     - Number of cells along the Z axis (32 bit integer)
   * - Bytes 21-24
     - Number of channels (32 bit integer, supported values: 1, 3 or 6)
   * - Bytes 25-48
     - Axis-aligned bounding box of the data stored in single precision (order:
       xmin, ymin, zmin, xmax, ymax, zmax)
   * - Bytes 49-*
     - Binary data of the volume stored in the specified encoding. The data
       are ordered so that the following C-style indexing operation makes sense
       after the file has been loaded into memory:
       :code:`data[((zpos*yres + ypos)*xres + xpos)*channels + chan]`
       where (xpos, ypos, zpos, chan) denotes the lookup location.

.. tabs::
    .. code-tab:: xml

        <medium type="heterogeneous">
            <volume type="grid" name="albedo">
                <string name="filename" value="my_volume.vol"/>
            </volume>
        </medium>

    .. code-tab:: python

        'type': 'heterogeneous',
        'albedo': {
            'type': 'grid',
            'filename': 'my_volume.vol'
        }

*/

/**
 * Interpolated 3D grid texture of scalar or color values.
 *
 * This plugin loads RGB data from a binary file. When appropriate,
 * spectral upsampling is applied at loading time to convert RGB values
 * to spectra that can be used in the renderer.
 *
 * Data layout:
 * The data must be ordered so that the following C-style (row-major) indexing
 * operation makes sense after the file has been mapped into memory:
 *     data[((zpos*yres + ypos)*xres + xpos)*channels + chan]}
 *     where (xpos, ypos, zpos, chan) denotes the lookup location.
 */
template <typename Float, typename Spectrum>
class GridVolume final : public Volume<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Volume, update_bbox, m_to_local, m_bbox, m_channel_count)
    MI_IMPORT_TYPES(VolumeGrid)

    GridVolume(const Properties &props) : Base(props) {
        std::string filter_type_str = props.string("filter_type", "trilinear");
        dr::FilterMode filter_mode;
        if (filter_type_str == "nearest")
            filter_mode = dr::FilterMode::Nearest;
        else if (filter_type_str == "trilinear")
            filter_mode = dr::FilterMode::Linear;
        else
            Throw("Invalid filter type \"%s\", must be one of: \"nearest\" or "
                  "\"trilinear\"!", filter_type_str);

        std::string wrap_mode_st = props.string("wrap_mode", "clamp");
        dr::WrapMode wrap_mode;
        if (wrap_mode_st == "repeat")
            wrap_mode = dr::WrapMode::Repeat;
        else if (wrap_mode_st == "mirror")
            wrap_mode = dr::WrapMode::Mirror;
        else if (wrap_mode_st == "clamp")
            wrap_mode = dr::WrapMode::Clamp;
        else
            Throw("Invalid wrap mode \"%s\", must be one of: \"repeat\", "
                  "\"mirror\", or \"clamp\"!",
                  wrap_mode_st);

        m_raw = props.get<bool>("raw", false);
        m_accel = props.get<bool>("accel", true);

        // Load volume data
        ref<VolumeGrid> volume_grid = nullptr;
        TensorXf* tensor = nullptr;
        {
            ScalarVector3u res;
            ScalarUInt32 channel_count = 0;

            if (props.has_property("grid")) {
                // Creates a Bitmap texture directly from an existing Bitmap object
                if (props.has_property("filename"))
                    Throw("Cannot specify both \"grid\" and \"filename\".");
                Log(Debug, "Loading volume grid from memory...");
                // Note: ref-counted, so we don't have to worry about lifetime
                ref<Object> other = props.object("grid");
                volume_grid = dynamic_cast<VolumeGrid *>(other.get());
                if (!volume_grid)
                    Throw("Property \"grid\" must be a VolumeGrid instance.");
                res = volume_grid->size();
                channel_count = volume_grid->channel_count();
            } else if(props.has_property("data")) {
                tensor = props.tensor<TensorXf>("data");
                if (tensor->ndim() != 3 && tensor->ndim() != 4)
                    Throw("Tensor->has %ul dimensions. Expected 3 or 4", tensor->ndim());
                res = { (uint32_t) tensor->shape(2), (uint32_t) tensor->shape(1), (uint32_t) tensor->shape(0) };
                channel_count = tensor->ndim() == 4 ? tensor->shape(3) : 1;

                if (channel_count != 1 && channel_count != 3 && channel_count != 6)
                    Throw("Tensor shape at index 3 is %lu invalid. Only volumes with 1, 3 or 6 "
                          "channels are supported!", to_string(), channel_count);
            } else {
                FileResolver *fs = Thread::thread()->file_resolver();
                fs::path file_path = fs->resolve(props.string("filename"));
                if (!fs::exists(file_path))
                    Log(Error, "\"%s\": file does not exist!", file_path);
                volume_grid = new VolumeGrid(file_path);
                res = volume_grid->size();
                channel_count = volume_grid->channel_count();
            }

            ScalarUInt32 size = dr::prod(res);

            // Apply spectral conversion if necessary
            if (is_spectral_v<Spectrum> && channel_count == 3 &&
                !m_raw) {
                if (tensor)
                    Throw("Spectral conversion of tensor input is not supported "
                          "and requires a volume grid");

                ScalarFloat *ptr = volume_grid->data();

                auto scaled_data =
                    std::unique_ptr<ScalarFloat[]>(new ScalarFloat[size * 4]);
                ScalarFloat *scaled_data_ptr = scaled_data.get();
                ScalarFloat max = 0.0;
                for (ScalarUInt32 i = 0; i < size; ++i) {
                    ScalarColor3f rgb = dr::load<ScalarColor3f>(ptr);
                    // TODO: Make this scaling optional if the RGB values are
                    // between 0 and 1
                    ScalarFloat scale = dr::max(rgb) * 2.f;
                    ScalarColor3f rgb_norm =
                        rgb / dr::maximum((ScalarFloat) 1e-8, scale);
                    ScalarVector3f coeff = srgb_model_fetch(rgb_norm);
                    max = dr::maximum(max, scale);
                    dr::store(scaled_data_ptr,
                              dr::concat(coeff, dr::Array<ScalarFloat, 1>(scale)));
                    ptr += 3;
                    scaled_data_ptr += 4;
                }
                m_max = (float) max;

                size_t shape[4] = {
                    (size_t) res.z(),
                    (size_t) res.y(),
                    (size_t) res.x(),
                    4
                };
                m_texture = Texture3f(TensorXf(scaled_data.get(), 4, shape),
                                      m_accel, m_accel, filter_mode, wrap_mode);
            } else if (volume_grid) {
                size_t shape[4] = {
                    (size_t) res.z(),
                    (size_t) res.y(),
                    (size_t) res.x(),
                    channel_count
                };
                m_texture = Texture3f(TensorXf(volume_grid->data(), 4, shape),
                                      m_accel, m_accel, filter_mode, wrap_mode);
                m_max = volume_grid->max();
                m_max_per_channel.resize(volume_grid->channel_count());
                volume_grid->max_per_channel(m_max_per_channel.data());
                m_channel_count = channel_count;
            } else if (tensor) {
                size_t shape[4] = {
                    (size_t) res.z(),
                    (size_t) res.y(),
                    (size_t) res.x(),
                    channel_count
                };
                m_texture = Texture3f(TensorXf(tensor->array(), 4, shape),
                                      m_accel, m_accel, filter_mode, wrap_mode);
                m_max = (float) dr::max_nested(dr::detach(m_texture.value()));
                m_channel_count = channel_count;
            }
        }

        if (props.get<bool>("use_grid_bbox", false)) {
            if (tensor)
                Throw("use_grid_bbox is unsupported with tensor input and requires a volume grid");
            m_to_local = volume_grid->bbox_transform() * m_to_local;
            update_bbox();
        }

        if (props.has_property("max_value")) {
            m_fixed_max = true;
            m_max = props.get<ScalarFloat>("max_value");
        }
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_parameter("data", m_texture.tensor(), +ParamFlags::Differentiable);
        Base::traverse(callback);
    }

    void parameters_changed(const std::vector<std::string> &keys) override {
        if (keys.empty() || string::contains(keys, "data")) {
            const size_t channels = nchannels();
            if (channels != 1 && channels != 3 && channels != 6)
                Throw("parameters_changed(): The volume data %s was changed "
                      "to have %d channels, only volumes with 1, 3 or 6 "
                      "channels are supported!", to_string(), channels);

            m_texture.set_tensor(m_texture.tensor());

            if (!m_fixed_max)
                m_max = (float) dr::max_nested(dr::detach(m_texture.value()));
        }
    }

    UnpolarizedSpectrum eval(const Interaction3f &it,
                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);

        const size_t channels = nchannels();
        if (channels == 3 && is_spectral_v<Spectrum> && m_raw)
            Throw("The GridVolume texture %s was queried for a spectrum, but "
                  "texture conversion into spectra was explicitly disabled! "
                  "(raw=true)",
                  to_string());
        else if (channels != 3 && channels != 1)
            Throw("The GridVolume texture %s was queried for a spectrum, but "
                  "has a number of channels which is not 1 or 3",
                  to_string());
        else {
            if (dr::none_or<false>(active))
                return dr::zeros<UnpolarizedSpectrum>();

            if constexpr (is_monochromatic_v<Spectrum>) {
                if (channels == 1)
                    return interpolate_1(it, active);
                else // 3 channels
                    return luminance(interpolate_3(it, active));
            }
            else{
                if (channels == 1)
                    return interpolate_1(it, active);
                else { // 3 channels
                    if constexpr (is_spectral_v<Spectrum>)
                        return interpolate_spectral(it, active);
                    else
                        return interpolate_3(it, active);
                }
            }
        }
    }

    Float eval_1(const Interaction3f &it, Mask active = true) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);

        const size_t channels = nchannels();
        if (channels == 3 && is_spectral_v<Spectrum> && !m_raw)
            Throw("eval_1(): The GridVolume texture %s was queried for a "
                  "scalar value, but texture conversion into spectra was "
                  "requested! (raw=false)",
                  to_string());
        else {
            if (dr::none_or<false>(active))
                return dr::zeros<Float>();

            if (channels == 1)
                return interpolate_1(it, active);
            else if (channels == 3)
                return luminance(interpolate_3(it, active));
            else // 6 channels
                return dr::mean(interpolate_6(it, active));
        }
    }

    void eval_n(const Interaction3f &it, Float *out, Mask active = true) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);

        return interpolate_per_channel<Float>(it, out, active);
    }

    Vector3f eval_3(const Interaction3f &it,
                    Mask active = true) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);

        const size_t channels = nchannels();
        if (channels != 3) {
            Throw("eval_3(): The GridVolume texture %s was queried for a 3D "
                  "vector, but it has %s channel(s)", to_string(), channels);
        } else if (is_spectral_v<Spectrum> && !m_raw) {
            Throw("eval_3(): The GridVolume texture %s was queried for a 3D "
                  "vector, but texture conversion into spectra was requested! "
                  "(raw=false)", to_string());
        } else {
            if (dr::none_or<false>(active))
                return dr::zeros<Vector3f>();

            return interpolate_3(it, active);
        }
    }

    dr::Array<Float, 6> eval_6(const Interaction3f &it,
                               Mask active = true) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);

        const size_t channels = nchannels();
        if (channels != 6)
            Throw("eval_6(): The GridVolume texture %s was queried for a 6D "
                  "vector, but it has %s channel(s)", to_string(), channels);
        else {
            if (dr::none_or<false>(active))
                return dr::zeros<dr::Array<Float, 6>>();

            return interpolate_6(it, active);
        }
    }

    ScalarFloat max() const override { return m_max; }

    void max_per_channel(ScalarFloat *out) const override {
        for (size_t i=0; i<m_max_per_channel.size(); ++i)
            out[i] = m_max_per_channel[i];
    }

    ScalarVector3i resolution() const override {
        const size_t *shape = m_texture.shape();
        return { (int) shape[2], (int) shape[1], (int) shape[0] };
    };

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "GridVolume[" << std::endl
            << "  to_local = " << string::indent(m_to_local, 13) << "," << std::endl
            << "  bbox = " << string::indent(m_bbox) << "," << std::endl
            << "  dimensions = " << resolution() << "," << std::endl
            << "  max = " << m_max << "," << std::endl
            << "  channels = " << m_texture.shape()[3] << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

protected:
    /**
     * \brief Returns the number of channels in the grid
     *
     * For object instances that perform spectral upsampling, the channel that
     * holds all scaling coefficients is omitted.
     */
    MI_INLINE size_t nchannels() const {
        const size_t channels = m_texture.shape()[3];
        // When spectral upsampling is requested, a fourth channel is added to
        // the internal texture data to handle scaling coefficients.
        if (is_spectral_v<Spectrum> && channels == 4 && !m_raw)
            return 3;

        return channels;
    }

    /**
     * \brief Evaluates the volume at the given interaction using spectral
     * upsampling
     */
    MI_INLINE UnpolarizedSpectrum interpolate_spectral(const Interaction3f &it,
                                                        Mask active) const {
        MI_MASK_ARGUMENT(active);

        Point3f p = m_to_local * it.p;

        if (m_texture.filter_mode() == dr::FilterMode::Linear) {
            dr::Array<Float, 4> d000, d100, d010, d110, d001, d101, d011, d111;
            dr::Array<Float *, 8> fetch_values;
            fetch_values[0] = d000.data();
            fetch_values[1] = d100.data();
            fetch_values[2] = d010.data();
            fetch_values[3] = d110.data();
            fetch_values[4] = d001.data();
            fetch_values[5] = d101.data();
            fetch_values[6] = d011.data();
            fetch_values[7] = d111.data();

            if (m_accel)
                m_texture.eval_fetch(p, fetch_values, active);
            else
                m_texture.eval_fetch_nonaccel(p, fetch_values, active);

            UnpolarizedSpectrum v000, v001, v010, v011, v100, v101, v110, v111;
            v000 = srgb_model_eval<UnpolarizedSpectrum>(dr::head<3>(d000), it.wavelengths);
            v100 = srgb_model_eval<UnpolarizedSpectrum>(dr::head<3>(d100), it.wavelengths);
            v010 = srgb_model_eval<UnpolarizedSpectrum>(dr::head<3>(d010), it.wavelengths);
            v110 = srgb_model_eval<UnpolarizedSpectrum>(dr::head<3>(d110), it.wavelengths);
            v001 = srgb_model_eval<UnpolarizedSpectrum>(dr::head<3>(d001), it.wavelengths);
            v101 = srgb_model_eval<UnpolarizedSpectrum>(dr::head<3>(d101), it.wavelengths);
            v011 = srgb_model_eval<UnpolarizedSpectrum>(dr::head<3>(d011), it.wavelengths);
            v111 = srgb_model_eval<UnpolarizedSpectrum>(dr::head<3>(d111), it.wavelengths);

            ScalarVector3i res = resolution();
            p = dr::fmadd(p, res, -.5f);
            Vector3i p_i = dr::floor2int<Vector3i>(p);

            // Interpolation weights
            Point3f w1 = p - Point3f(p_i),
                    w0 = 1.f - w1;

            Float f00 = dr::fmadd(w0.x(), d000.w(), w1.x() * d100.w()),
                  f01 = dr::fmadd(w0.x(), d001.w(), w1.x() * d101.w()),
                  f10 = dr::fmadd(w0.x(), d010.w(), w1.x() * d110.w()),
                  f11 = dr::fmadd(w0.x(), d011.w(), w1.x() * d111.w());
            Float f0 = dr::fmadd(w0.y(), f00, w1.y() * f10),
                  f1 = dr::fmadd(w0.y(), f01, w1.y() * f11);
            Float scale = dr::fmadd(w0.z(), f0, w1.z() * f1);

            UnpolarizedSpectrum v00 = dr::fmadd(w0.x(), v000, w1.x() * v100),
                                v01 = dr::fmadd(w0.x(), v001, w1.x() * v101),
                                v10 = dr::fmadd(w0.x(), v010, w1.x() * v110),
                                v11 = dr::fmadd(w0.x(), v011, w1.x() * v111);
            UnpolarizedSpectrum v0 = dr::fmadd(w0.y(), v00, w1.y() * v10),
                                v1 = dr::fmadd(w0.y(), v01, w1.y() * v11);
            UnpolarizedSpectrum result = dr::fmadd(w0.z(), v0, w1.z() * v1);

            result *= scale;

            return result;
        } else {
            dr::Array<Float, 4> v;
            if (m_accel)
                m_texture.eval(p, v.data(), active);
            else
                m_texture.eval_nonaccel(p, v.data(), active);

            return v.w() * srgb_model_eval<UnpolarizedSpectrum>(dr::head<3>(v), it.wavelengths);
        }
    }

    TensorXf local_majorants(ScalarVector3i resolution_factor,
                             ScalarFloat value_scale) const override {
        using Value   = mitsuba::DynamicBuffer<Float>;
        using Index   = mitsuba::DynamicBuffer<UInt32>;
        using Index3i = mitsuba::Vector<mitsuba::DynamicBuffer<Int32>, 3>;

        if (m_accel)
            NotImplementedError("local_majorants() with m_accel");

        if (m_texture.shape()[3] != 1)
            NotImplementedError("local_majorants() when Channels != 1");

        if constexpr (dr::is_jit_v<Float>) {
            dr::eval(m_texture.value());
            dr::sync_thread();
        }

        // This is the real (user-facing) resolution, but recall
        // that the layout in memory is (Z, Y, X, C).
        const ScalarVector3i full_resolution = resolution();
        if ((full_resolution.x() % resolution_factor.x()) != 0 ||
            (full_resolution.y() % resolution_factor.y()) != 0 ||
            (full_resolution.z() % resolution_factor.z()) != 0) {
            Throw("Supergrid construction: grid resolution %s must be "
                  "divisible by %s",
                  full_resolution, resolution_factor);
        }
        const ScalarVector3i resolution(full_resolution / resolution_factor);

        Log(Debug,
            "Constructing supergrid of resolution %s from full grid of "
            "resolution %s",
            resolution, full_resolution);
        size_t n     = dr::prod(resolution);
        Value result = dr::full<Value>(-dr::Infinity<Float>, n);

        // Z is the slowest axis, X is the fastest.
        auto [Z, Y, X] =
            dr::meshgrid(dr::arange<Index>(resolution.z()),
                         dr::arange<Index>(resolution.y()),
                         dr::arange<Index>(resolution.x()), /*index_xy*/ false);
        Index3i cell_indices = resolution_factor * Index3i(X, Y, Z);

        // We have to include all values that participate in interpolated
        // lookups if we want the true maximum over a region.
        ScalarVector3i begin_offset, end_offset;
        switch (m_texture.filter_mode()) {
            case dr::FilterMode::Nearest: {
                begin_offset = { 0, 0, 0 };
                end_offset   = resolution_factor;
                break;
            }
            case dr::FilterMode::Linear: {
                begin_offset = { -1, -1, -1 };
                end_offset   = resolution_factor + 1;
                break;
            }
        };

        Value scaled_data = value_scale * dr::detach(m_texture.value());

        // TODO: any way to do this without the many operations?
        for (int32_t dz = begin_offset.z(); dz < end_offset.z(); ++dz) {
            for (int32_t dy = begin_offset.y(); dy < end_offset.y(); ++dy) {
                for (int32_t dx = begin_offset.x(); dx < end_offset.x(); ++dx) {
                    // Ensure valid lookups with clamping
                    Index3i offset_indices =
                        Index3i(dr::clamp(cell_indices.x() + dx, 0,
                                          full_resolution.x() - 1),
                                dr::clamp(cell_indices.y() + dy, 0,
                                          full_resolution.y() - 1),
                                dr::clamp(cell_indices.z() + dz, 0,
                                          full_resolution.z() - 1));
                    // Linearize indices
                    const Index idx = offset_indices.x() +
                                      offset_indices.y() * full_resolution.x() +
                                      offset_indices.z() * full_resolution.x() *
                                          full_resolution.y();

                    Value values = dr::gather<Value>(scaled_data, idx);
                    result       = dr::maximum(result, values);
                }
            }
        }

        size_t shape[4] = { (size_t) resolution.z(), (size_t) resolution.y(),
                            (size_t) resolution.x(), 1 };
        return TensorXf(result, 4, shape);
    }

    ref<Volume<Float, Spectrum>>
    get_majorant_grid(ScalarVector3i resolution_factor,
                      ScalarFloat value_scale) const override {
        // Compute majorants
        TensorXf majorants = local_majorants(resolution_factor, value_scale);
        dr::eval(majorants);
        const auto &shape = majorants.shape();

        // Initialize a volume texture with these
        ref<VolumeGrid> majorant_grid = new VolumeGrid(
            ScalarVector3u(shape[2], shape[1], shape[0]),
            (uint32_t) shape[3]
        );
        memcpy(majorant_grid->data(), majorants.data(), majorant_grid->buffer_size());

        // These are technically useless, but we set them anyway for consistency
        majorant_grid->set_max(m_max);
        majorant_grid->set_max_per_channel((ScalarFloat*) &m_max_per_channel[0]);

        // Initialize a gridvolume plugin with those data
        Properties props("gridvolume");
        props.set_string("filter_type", "nearest");
        props.set_string("wrap_mode", "clamp");
        props.set_transform("to_world", m_to_local.inverse());
        props.set_object("grid", majorant_grid.get());
        ref<Volume<Float, Spectrum>> result =
            PluginManager::instance()->create_object<Volume<Float, Spectrum>>(props);

        return result.get();
    }

    std::tuple<Float, Vector3f, Vector3f>
    prepare_majorant_grid_traversal(const Ray3f &ray, Float mint, Float maxt,
                                    Mask /*active*/) const override {
        const auto extents = m_bbox.extents();
        Ray3f local_ray(
            /* o */ (ray.o - m_bbox.min) / extents,
            /* d */ ray.d / extents, ray.time, ray.wavelengths);
        const ScalarVector3i res  = resolution();
        Vector3f local_voxel_size = 1.f / res;

        // The ID of the first and last voxels hit by the ray
        Vector3f current_voxel_coords(
            dr::clamp(local_ray(mint) / local_voxel_size, 0.f, res - 1.f)
        );
        Vector3f last_voxel_coords(
            dr::clamp(local_ray(maxt) / local_voxel_size, 0.f, res - 1.f)
        );
        Vector3i current_voxel(dr::floor(current_voxel_coords));
        Vector3i last_voxel(dr::floor(last_voxel_coords));
        Log(Info, "current_voxel = %s", current_voxel);
        Log(Info, "last_voxel = %s", last_voxel);

        // By definition, current and last voxels should be valid voxel indices.
        // current_voxel = dr::clamp(current_voxel, 0, res - 1);
        // last_voxel    = dr::clamp(last_voxel, 0, res - 1);
        // if (dr::any_nested(last_voxel < 0) || dr::any_nested(last_voxel > res - 1))
        //     Log(Error, "invalid voxel index detected: %s", last_voxel);

        // Increment (in number of voxels) to take at each step
        Vector3i step = dr::select(local_ray.d >= 0, 1, -1);

        // Distance along the ray to the next voxel border from the current
        // position
        Vector3f next_voxel_boundary =
            (current_voxel + step) * local_voxel_size;
        next_voxel_boundary +=
            dr::select(dr::neq(current_voxel, last_voxel) && (local_ray.d < 0),
                       local_voxel_size, 0);

        // Value of ray parameter until next intersection with voxel-border
        // along each axis
        auto ray_nonzero  = dr::neq(local_ray.d, 0);
        Vector3f dda_tmax = dr::select(
            ray_nonzero, (next_voxel_boundary - local_ray.o) / local_ray.d,
            dr::Infinity<Float>);

        // How far along each component of the ray we must move to move by one
        // voxel
        Vector3f dda_tdelta =
            dr::select(ray_nonzero, step * local_voxel_size / local_ray.d,
                       dr::Infinity<Float>);

        // Current ray parameter throughout DDA traversal
        Float dda_t = mint;

        // Note: `t` parameters on the reparametrized ray yield locations on the
        // normalized majorant supergrid in [0, 1]^3. But they are also directly
        // valid parameters on the original ray, yielding positions in the
        // bbox-aligned supergrid.
        return { dda_t, dda_tmax, dda_tdelta };
    }

    Float traverse_majorant_grid(
            const Float desired_tau, const Ray3f &ray,
            const Float mint, const Float maxt, Mask active) const override {

        // 1. Prepare for DDA traversal
        // Adapted from:
        // https://github.com/francisengelmann/fast_voxel_traversal/blob/9664f0bde1943e69dbd1942f95efc31901fbbd42/main.cpp
        // TODO: allow precomputing all this (but be careful when ray origin is
        // updated)
        auto [dda_t, dda_tmax, dda_tdelta] =
            prepare_majorant_grid_traversal(ray, mint, maxt, active);
        MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();

        // 2. Traverse the medium with DDA until we reach the desired
        // optical depth
        Mask active_dda = active;
        Mask reached    = false;
        Float tau_acc   = 0.f;
        size_t i = 0;
        dr::Loop<Mask> dda_loop("Volume::traverse_majorant_grid");
        dda_loop.put(active_dda, reached, dda_t, dda_tmax, tau_acc, mei);
        dda_loop.init();
        Log(Info, "Entering DDA loop");
        Log(Info, "active_dda = %s", active_dda);


        while (dda_loop(dr::detach(active_dda))) {
            Log(Info, "---- i = %s ----", i);
            Log(Info, "reached = %s", reached);
            Log(Info, "dda_t = %s", dda_t);
            Log(Info, "dda_tmax = %s", dda_tmax);
            Log(Info, "tau_acc = %s", tau_acc);

            // Figure out which axis we hit first.
            // `t_next` is the ray's `t` parameter when hitting that axis.
            Float t_next = dr::min(dda_tmax);
            Log(Info, "t_next = %s", t_next);
            Vector3f tmax_update;
            Mask got_assigned = false;
            for (size_t k = 0; k < 3; ++k) {
                Mask active_k = dr::eq(dda_tmax[k], t_next);
                tmax_update[k] =
                    dr::select(!got_assigned && active_k, dda_tdelta[k], 0);
                got_assigned |= active_k;
            }

            // Lookup and accumulate majorant in current cell.
            dr::masked(mei.t, active_dda) = .5f * (dda_t + t_next);
            dr::masked(mei.p, active_dda) = ray(mei.t);
            Float majorant = eval_1(mei, active_dda);
            Log(Info, "majorant = %s", majorant);
            Float tau_next = tau_acc + majorant * (t_next - dda_t);
            Log(Info, "tau_next = %s", tau_next);

            // For rays that will stop within this cell, figure out
            // the precise `t` parameter where `desired_tau` is reached.
            Float t_precise = dda_t + (desired_tau - tau_acc) / majorant;
            reached |= active_dda && (majorant > 0) && (t_precise < maxt) &&
                       (tau_next >= desired_tau);
            dr::masked(dda_t, active_dda) =
                dr::select(reached, t_precise, t_next);

            // Prepare for next iteration
            active_dda &= !reached && (t_next < maxt);
            dr::masked(dda_tmax, active_dda) = dda_tmax + tmax_update;
            dr::masked(tau_acc, active_dda)  = tau_next;
            i++;
        }
        // Adopt the stopping location, making sure to convert to the main
        // ray's parametrization.
        Log(Info, "Exiting DDA loop");

        Log(Info, "---- i = %s ----", i);
        Log(Info, "reached = %s", reached);
        Log(Info, "dda_t = %s", dda_t);
        Log(Info, "dda_tmax = %s", dda_tmax);
        Log(Info, "tau_acc = %s", tau_acc);

        return dr::select(reached, dda_t, dr::Infinity<Float>);
    }

    /**
     * \brief Evaluates the volume at the given interaction
     *
     * Should only be used when the volume data has exactly 1 channel.
     */
    MI_INLINE Float interpolate_1(const Interaction3f &it, Mask active) const {
        MI_MASK_ARGUMENT(active);

        Point3f p = m_to_local * it.p;
        Float result;
        if (m_accel)
            m_texture.eval(p, &result, active);
        else
            m_texture.eval_nonaccel(p, &result, active);

        return result;
    }

    /**
     * \brief Evaluates the volume at the given interaction
     *
     * Should only be used when the volume data has exactly 3 channels.
     */
    MI_INLINE Color3f interpolate_3(const Interaction3f &it,
                                     Mask active) const {
        MI_MASK_ARGUMENT(active);

        Point3f p = m_to_local * it.p;
        Color3f result;
        if (m_accel)
            m_texture.eval(p, result.data(), active);
        else
            m_texture.eval_nonaccel(p, result.data(), active);

        return result;
    }

    /**
     * \brief Evaluates the volume at the given interaction
     *
     * Should only be used when the volume data has exactly 6 channels.
     */
    MI_INLINE dr::Array<Float, 6> interpolate_6(const Interaction3f &it,
                                                 Mask active) const {
        MI_MASK_ARGUMENT(active);

        Point3f p = m_to_local * it.p;
        dr::Array<Float, 6> result;
        if (m_accel)
            m_texture.eval(p, result.data(), active);
        else
            m_texture.eval_nonaccel(p, result.data(), active);

        return result;
    }


    /**
     * \brief Evaluates the volume at the given interaction
     *
     * Should only be used when the volume contains data with multiple channels
     */
    template <typename VALUE>
    MI_INLINE void interpolate_per_channel(const Interaction3f &it, Float *out, Mask active) const {
        MI_MASK_ARGUMENT(active);

        Point3f p = m_to_local * it.p;
        if (m_accel)
            m_texture.eval(p, out, active);
        else
            m_texture.eval_nonaccel(p, out, active);
    }

protected:
    Texture3f m_texture;
    bool m_accel;
    bool m_raw;
    bool m_fixed_max = false;
    ScalarFloat m_max;
    std::vector<ScalarFloat> m_max_per_channel;
};

MI_IMPLEMENT_CLASS_VARIANT(GridVolume, Volume)
MI_EXPORT_PLUGIN(GridVolume, "GridVolume texture")

NAMESPACE_END(mitsuba)
