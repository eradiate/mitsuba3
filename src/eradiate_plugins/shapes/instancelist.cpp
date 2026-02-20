#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/shapegroup.h>
#include <mitsuba/core/transform.h>
#include <drjit/tensor.h>
#include <format>

NAMESPACE_BEGIN(mitsuba)
/**!

.. _shape-instancelist:

InstanceList (:monosp:`instancelist`)
-------------------------------------------------

.. pluginparameters::

 * - (Nested plugin)
   - :paramtype:`shapegroup`
   - A reference to a shape group that should be instantiated.

 * - transfroms
   - |tensor|
   - Specifies a buffer of object-to-world transformation of shape [N,4,4] with 
     N the being the number of instances. 
   - |exposed|

This is a utility plugins that allows to create many instances of a Shapegroup
from a buffer of object-to-world transforms. For details on how to create 
instances, refer to the :ref:`shape-shapegroup` and :ref: `instances` plugin.
 */

template <typename Float, typename Spectrum>
class InstanceList final: public Shape<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Shape, id)
    MI_IMPORT_TYPES()

    using ShapeGroup_ = ShapeGroup<Float, Spectrum>;

    InstanceList(const Properties &props) : Base(props) {
        Log(Debug, "Start Constructor");

        m_shapegroup = nullptr; 
        for (auto &prop : props.objects()) {
            if (Base *shape = prop.try_get<Base>()) {
                if (shape->is_shape_group()) {
                    if (m_shapegroup)
                        Throw("Only a single shapegroup can be specified per instance list.");
                    m_shapegroup = (ShapeGroup_ *) shape;
                }
            }
        }

        if (!m_shapegroup)
            Throw("A reference to a 'shapegroup' must be specified!");

        m_tensor = const_cast<TensorXf*>(&props.get_any<TensorXf>("transforms"));
        if (m_tensor->ndim() != 3)
            Throw("Tensor->has %ul dimensions. Expected 3", m_tensor->ndim());

        m_instance_count = (uint32_t) m_tensor->shape(0);
    }

    std::vector<ref<Object>> expand() const override {
        std::vector<ref<Object>> shapes;
        auto pmgr = PluginManager::instance();
        
        std::string id_prefix = std::string(id());
        for (uint32_t i = 0; i < m_instance_count; ++i) {
            ScalarAffineTransform4f to_world( dr::gather<ScalarMatrix4f>(m_tensor->array(), i) );
            Properties props_instance("instance");
            props_instance.set("shapegroup", (Object *)m_shapegroup);
            props_instance.set("to_world", to_world);

            
            ref<Object> shape = (Object *)pmgr->create_object<Base>(props_instance);
            
            if(id_prefix != ""){
                std::string id = id_prefix + "_" + std::to_string(i);
                shape->set_id(id);
            }

            shapes.push_back(shape);
        }
        return shapes;
    }

    ScalarBoundingBox3f bbox() const override {
        return dr::zeros<ScalarBoundingBox3f>();
    }

    MI_DECLARE_CLASS()
private:
    uint32_t m_instance_count;
    TensorXf *m_tensor;
    ShapeGroup_ *m_shapegroup;
};

MI_EXPORT_PLUGIN(InstanceList)
NAMESPACE_END(mitsuba)