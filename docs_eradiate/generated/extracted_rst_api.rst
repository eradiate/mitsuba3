.. py:data:: mitsuba.ERD_MI_VERSION
    :type: str
    :value: 0.4.3

.. py:data:: mitsuba.ERD_MI_VERSION_MAJOR
    :type: int
    :value: 0

.. py:data:: mitsuba.ERD_MI_VERSION_MINOR
    :type: int
    :value: 4

.. py:data:: mitsuba.ERD_MI_VERSION_PATCH
    :type: int
    :value: 3

.. py:class:: mitsuba.ExtremumSegment

    Stores the extremum (minorant/majorant) data for a ray segment.

    Used as the output type of ExtremumStructure traversal. Tracks the
    segment's entry/exit distances and the local extinction coefficient
    bounds within that interval.


    .. py:method:: __init__()


    .. py:method:: __init__(self, other)

        Copy constructor

        Parameter ``other`` (:py:obj:`mitsuba.ExtremumSegment`):
            *no description available*

    .. py:method:: __init__(self, mint, maxt, minorant, majorant)

        Construct from entry/exit distances and a combined extremum vector.

        Parameter ``mint`` (float):
            Segment entry distance

        Parameter ``maxt`` (float):
            Segment exit distance

        Parameter ``value`` (:py:obj:`mitsuba.Vector2f`):
            Extremum vector [minorant, majorant]

        Parameter ``minorant`` (float):
            *no description available*

        Parameter ``majorant`` (float):
            *no description available*

    .. py:method:: __init__(self, mint, maxt, value)

        Construct from entry/exit distances and a combined extremum vector.

        Parameter ``mint`` (float):
            Segment entry distance

        Parameter ``maxt`` (float):
            Segment exit distance

        Parameter ``value`` (:py:obj:`mitsuba.Vector2f`):
            Extremum vector [minorant, majorant]

    .. py:method:: mitsuba.ExtremumSegment.assign(self, arg)

        Parameter ``arg`` (:py:obj:`mitsuba.ExtremumSegment`, /):
            *no description available*

        Returns → None:
            *no description available*

    .. py:method:: mitsuba.ExtremumSegment.majorant()

        Majorant value over the segment. Accessor to the second element of
        ``value``.

        Returns → float:
            *no description available*

    .. py:property:: mitsuba.ExtremumSegment.maxt

        Segment exit distance along ray

    .. py:method:: mitsuba.ExtremumSegment.minorant()

        Minorant value over the segment. Accessor to the first element of
        ``value``.

        Returns → float:
            *no description available*

    .. py:property:: mitsuba.ExtremumSegment.mint

        Segment entry distance along ray

    .. py:method:: mitsuba.ExtremumSegment.reset()

        Mark the extremum segment as invalid.

        This operation sets segment's minimum and maximum distances to
        :math:`\infty` and :math:`-\infty`, respectively.

        Returns → None:
            *no description available*

    .. py:method:: mitsuba.ExtremumSegment.valid()

        Check whether this is a valid segment

        A segment is considered valid when

        .. code-block:: c

            segment.mint < segment.maxt


        Returns → bool:
            *no description available*

    .. py:property:: mitsuba.ExtremumSegment.value

        Extremum data stored as [minorant, majorant]

    .. py:method:: mitsuba.ExtremumSegment.zero_(self, size=1)

        Overloaded function.

        1. ``zero_(self, size: int = 1) -> None``


        2. ``zero_(self, arg: int, /) -> None``

        This callback method is invoked by dr::zeros<>, and takes care of
        fields that deviate from the standard zero-initialization convention.
        In ExtremumSegment, the ``mint`` and ``maxt`` fields are set to + and
        - infinity respectively to to mark invalid intersection records.

        Parameter ``size`` (int):
            *no description available*

        Returns → None:
            *no description available*

.. py:class:: mitsuba.ExtremumStructure

    Base class: :py:obj:`mitsuba.Object`

    Abstract base class for extremum structures

    ExtremumStructure provides an interface for spatial data structures
    that store local extrema (majorant/minorant) of volumetric extinction
    coefficients. This enables efficient delta tracking with locally-
    adaptive majorants.

    To minimize virtual function overhead, the ``traverse_extremum()``
    method encapsulates the entire traversal loop internally, requiring
    only a single virtual call per distance sample.

    .. py:method:: __init__(self, props)

        Parameter ``props`` (:py:obj:`mitsuba.Properties`):
            *no description available*


    .. py:method:: mitsuba.ExtremumStructure.bbox()

        Return the bounding box of the extremum structure

        Returns → :py:obj:`mitsuba.BoundingBox3f`:
            *no description available*

    .. py:method:: mitsuba.ExtremumStructure.eval_1(self, it, active=True)

        Evaluate the minorant and majorant at a medium interaction point.

        This method performs point evaluation at interaction point specified
        in local space.

        Parameter ``it`` (:py:obj:`mitsuba.Interaction3f`):
            Interaction interaction point in local space

        Parameter ``active`` (bool):
            Mask for active lanes

        Returns → tuple[float, float]:
            The minorant and majorant values at the medium interaction point.
            Clamped values outside bounds.

.. py:class:: mitsuba.Medium

    Base class: :py:obj:`mitsuba.Object`

    .. py:method:: mitsuba.Medium.get_majorant(self, mi, active=True)

        Returns the medium's majorant used for delta tracking

        Parameter ``mi`` (:py:obj:`mitsuba.MediumInteraction3f`):
            *no description available*

        Parameter ``active`` (bool):
            Mask to specify active lanes.

        Returns → :py:obj:`mitsuba.Color3f`:
            *no description available*

    .. py:method:: mitsuba.Medium.get_scattering_coefficients(self, mi, active=True)

        Returns the medium coefficients Sigma_s, Sigma_n and Sigma_t evaluated
        at a given MediumInteraction mi

        Parameter ``mi`` (:py:obj:`mitsuba.MediumInteraction3f`):
            *no description available*

        Parameter ``active`` (bool):
            Mask to specify active lanes.

        Returns → tuple[:py:obj:`mitsuba.Color3f`, :py:obj:`mitsuba.Color3f`, :py:obj:`mitsuba.Color3f`]:
            *no description available*

    .. py:method:: mitsuba.Medium.has_spectral_extinction()

        Returns whether this medium has a spectrally varying extinction

        Returns → bool:
            *no description available*

    .. py:method:: mitsuba.Medium.intersect_aabb(self, ray)

        Intersects a ray with the medium's bounding box

        Parameter ``ray`` (:py:obj:`mitsuba.Ray3f`):
            *no description available*

        Returns → tuple[bool, float, float]:
            *no description available*

    .. py:method:: mitsuba.Medium.is_homogeneous()

        Returns whether this medium is homogeneous

        Returns → bool:
            *no description available*

    .. py:property:: mitsuba.Medium.m_has_spectral_extinction

        (self) -> bool

    .. py:property:: mitsuba.Medium.m_is_homogeneous

        (self) -> bool

    .. py:property:: mitsuba.Medium.m_sample_emitters

        (self) -> bool

    .. py:method:: mitsuba.Medium.phase_function()

        Return the phase function of this medium

        Returns → :py:obj:`mitsuba.PhaseFunction`:
            *no description available*

    .. py:method:: mitsuba.Medium.sample_interaction(self, ray, sample, channel, active)

        Sample a free-flight distance in the medium.

        This function samples a (tentative) free-flight distance according to
        an exponential transmittance. It is then up to the integrator to then
        decide whether the MediumInteraction corresponds to a real or null
        scattering event.

        Parameter ``ray`` (:py:obj:`mitsuba.Ray3f`):
            Ray, along which a distance should be sampled

        Parameter ``sample`` (float):
            A uniformly distributed random sample

        Parameter ``channel`` (int):
            The channel according to which we will sample the free-flight
            distance. This argument is only used when rendering in RGB modes.

        Parameter ``active`` (bool):
            Mask to specify active lanes.

        Returns → :py:obj:`mitsuba.MediumInteraction3f`:
            This method returns a MediumInteraction. The MediumInteraction
            will always be valid, except if the ray missed the Medium's
            bounding box.

    .. py:method:: mitsuba.Medium.sample_interaction_analytical(self, ray, it, sample, channel, active)

        Sample a free-flight distance in the medium analytically.

        This function samples a (tentative) free-flight distance according to
        an exponential transmittance. It is then up to the integrator to then
        decide whether the MediumInteraction corresponds to a real or null
        scattering event.

        Parameter ``ray`` (:py:obj:`mitsuba.Ray3f`):
            Ray, along which a distance should be sampled

        Parameter ``it`` (:py:obj:`mitsuba.Interaction3f`):
            The boundary interaction that the sampled distance cannot exceed.

        Parameter ``sample`` (float):
            A uniformly distributed random sample

        Parameter ``channel`` (int):
            The channel according to which we will sample the free-flight
            distance. This argument is only used when rendering in RGB modes.

        Parameter ``active`` (bool):
            Mask to specify active lanes.

        Returns → tuple[:py:obj:`mitsuba.MediumInteraction3f`, :py:obj:`mitsuba.Color3f`, :py:obj:`mitsuba.Color3f`]:
            This method returns a MediumInteraction. The MediumInteraction
            will always be valid, except if the ray missed the Medium's
            bounding box.

    .. py:method:: mitsuba.Medium.transmittance_eval_analytical(self, ray, it, active)

        Compute the analytical transmittance along a ray to an interaction.

        Parameter ``ray`` (:py:obj:`mitsuba.Ray3f`):
            Ray, along which to compute the transmittance, use mint

        Parameter ``si``:
            Interaction that marks the end of the segment along which to
            compute the transmittance.

        Parameter ``it`` (:py:obj:`mitsuba.Interaction3f`):
            *no description available*

        Parameter ``active`` (bool):
            Mask to specify active lanes.

        Returns → :py:obj:`mitsuba.Color3f`:
            The transmittance along a ray

    .. py:method:: mitsuba.Medium.transmittance_eval_pdf(self, mi, si, active)

        Compute the transmittance and PDF

        This function evaluates the transmittance and PDF of sampling a
        certain free-flight distance The returned PDF takes into account if a
        medium interaction occurred (mi.t <= si.t) or the ray left the medium
        (mi.t > si.t)

        The evaluated PDF is spectrally varying. This allows to account for
        the fact that the free-flight distance sampling distribution can
        depend on the wavelength.

        Parameter ``mi`` (:py:obj:`mitsuba.MediumInteraction3f`):
            *no description available*

        Parameter ``si`` (:py:obj:`mitsuba.SurfaceInteraction3f`):
            *no description available*

        Parameter ``active`` (bool):
            Mask to specify active lanes.

        Returns → tuple[:py:obj:`mitsuba.Color3f`, :py:obj:`mitsuba.Color3f`]:
            This method returns a pair of (Transmittance, PDF).

    .. py:method:: mitsuba.Medium.use_emitter_sampling()

        Returns whether this specific medium instance uses emitter sampling

        Returns → bool:
            *no description available*

