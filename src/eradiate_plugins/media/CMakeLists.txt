set(MI_PLUGIN_PREFIX "media")

add_plugin(piecewise piecewise.cpp)
add_plugin(piecewise_spherical piecewise_spherical.cpp)
add_plugin(piecewise_spherical_2 piecewise_spherical_2.cpp)

set(MI_PLUGIN_TARGETS "${MI_PLUGIN_TARGETS}" PARENT_SCOPE)
