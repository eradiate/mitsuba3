set(MI_PLUGIN_PREFIX "integrators")

add_plugin(eovolpath eovolpath.cpp)
add_plugin(stokes_moment stokes_moment.cpp)

set(MI_PLUGIN_TARGETS "${MI_PLUGIN_TARGETS}" PARENT_SCOPE)
