set(MI_PLUGIN_PREFIX "media")

add_plugin(piecewise piecewise.cpp)
add_plugin(eoheterogeneous eoheterogeneous.cpp)
add_plugin(multicomponent multicomponent.cpp)

set(MI_PLUGIN_TARGETS "${MI_PLUGIN_TARGETS}" PARENT_SCOPE)
