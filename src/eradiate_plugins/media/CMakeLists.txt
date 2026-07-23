set(MI_PLUGIN_PREFIX "media")

add_plugin(piecewise piecewise.cpp)
add_plugin(eoheterogeneous eoheterogeneous.cpp)
add_plugin(multicomponent multicomponent.cpp)
add_plugin(repeat repeat.cpp)

set(MI_PLUGIN_TARGETS "${MI_PLUGIN_TARGETS}" PARENT_SCOPE)
