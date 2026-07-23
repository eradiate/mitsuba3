set(MI_PLUGIN_PREFIX "shapes")

add_plugin(arectangle arectangle.cpp)
add_plugin(instancelist instancelist.cpp)
add_plugin(instancelist_legacy instancelist_legacy.cpp)

set(MI_PLUGIN_TARGETS "${MI_PLUGIN_TARGETS}" PARENT_SCOPE)
