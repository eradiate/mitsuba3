set(MI_PLUGIN_PREFIX "extremumstructure")

add_plugin(extremum_grid extremum_grid.cpp)
add_plugin(extremum_spherical extremum_spherical.cpp)
add_plugin(extremum_global extremum_global.cpp)
add_plugin(extremum_overlap extremum_overlap.cpp)
add_plugin(extremum_repeat extremum_repeat.cpp)

set(MI_PLUGIN_TARGETS "${MI_PLUGIN_TARGETS}" PARENT_SCOPE)
