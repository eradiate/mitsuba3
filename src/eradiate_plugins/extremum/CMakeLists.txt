set(MI_PLUGIN_PREFIX "extremumstructure")

add_plugin(extremum_grid extremum_grid.cpp)
add_plugin(extremum_irregular extremum_irregular.cpp)
add_plugin(extremum_irregular2 extremum_irregular2.cpp)
add_plugin(extremum_spherical extremum_spherical.cpp)
add_plugin(extremum_global extremum_global.cpp)

set(MI_PLUGIN_TARGETS "${MI_PLUGIN_TARGETS}" PARENT_SCOPE)
