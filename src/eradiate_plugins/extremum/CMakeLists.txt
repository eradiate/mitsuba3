set(MI_PLUGIN_PREFIX "extremumstructure")

# #ERADIATE_CHANGE: Add extremum_grid plugin for local majorant acceleration
add_plugin(extremum_grid extremum_grid.cpp)
add_plugin(extremum_spherical extremum_spherical.cpp)

set(MI_PLUGIN_TARGETS "${MI_PLUGIN_TARGETS}" PARENT_SCOPE)
