# Creates full-range overview plots (collection assumed open).

set_var COLLDIR  = /home/brandon/vestawd/omz/wd/r_wd/Profiles
set_var COLLNAME = Argo_Profile_Collection_2022-03-11T15-20-13
set_var VIEWSDIR = %COLLDIR%/%COLLNAME%.Data/views/_overviews

load_view %VIEWSDIR%/OVERVIEW_AllStationsMap.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_AllStationsMap.gif, 300

load_view %VIEWSDIR%/OVERVIEW_vars_0002-0011.xview
set_axis_full_ranges -1
set_axis_full_ranges 0
save_view %VIEWSDIR%/OVERVIEW_vars_0002-0011.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_vars_0002-0011.gif, 300

