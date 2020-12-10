# Reset the rebuild parameters
##
#
#
#
## ReBuild settings
rebuild_time = 10
reference_time = 100
management_state = "none" # "none", "rebuild", "closed"
start_rebuild_year = 0 # if in rebuild phase this is how many years are left in the rebuild phase.
n_proj_sims = 100# projections into the future used for both B_ref and rebuild F.
catch_tonnes_to_search = 250 # consider future catches in 250 increments
is_final_year = 0 ## when this gets triggered == 1, then we are in final year and don't have to bother about projections