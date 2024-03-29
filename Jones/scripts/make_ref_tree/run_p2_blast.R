#!/usr/bin/env Rscript
# Running the p2_blast.sh with command line args
# if this doesn't work, you might need to run the scripts 
# in balance_tree_exploration/r_libraries/cml_scripts/creat_ref_tree_blst_db

rm(list = ls()) #clear workspace

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Jones"
cml_scripts <- file.path(home_dir, "r_libraries", "cml_scripts")
r_script <- file.path(cml_scripts, "make_ref_tree", "p2_blast.sh")

##-Make args for cml script-----------------------------------------##
my_args <- paste( home_dir, project)

sys_command <- paste(r_script, my_args)

tryCatch(
  { 
  system(sys_command,
         intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL,
         minimized = FALSE, invisible = TRUE, timeout = 0)
  },
  error=function(cond) {
    print('Opps, an error is thrown, did you run balance_tree_exploration/r_libraries/cml_scripts/creat_ref_tree_blst_db?')
    message(cond)
  },
  warning=function(cond) {
    print('Oppsa warning is thrown')
    message(cond)
    # Choose a return value in case of warning
    #return(NULL)
  }
)
print("Subscript complete!")

