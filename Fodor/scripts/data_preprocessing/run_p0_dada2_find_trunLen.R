#!/usr/bin/env Rscript
# Running the p0_dada2_find_trunLen.R with command line args
# Note: download result with scp amy@hpc.uncc.edu:~/git/balance_tree_exploration/McDonald/output/graphics/plotQualF.png .


rm(list = ls()) #clear workspace

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Fodor"

r_script <- file.path(home_dir, "r_libraries", "cml_scripts", "data_preprocessing","p0_dada2_find_trunLen.R")

my_args <- paste(
  "-d", home_dir,
  "-p", project,
  "-f", "_1.fastq.gz",
  "-r", "_2.fastq.gz"
)

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
    print('Opps, an error is thrown')
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

