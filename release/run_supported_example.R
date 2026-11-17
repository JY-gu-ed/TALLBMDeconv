# Top-level backend initialisation avoids the documented nested lazy-load warning.
# This does not suppress warnings or change any inference parameter.
.libPaths(c(Sys.getenv('TALLBMDECONV_LIBRARY','/work/R_library'),.libPaths()))
loadNamespace('BayesPrism')
library(TALLBMDeconv)
Sys.setenv(TZ='UTC')
args<-commandArgs(TRUE)
destination<-if(length(args))args[1]else'new_synthetic_run'
if(dir.exists(destination))stop('Use a NEW example output directory')
dir.create(destination,recursive=TRUE)
entry<-normalizePath('PUBLIC_EXAMPLE.R');setwd(destination);source(entry)
cat('PASS_SUPPORTED_STARTUP_AND_SYNTHETIC_WORKFLOW\n')
