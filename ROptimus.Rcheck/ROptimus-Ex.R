pkgname <- "ROptimus"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "ROptimus-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('ROptimus')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Optimus")
### * Optimus

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Optimus
### Title: Acceptance Ratio Simulated Annealing and Acceptance Ratio
###   Replica Exchange Monte Carlo Optimisation Engine
### Aliases: Optimus

### ** Examples


K <- IJ_ORIG
K$j <- sample(x=K$j, size=nrow(K), replace=FALSE)
out.dir <- tempdir()

Optimus(NCPU=1, OPTNAME="IJ.NEW.OPTI.SA", NUMITER=500, CYCLES=2, DIR=out.dir,
        DUMP.FREQ=10, LONG=FALSE, OPT.TYPE="SA", K.INITIAL=K,
        rDEF=ex.r.fun, mDEF=ex.m.fun, uDEF=ex.u.fun,
        DATA=list(IJ_ORIG=IJ_ORIG, gaplimit=50, numContacts=nrow(IJ_ORIG)))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Optimus", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("OptimusExamples")
### * OptimusExamples

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: OptimusExamples
### Title: Generate script for reproducing a tutorial
### Aliases: OptimusExamples

### ** Examples


out.dir <- tempdir()
OptimusExamples(dir=out.dir, example=1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("OptimusExamples", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
