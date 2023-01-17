pkgname <- "ROptimus"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ROptimus')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Optimus")
### * Optimus

flush(stderr()); flush(stdout())

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



cleanEx()
nameEx("OptimusExamples")
### * OptimusExamples

flush(stderr()); flush(stdout())

### Name: OptimusExamples
### Title: Generate script for reproducing a tutorial
### Aliases: OptimusExamples

### ** Examples


out.dir <- tempdir()
OptimusExamples(dir=out.dir, example=1)



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
