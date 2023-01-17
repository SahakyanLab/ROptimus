## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## 15-Dec-2022 First submission, new release

Comments received from CRAN and subsequent modifications:

1. Please do not start the description with "This package", package name, title or similar.

--> Modified to A general-purpose optimisation engine that supports...

2. You write information messages to the console that cannot be easily suppressed. It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object.Instead of print()/cat() rather use message()/warning() or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console.(except for print, summary, interactive functions)

--> Changed all instances of print() to message() and one print() line in Optimus() to stop()

3. Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()). This is not allowed by CRAN policies.Please omit any default path in writing functions. In your examples/vignettes/tests you can write to tempdir().

--> Removed the default value (i.e. current working directory) for the DIR parameter in all affected functions to require users to supply a path 

4. Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited. e.g.:

...

oldpar <- par(no.readonly = TRUE) # code line i
on.exit(par(oldpar)) # code line i + 1
...
par(mfrow=c(5,1)) # somewhere after
...
e.g.: R/OptimusRE.R ; R/OptimusSA.R

If you're not familiar with the function, please check ?on.exit. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.

--> Added on.exit code as suggested

5. Please do not set a seed to a specific number within a function.

--> This is because of the OptimusExamples() function, which produces an R script that the user can run separately at their discretion if they want to replicate any of the examples discussed in the vignette. All the R scripts that can be produced using the function contain the line:

...
set.seed(845) (Examples 1 to 4) or set.seed(seed) where seed is 840 (Example 5)
...

To make sure that this line will never be applied every time OptimusExamples() is called, the run parameter, which
if TRUE runs the script after generating, was removed. 

