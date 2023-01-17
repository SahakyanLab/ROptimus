# ROptimus
ROptimus: a parallel general purpose adaptive optimisation engine

A general-purpose optimisation engine that supports:

i) Monte Carlo optimisation with Metropolis criterion [Metropolis et al. (1953) <doi:10.1063/1.1699114>, Hastings (1970) <doi:10.1093/biomet/57.1.97>] and Acceptance Ratio Simulated Annealing [Kirkpatrick et al. (1983) <doi:10.1126/science.220.4598.671>, Černý (1985) <doi:10.1007/BF00940812>] 
	on multiple cores, and 
  
ii) Acceptance Ratio Replica Exchange Monte Carlo Optimisation. In each case, the system pseudo-temperature is dynamically adjusted such that the observed acceptance ratio is kept near to the desired (fixed or changing) acceptance ratio.

Latest stable version is available in CRAN https://cran.r-project.org/web/packages/ROptimus/index.html
