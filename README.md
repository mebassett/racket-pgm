# racket-pgm

Work in Progress.  This is not a license.

An ad-hoc, informally specified, bug-ridden, slow implementation of half of 
python's libpgm.

The author felt that it would be instructive (to him) to implement 
representation, inference and learning for Bayesian networks, and that it would 
be fun to do so in typed racket.

At present this implements variable-elimination inference (for conditional 
probability queries P( Y | E = e)) for a Bayesian network with discrete 
variables only.

The author would this project/package to be able to do:

* continuous random variables - id est, a linear Gaussian network. 
* VE with canonical factors
* conditional linear Gaussian networks.
* related VE (is that even possible?)
* MLE parameter learning over complete data.
* EM-MLE parameter learning over incomplete data.

The author aware that VE for conditional linear Gaussian would require that the
space of intermediate factors be tables whose values are mixtures of canonical 
factors, and that this would have an exponential space requirement.  The author, 
however, has access to machines with a terabyte of ram and does not care.  At the 
worst, this project could handle fun-sized networks without choking.

The author would like to implement other inference algorithms besides VE, as he
feels they would be equally instructive.  Alas, the author does not currently
understand many of them, while others are just frightening.  These including:
* Message Passing
* believe propagation
* forward sampling
* MCMC
* variational inference

Additionally, the author would like to implement Bayesian learning methods, but
he is afraid of them especially in the case of missing or incomplete data.

If the reader scoffs at the lack of intellectual intrepidness of the author, then the
reader is of course invited to contribute. :)

Example datasets: 


    Lichman, M. (2013). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. Irvine, CA: University of California, School of Information and Computer Science.

