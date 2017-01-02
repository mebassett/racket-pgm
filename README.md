# racket-pgm

Work in Progress.  This is not a license.

An ad-hoc, informally specified, bug-ridden, slow implementation of half of python's libpgm.

The author felt that it would be instructive (to him) to implement representation, inference and learning for a pgm, and that it would be fun to do so in typed racket.

At present this implements variable-elimination inference (for conditional 
probability queries P( Y | E = e)) for a bayesian network with discrete 
variables.

I'd like to be able to do:

* some form of learning
* learning with hidden variables
* BN's with all gaussian random variables and related inference
* representation hybrid networks, or at least conditional linear gaussian networks
* inference and learning for the above (ha!)



