# racket-pgm

Work in Progress.  This is not a license.

An ad-hoc, informally specified, bug-ridden, slow implementation of half of 
python's libpgm.

The author felt that it would be instructive (to him) to implement 
representation, inference and learning for Bayesian networks, and that it would 
be fun to do so in typed racket.

At present this implements variable-elimination inference (for conditional 
probability queries P( Y | E = e)) for a conditional gaussian Bayesian network.
That means no discrete variable can have continuous parents.

The internal representation for factors in variable elimination is tables of
mixtures of canonical representations of normal distributions.  The author is 
aware that these have an exponential space requirement and the library would
be unusable for many models.  This project is thus a "toy library" and only 
useful for fun-sized networks.

This project also implements Maximum Liklihood parameter learning from data.  

The author would this project/package to be able to do:

* EM-MLE parameter learning over incomplete data.
* Forward sampling of networks (and related approx inference).

The author would like to implement other inference algorithms besides VE, as he
feels they would be equally instructive.  Alas, the author does not currently
understand many of them, while others are just frightening.  These including:
* believe propagation
* MCMC
* variational inference

Additionally, the author would like to implement Bayesian learning methods, but
he is afraid of them especially in the case of missing or incomplete data.

If the reader scoffs at the lack of intellectual intrepidness of the author, then the
reader is of course invited to contribute. :)

Example datasets come from: 

    Lichman, M. (2013). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. Irvine, CA: University of California, School of Information and Computer Science.
    
    Myself (2017). Trying to generate a dataset that recreates via MLE parameter learning
    the netowrk given in Koller, Friedman, Probabilistic Graphical Models.
