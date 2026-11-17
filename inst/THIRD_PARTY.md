# Algorithm attribution

BayesPrism 2.2.3, Tinyi Chu and Charles G. Danko, Danko-Lab/BayesPrism,
GPL-2 OR GPL-3, is an actual required dependency. Its package is not copied or
renamed into this package. The C++ initial-allocation execution kernel implements
the same multinomial allocation / symmetric-Dirichlet Gibbs updates and is tested
against `BayesPrism::sample.Z.theta_n`. The optional typed-observation adapter is
the existing audited project implementation, including the exact Metropolis
proposal correction; it is not described as native upstream BayesPrism flow support.

Attribution: Chu et al., Cell type and gene expression deconvolution with
BayesPrism enables Bayesian integrative analysis across bulk and single-cell
RNA sequencing in oncology. Nature Cancer (2022).
https://doi.org/10.1038/s43018-022-00356-3

The wrappers, reference management, observation audit, hierarchy/source exports
and benchmark code are project code. The GPL-3 source distribution preserves
this attribution and does not imply an independent new statistical model.
