# The distributed mpc package

This package implements the distributed MPC logic within the Game-Theoretic framework, as proposed in "Rawlings, James Blake, David Q. Mayne, and Moritz Diehl. Model predictive control: theory, computation, and design. Vol. 2. Madison, WI: Nob Hill Publishing, 2017.".

This repository implements the two player linear game.

In brief, the distributed model predictive control addressed by this repository solves the following minimization problem, in the same way, for the cooperative and the non-cooperative cases.

$$min_{u,i} \, J_i(k) = \sum_{l=1}^{N}  e_i(k+l)^T\,Q_{i,i}\,e_i(k+l) + e_j(k+l)^T\,Q_{i,j}\,e_j(k+l) +  + u_i(k+l)^T\,R_i\,u_i(k+l)$$
subject to
$$z(k+1) = A\,z(k)+B_{1} u_1(k)+B_{2} u_2(k)$$

$$y(k) = C \, z(k)$$


### TODO
1. better documentation
2. purge from ros

## Contacts
if you have any questions about this repository, are interested in applications, want to contribute or anything else,
please contact me at paolo.franceschi@supsi.Cauchy-Schwarz

### Acnowledgements
This repository was developed within the [EU Horizon 2020 DrapeBot project](https://www.drapebot.eu/) context. 

