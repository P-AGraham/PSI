## Relevant forum posts

# $\mathbb{Z}_2$ symmetry
Saying a quantum system has a *discrete* $\mathbb{Z}_2$ symmetry is the same as saying that there exists a symmetry operator $S$ which squares to identity and commutes with the hamiltonian $H$ while not being a generator of a continuous symmetry (which could include $\mathbb{Z}_2$ as a subgroup)
https://physics.stackexchange.com/questions/497695/naming-symmetries-in-quantum-systems-e-g-mathbbz-2-or-u1


# Mermin-Wagner theorem
https://physics.stackexchange.com/questions/478518/mermin-wagner-and-graphene


# Spontaneous symmetry breaking
SSB is usually associated to a degeneracy in the ground state of a system having a continuous symmetry. At low temperature, the system can localize the expectation value for its dynamical symmetry breaking field in the field value sub-manifold minimizing the potential. The Mermin-Wagner theorem theorem describes if this state can stay localised in this sub-manifold or if the goldstone bosons (excitations of the field specifically in the sub-manifold where the potential is flat leaving only the massless kinetic term in the lagrangian) fluctuate too much leading to an impossibility of long range ordering (order means symmetry breaking) at a specific point in the sub-manifold. Quantum tunneling trough the potential barrier defining the ground state sub-manifold is possible, but it is a supressed effect with growing system size: the average value for the field at a given point $x$ is the same everywhere for a transationally symmetric theory. The bigger the volume is the longer it will take for a field prepared in a a localized region of the sub-manifold to tunnel to the other side in a substancial amount of points to make a difference in the expectation value.  
https://philsci-archive.pitt.edu/14983/1/ssbfinite.pdf

## Papers


## Ideas

1. Instead of icosaedron use stars (multiple covers of the 2-sphere) https://en.wikipedia.org/wiki/Star_polyhedron
2. Write the set of all finite dimensionnal systems on which SO(3) is represented. Use these states to form models with Z2 symmetry breaking (ordering field). 
3. Fuzzyfy the fuzzy sphere: put a magnetic monopole at the center of the fuzzy sphere and implement it in the hamiltonian with the same potnetial vector but where the position operator has a truncated representation (containing only monopole harmonics). This amounts to projecting the Hamiltonian for the original sphere to the monopole harmonic basis. Since it is an eigenbasis of the hamiltonian, this projection does not change the eigen states and the fuzification is an idempotent map in this context.


# Code 

1. The $\mathbb{Z}_2$ symmetry operator $S$ has to be diagonalized simultaneously to $H$ and $L^2$ which can be acheived by instead diagonalizing $H + s S + t L^2$ (see https://math.stackexchange.com/questions/4387863/simultaneous-diagonalization-and-repeated-eigenvalues) with choices for $s$ and $t$ that lift all degeneracies in the spectrum. We know the spectrum of $L^2$ ($\ge 1$), $S$ ($\in \{-1, 1\}$) and $H$. We can pick $s = 1/2$ and $t = 1$ if we are certain that the hamilonian has no half integer gap between its eigenvalues. We also want to preserve the order of the eigenvalues which is not necessarly done if $S$ shifts by $1/2$. The logic goes as follows: identify smalest non-zero gap in $H$ (noted $a$), $L^2$ ($j(j+1) - (j-1)j = 2j \to 2$) and $Z$ (2). Then the smalest non-zero gap in $H + tL^2$ is $a + 2t$ where we pick $2t < a$. Then the smallest non-zero gap in $H + t L^2 + s S$ is $a + 2t + 2s$ and we want to have $2s < a + 2t$. A consistant choice is $t = a/4$ and $s = (a/2 + t)/2 = (a/2 + t)/2  = 2a/8 + a/8 = 3a/8$. So we diagonalize $H$, find $a$, calculate $H + s S + t L^2$, diagonalize it than remove the eigenvalues of $H$ to find eigenvalues of $s S + t L^2$ which are $V = \pm s + t j(j + 1)$. To extract $j$ try $j^2 + j - (V \mp s)/t = 0$: $j = -1/2 \pm \sqrt{1 + 4 (V \mp s)/t}/2 \to -1/2 + \sqrt{1 + 4 (V \mp s)/t}/2$ (positive root because $j > 0$). Now we need one of the solutions to be an integer. To simplify, we use $V/a$ instead $t \to 1/4$. We then have $1 + 16 (V \mp 3/8)  = 1 \mp 6 + 16 V = (\text{odd})^2$. Supposing the two solutions are valid, we get $(\text{odd})^2 - (\text{odd})^2 = 12 = (2 k + 1)^2 - (2 m + 1)^2 = 4 (k^2 + k - m^2 - m)\implies 3 = (k^2 + k - m^2 - m)$. 