# mfreducible

mfreducible is a gp code that compute reducible residual modular representations.

## Installation

The main function is mfreducible and is situated in the file ``mfred.gp``. Make sure that the value of ``path`` is set to directory where the files are.

````
default(path,"my_path_to_the_files/");
read("mfred.gp");
````

## `mfreducible(f,{flag = 0})`

Compute the list of the prime ideals of the ring of intergers of the coefficient field of the modular newform `f` for which the residual representation of `f` is reducible, and for each prime ideal, the decomposition of the represention in two characters.

````
gp > Red = mfreducible(mfDelta());
gp > #Red
%4 = 5
gp > Red[1]
%5 = [[2, y, 0], [[0, 1], [0, 1], 0, 0]]
gp > Red[2]
%6 = [[3, y, 0], [[0, 1], [0, 1], 0, 1]]
gp > Red[3]
%7 = [[5, y, 0], [[0, 1], [0, 1], 1, 2]]
gp > Red[4]
%8 = [[7, y, 0], [[0, 1], [0, 1], 1, 4]]
gp > Red[5]
%9 = [[691, y, 0], [[0, 1], [0, 1], 0, 11]]
````

If `flag = 0`, the computations are only done for the prime ideals which residue characteristic **does not** divide the index of the irreducible polynomial defining the coefficient field of `f`. When the degree of this polynomial is huge, it speeds up the computations but output only a partial information. However, the output includes also the prime numbers for which the computations could'n have been done, so that the exact list of reducible primes is a subset of the one given by the algorithm.

## Format of the output

### Prime ideals

Let `lambda` be a prime ideal the coefficient field of `f` represented by an irreducible polynomial `Pf`. When the residue characteristic of `lambda` does not divide the index of `Pf`, lambda is given as a  3-component vector `[l,P,gen]` where

* `l` is the residue characteristic of the prime ideal;
* `P` is a polynomial which reduction modulo generates the residue field of `lambda`;
* `gen` is a generator of the residue field.

When the residue characteristic of `lambda` divides the index of `Pf`, `lambda` is a `modpr` PARI format.

### Dirichlet characters

The **primitive** characters given in the output of `mfreducible` are given by a 2-component vector `[a,c]` where

* `c` is the conductor of the character;
* `a` is the Conrey label modulo `c` of the character.

### Output of `mfreducible`

The elements of `mfreducible(f)` are 2-component vectors `[lambda,decomp]` where

* `lambda` is prime ideal as describe above;
* `decomp = [eps1,eps2,m1,m2]` where `eps1` and `eps2` are two primitive Dirichlet characters as above and `m1` and `m2` are two integers.

The residual representation of `f` modulo `lambda` is then isomorphic to `bar(eps1)chi_l^m1 + bar(eps2)chi_l^m2`, where `chi_l` is the cyclotomic character modulo `l`, and `bar(eps1)` and `bar(eps2)` denotes the reductions of `eps1` and `eps2` modulo `lambda` respectively.
