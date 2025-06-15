# Ground state of the anisotropic Shastry-Sutherland model


Computes the ground state wavefunction for the anisotropic Shastry-Sutherland model with interaction $J_1, J_2, \ldots $. We assume the ground state to have a stripped $N\times N$ unit cell in the dimer basis. We illustrate the unit cell for the case $N = 6$ in the /Figures folder.

The optimisation of the local tensors is done using the simple update. 

By mapping the dimer to local tensors, the wavefunction is written as an infinite tensor network on the square lattice. The overlap and observables are then computed by using the CTMRG algorithm.

To start the environement with julia 1.8
```
julia --project=.
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```


The tensors are stored in the /Results folder in a jld2 format, while the observables are stored in the /Results folder in .mat format.

Two different models are possible, the $XY$ model

$$
H = \sum J_1 (S_x S_x + \Delta (S_y S_y + S_z S_z)) + h_x S_x + h_z S_z
$$

or $XYZ$ model

$$
H = \sum J_1 (Jx S_x S_x + Jy S_y S_y + Jz S_z S_z) + h_x S_x + h_z S_z.
$$

To run the simulations, type 

```
include("main.jl")
```
this perform the imaginary time evolution with time-step $d\beta$ and parameters defined in the same file. If the simulation with those parameters has already been performed for bond dimension $D$ or $D-1$, we will use the results of the previous simulation for the initial state. Otherwise, the simulations will start with a random product state. 

To run the code on HPC environments, run

```
sbatch SSM.run
```

or

```
sbatch submit_bash.sh
```
