# Ground state of the anisotropic Shastry-Sutherland model


Computes the ground state wavefunction for the the anisotropic Shastry-Sutherland interaction $J_1, J2, ..$. We assume the ground state to have a $N\times N$ unit cell in the dimer basis. 

The optimisation of the local tensors is done using the simple update. 

By mapping the dimer to local tensors, the wavefunction is written as an infinite tensor network on the square lattice. The overlap and observables are then computed by using the CTMRG algorithm.

To start the environement 
```
julia --project=.
using Pkg
Pkg.instantiate("ShastryGS")
```


The tensors are stored in the /Results folder in a jld2 format, while the observables are stored in the /Results folder in .mat format.

Two different models are possible, the $XY$ model

$$
H = \sum J_1 (S_x S_x + \Delta (S_y S_y + S_z S_z)) + h_x S_x + h_z S_z
$$

or $XYZ$ model

$$
H = \sum J_1 (Jx S_x S_x + Jy S_y S_y + Jz S_z S_z) + h_x S_x + h_z S_z
$$

To run the simulations, 

```
include("main.jl")
```

or 
```
sbatch SSM.run
```

for HPC environments
