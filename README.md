# Ground state of the anisotropic Shastry-Sutherland model

Computes the ground state wavefunction for the anisotropic Shastry-Sutherland model with iPEPS and simple update in the dimer basis. 

$N$ is the length of the unit cell in the dimer basis, it can take N = (2,6,8). $D$ is the bond dimension of the local tensor, $J_2$ is the intra dimer interaction, $J_1$ the interdimer interaction. 

Two different models ar possible, the XY 

$$
H = /sum J_1 (S_x S_x + \Delta (S_y S_y + S_z S_z)) + h_x S_x + h_z S_z
$$
and XYZ model

$$
H = /sum J_1 (Jx S_x S_x + Jy S_y S_y + Jz S_z S_z) + h_x S_x + h_z S_z
$$

