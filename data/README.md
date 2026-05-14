### last timestep data: clusterN_last_timestep.npz
contains 2D arrays for velocity $u$ (variable 'utot') and $w$ (variable 'w') and buoyancy $b$ (variable 'b')

### grid data: clusterN_grid.npz

x_array: array of x-coordinate values

z_array: array of z-coordinate values

bh_array: location of the immersed boundary (keel/topography)

Z_masked: grid masking due to ice keel/topography

Nx, Nz: number of grid points in horizontal and vertical directions

L, H: horizontal and vertical extents of the domain

J, eta, Fr, J, chi: nondimensional parameters for the simulation set up

U0, k0, h0, z0, N: dimensional parameters for the simulation set up
