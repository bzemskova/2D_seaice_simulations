# module cuda/12.6
# module julia

#using Pkg; Pkg.add("CUDA")
using Oceananigans
using Oceananigans.Units
using SpecialFunctions
using CUDA

# Define parameters

# domain
Nx, Nz = 4096, 4096
const H = 350meters
const L = 1200meters

# non-dimensional parameters
const zeta = 0.033
const eta = 55.
const Ri = 3.6
const chi = 0.52
const J = 0.21

# fixed dimensional values
const width = 40meters      # \sigma parameter of the versoria function to define ice keel [m]
const U0 = 0.1              # The steady component of the flow in m/s
const k0 = 0.5*pi/width     # wavenumber of the ice keel (relating sinusoidal and versoria functions) [1/m]

# calculate other parameters based on fixed values and non-dimensional parameters
const N = U0*k0/chi         # buoyancy frequency [1/s]
const deltaB = Ri*k0*(U0^2) # buoyancy jump across the pycnocline [m/s^2]
const h₀ = pi*zeta/k0       # ice keel height [m]
const z0 = -eta*h₀          # mixed layer depth [m]

# other stratification parameters
const mu = 0.5  
const N² = N^2  # [s⁻²] buoyancy frequency / stratification

# simulation time parameters
Δt = 0.006seconds
damping_rate = 1/(20*Δt)
stop_time = 6hours
save_fields_interval = 5minutes
chkpt_iterations = 100000
filename = "cluster0_4096_4096"
chkpt_prefix = "chkpt_cluster0"

###################################################################
# Setting up the grid
###################################################################

#The following code is used to create a non-linear grid spacing with
#faces znew_faces. This grid spacing creates a higher
#resolution near the hill.

#zscaler(x)=erf((x-z0+(H/2))/(H/8))+2

#zscale_sum=0
#for k in 1:Nz
 #   global zscale_sum=zscale_sum+zscaler(-H+(k-1)*(H/Nz))
#end

#zscaling_factor=H/zscale_sum

#zface_array=Array{Float64}(undef,Nz+1)
#zface_array[1]=-H
#for i in 2:(Nz+1)
 #   zface_array[i]=zface_array[i-1]+zscaling_factor*zscaler(-H+(i-2)*(H/Nz))
#end

#znew_faces(k)=zface_array[trunc(Int,k)]


underlying_grid = RectilinearGrid(GPU(),size = (Nx, Nz),
                                  x = (-L, L),
                                  z = (-H, 0),
                                  halo = (4, 4),
                                  topology = (Periodic, Flat, Bounded))


ice(x) = h₀ * width^2/(width^2 + 4*x^2)
bottom(x) = - H + ice(x)

#grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom))
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom))

##########################################################################
# Setting up forcing, ICs, sponge layer
##########################################################################

coriolis = FPlane(latitude = 70)

# Forcing require to maintain steady flow
const mean_forcing_amplitude = (coriolis.f)*U0         
@inline mean_forcing(x, z, t) = mean_forcing_amplitude

# Initial conditions
uᵢ(x, z) = U0
bᵢ(x, z) = 0.5*(deltaB - N² * ((-H-z0)-z)) *(1-tanh(((-H-z0)-z)/mu))

function bottom_mask_func(x, z)
    sponge_width = -H/20
    return exp(-z^2 / (2 * sponge_width^2)) #+ (erf((x-L/1.3)/(L/16)) + erf((-x-L/1.3)/(L/16)) +2)/2
end

function target_b_func(x, z, t)
    return 0.5*(deltaB - N² * ((-H-z0)-z)) *(1-tanh(((-H-z0)-z)/mu))
end

u_sponge = Relaxation(rate=damping_rate, mask=bottom_mask_func, target = U0)
w_sponge = Relaxation(rate=damping_rate, mask=bottom_mask_func)
b_sponge = Relaxation(rate=damping_rate, mask=bottom_mask_func, target=target_b_func)

model = NonhydrostaticModel(; grid, coriolis,
                            advection = Centered(order=4),
                            closure = ScalarDiffusivity(ν=1e-3, κ=1e-3),
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            forcing = (; u = u_sponge, v=mean_forcing, w=w_sponge, b=b_sponge))


###############################################################################################
# Run the simulation
##############################################################################################


set!(model, u=uᵢ, b=bᵢ)
simulation = Simulation(model; Δt, stop_time)

b = model.tracers.b
u, v, w = model.velocities
N2 = ∂z(b)
#U = Field(Average(u))
#u′ = u - U
simulation.output_writers[:fields] = JLD2Writer(model, (; u, w, b, N2);
                                                filename,
                                                schedule = TimeInterval(save_fields_interval),
                                                overwrite_existing = true)

simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=IterationInterval(chkpt_iterations), prefix=chkpt_prefix)
run!(simulation)
