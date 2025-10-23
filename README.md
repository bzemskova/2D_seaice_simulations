# 2D_seaice_simulations
Code repository to support manuscript titled "Nondimensional parameter regimes of Arctic ice keel-ocean flow interactions"

### Model set-up

Rotating incompressible Navier-Stokes equations with Boussinesq approximation are solved in a two-dimensional $$(x,z)$$ domain using Oceananigancs. In the current implementation of Oceananigans, an immersed boundary can only be specified along the bottom, and not the top boundary. However, in the Boussinesq approximation, flows are symmetric when flipped upside down if the sign of buoyancy is also reversed. We use this property to implement an "ice keel" along the bottom of the domain and reverse the sign of the initial vertical buoyancy profile.

- The sea ice keel shape is modeled as a Versoria shape:

  $$h(x) = \frac{h\sigma^2}{\sigma^2+ 4x^2}$$

where $$h$$ is the depth/height of the keel and $$\sigma$$ is the width of the keel

- 
  
