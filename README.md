# LaplaceBIE.jl

Boundary integral equation methods for solving Laplace's equation on the interface.

## Applications

   + Calculate the surface deforming force and find equilibrium shape either in electric or magnetic field. See `examples/droplet.jl`.
   + Calculate how homogenous object perturbs the field and energy for it to be put into homogenous field.
   + Calculate the force exerted on the body due to a field gradient. At the moment requires modifications in the code (making H0 vertex dependant) and some test cases if that works.

