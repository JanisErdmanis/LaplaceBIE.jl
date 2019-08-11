# # Sphere in a homogenous field

# ## Analytics

# To test numerics we can use analytical solution for ellipsoid in a homogenous field.
# In such case the internal field is given through demagnetization coefficients:
# ```math
# H_i = \frac{H_{0i}}{(1 + (\mu-1)n_i)}
# ```

using QuadGK

function ellipsoid_demagnetization_coefficients(a,b,c)

    UP_BOUND = 1000

    Ru2(u) = (u+a^2)*(u+b^2)*(u+c^2)

    nx = 1/2 * a*b*c * quadgk(s -> 1/(s+a^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]
    ny = 1/2 * a*b*c * quadgk(s -> 1/(s+b^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]
    nz = 1/2 * a*b*c * quadgk(s -> 1/(s+c^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]

    return [nx, ny, nz]
end

function EllipsoidField(a,b,c,mu,H0)

    H0x, H0y, H0z = H0
    nx, ny, nz = ellipsoid_demagnetization_coefficients(a,b,c)

    Hix = H0x/(1 + (mu-1)*nx)
    Hiy = H0y/(1 + (mu-1)*ny)
    Hiz = H0z/(1 + (mu-1)*nz)

    return [Hix,Hiy,Hiz]
end

# ## Numerics

# First a spherical mesh is loaded from `SoftSurfaceDynamics.jl` package. We also need normals for the mesh for boundary integral equations which crudelly can also be calculated with `SoftSurfaceDynamics.normals` method. However here more sensible approach is to take normal vectors the vertices themselves.

using SoftSurfaceDynamics 

msh = unitsphere(3)
vertices, faces = msh.vertices, msh.faces
n = vertices 
## The question is whether that is because of normal calculation method fault or due to sphere being incorrectly orriented.


# Now we can define our problem with a field $H_0 = e_z$ and $\epsilon=15$. A corresponding potential for that is $\psi = z$ which we need to pass to the solver.

using LaplaceBIE
using LinearAlgebra
using GeometryTypes
using SurfaceTopology

ϵ = 15
H0 = [0,0,1]

freepotential(x) = dot(x,H0)
freefield(x) = H0

# For finding solution of the problem we execute a line:

ψ = surfacepotential(vertices,n,faces,ϵ,freepotential) 

# which solves a regularized boundary integral equation. And gives the potential on the surface. That then allows to obtain solution everywhere on the space very easally except near the boundary.

# The library offers also to recalculate the pontential gradient on the surface which is needed for example to calculate the surface force due to $(M \cdot \nabla)H$. 

# The tnagential derivatives can be easally claculated from the potential by a finite differences which can be executed.
P∇ψ = tangentderivatives(vertices,n,faces,ψ)

# The normal derivatives are calculated from Biot-Savarat law treating tangential field as raising from a surface current.

n∇ψ = normalderivatives(vertices,n,faces,P∇ψ,ϵ,freefield) 


# ## Comparission

# Since normals are the same as vertices for the nit spehre we expect normal field and potential to be the same. As can be seen from the results they indeed are:

Hin = EllipsoidField(1,1,1,ϵ,H0)

x = [v[3] for v in vertices]
sp = sortperm(x)

for xkey in sp[1:20:end] 
    psit = dot(Hin,vertices[xkey])
    println("$(vertices[xkey][3]) \t $psit \t $(ψ[xkey]) \t $(n∇ψ[xkey])")
end
