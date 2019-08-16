# # Sphere in a point charge field.

# ## Analytics and the problem

# In general the LaplaceBIE allows to find field on the object if the potential in absence of objects is known. To demonstrate that let's consider a field on the dielectric sphere (radius $r_1$) due to a point charge at distance $\zeta$ away. The solution inside the sphere is given as series [Straton page 221]:
# ```math
# \psi = \sum_{n=0}^{\infty} \frac{2n + 1}{\epsilon n + n + 1} \frac{r_1^n}{\zeta^{n+1}} L_n(\cos \theta)
# ```
# where $L_n$ is Legendre polynomial and $\theta$ angle from line connecting point charge and spehre.

# The normal derivative can also be given for the sphere which coincides with derivative in the radial direction
# ```math
# \frac{\partial \psi}{\partial n} = \sum_{n=0}^{\infty} \frac{n(2n + 1)}{\epsilon n + n + 1} \frac{r_1^{n-1}}{\zeta^{n+1}} L_n(\cos \theta)
# ```

using Jacobi

function ψt(cosθ,ϵ,ζ,r1)
    s = 0
    for n in 0:25
        s += (2*n + 1)/(n*ϵ + n + 1)*r1^n/ζ^(n+1)*legendre(cosθ,n)
    end
    return s
end

function ∇ψn(cosθ,ϵ,ζ,r1)
    s = 0
    for n in 1:25
        s += n*(2*n + 1)/(n*ϵ + n + 1)*r1^(n-1)/ζ^(n+1) * legendre(cosθ,n)
    end
    return s
end

# ## Numerics

# First we need to load a sphere mesh. 

using LinearAlgebra
using GeometryTypes
using SurfaceTopology
using LaplaceBIE

include("sphere.jl")
msh = unitsphere(2)
vertices, faces = msh.vertices, msh.faces
n = vertices

# The next step is to set up the free field solution which we do as follows

ζ = 1.2
r1 = 1
y = [0,0,ζ]

freepotential(x) = 1/norm(x.-y)
freefield(x) = -(x.-y)/norm(x.-y)^3

# Then in a presence of a sphere with $\epsilon=10$ we have a solution:

ϵ = 10
ψ = surfacepotential(vertices,n,faces,ϵ,freepotential) 
P∇ψ = tangentderivatives(vertices,n,faces,ψ)
n∇ψ = normalderivatives(vertices,n,faces,P∇ψ,ϵ,freefield) 

# ## Comparision

using Winston

cosθ = range(-1,1,length=100)
cosθs = [x[3]/1 for x in vertices]

# Potential test

scatter(cosθs,ψ)
oplot(cosθ,ψt.(cosθ,ϵ,ζ,r1) )

# Test of normal derivatives derivatives

scatter(cosθs,n∇ψ)
oplot(cosθ,∇ψn.(cosθ,ϵ,ζ,r1) )
