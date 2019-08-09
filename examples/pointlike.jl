using Jacobi
using Winston

### Probably the best option is to compare the formula in the limit. 

function ψt(cosθ,ϵ,ζ,r1)
    s = 0
    for n in 0:25
        #s += -(ϵ-1)*n/(n*(ϵ+1) + 1)*ab^n*legendre(cosθ,n)
        s += (2*n + 1)/(n*ϵ + n + 1)*r1^n/ζ^(n+1)*legendre(cosθ,n)
    end
    return s
end

function ∇ψn(cosθ,ϵ,ζ,r1)
    s = 0
    for n in 1:25
        #s += n^2/(n*(ϵ+1) + 1) * (ab)^(n-1) * legendre(cosθ,n)
        s += n*(2*n + 1)/(n*ϵ + n + 1)*r1^(n-1)/ζ^(n+1) * legendre(cosθ,n)
    end
    return s
end

### The problem perhaps is that when one calculates the potential dot(H0[xkey],x) is no longer a potential!!! 

ϵ = 10

#plot(cosθ,∇ψn.(cosθ,ϵ,ab))

### Now I need to do the same thing nummerically

using LinearAlgebra
using GeometryTypes
using SurfaceTopology
using SoftSurfaceDynamics # used oly to load meshes
using LaplaceBIE

msh = unitsphere(2)
vertices, faces = msh.vertices, msh.faces

n = SoftSurfaceDynamics.normals(vertices,faces) 

ζ = 1.2
r1 = 1

y = [0,0,ζ]

freepotential(x) = 1/norm(x.-y)
freefield(x) = -(x.-y)/norm(x.-y)^3

ψ = surfacepotential(vertices,n,faces,ϵ,freepotential) # H0 could be a coordinate dpeendant function 
P∇ψ = tangentderivatives(vertices,n,faces,ψ)
n∇ψ = normalderivatives(vertices,n,faces,P∇ψ,ϵ,freefield) # H0 could be a coordinate dpeendant function

cosθs = [x[3]/1 for x in vertices] ### if one puts point like charge on the z axis

### Let's begin with potential
scatter(cosθs,ψ)

# Ptential test
# cosθ = range(-1,1,length=100)
# oplot(cosθ,ψt.(cosθ,ϵ,ζ,r1) )

### Normal derivatives
scatter(cosθs,n∇ψ)
oplot(cosθ,∇ψn.(cosθ,ϵ,ζ,r1) )
