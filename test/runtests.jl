#ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

#using SurfaceGeometry
#using Storage

using GeometryTypes
using SurfaceTopology
import LaplaceBIE
using LinearAlgebra

@info "Testing interface"

t = ( 1 + sqrt( 5 ) ) / 2;

vertices = Point{3,Float64}[
    [ -1,  t,  0 ], [  1, t, 0 ], [ -1, -t,  0 ], [  1, -t,  0 ],
    [  0, -1,  t ], [  0, 1, t ], [  0, -1, -t ], [  0,  1, -t ],
    [  t,  0, -1 ], [  t, 0, 1 ], [ -t,  0, -1 ], [ -t,  0,  1 ]
]

faces = Face{3,Int64}[
    [1, 12, 6], [1, 6, 2], [1, 2, 8], [1, 8, 11], [1, 11, 12], [2, 6, 10], [6, 12, 5], 
    [12, 11, 3], [11, 8, 7], [8, 2, 9], [4, 10, 5], [4, 5, 3], [4, 3, 7], [4, 7, 9],  
    [4, 9, 10], [5, 10, 6], [3, 5, 12], [7, 3, 11], [9, 7, 8], [10, 9, 2] 
]

### Why GeometryTypes normal calculation did not work?
n = Point{3,Float64}[]
for v in 1:length(vertices)
    s = Point(0,0,0)
    for (v1,v2) in EdgeRing(v,faces)
        s += cross(vertices[v2],vertices[v1])
    end
    normal = s ./ norm(s)
    push!(n,normal)
end

# msh = HomogenousMesh(vertices,faces)
# n = normals(msh)

hmag = 10.
H0 = [1.,0,0]

vareas = LaplaceBIE.vertexareas(vertices,faces)
psi = LaplaceBIE.surfacepotential(vertices,n,faces,hmag,H0)
# ∇ψ = LaplaceBIE.HField(vertices,n,faces,psi) ### Takes a tangential derivative of the field on the surface
P∇ψ = LaplaceBIE.tangentderivatives(vertices,n,faces,psi)
n∇ψ = LaplaceBIE.normalderivatives(vertices,n,faces,P∇ψ,hmag,H0)

∇ψ = LaplaceBIE.surfacefield(vertices,n,faces,psi,hmag,H0)

energy = LaplaceBIE.fieldenergy(vertices,n,faces,psi,hmag,H0)

@info "Testing precission"

# using JLD
# 

# a,b,c = 2,1,1
# @load "meshes/211-0.2matlab.jld"

# # (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

# mup = 10

# # a,b,c = 2,1/4,1/4
# # (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

# # a,b,c = 2,1/6,1/6
# # (points,faces)=EllipsoidMeshLoad(a,b,c,0.05)

# # a,b,c = 1/4,1,1
# # points,faces = EllipsoidMeshLoad(a,b,c,0.15)

# function ellipsoid_demagnetization_coefficients(a,b,c)

#     UP_BOUND = 1000

#     Ru2(u) = (u+a^2)*(u+b^2)*(u+c^2)

#     nx = 1/2 * a*b*c * quadgk(s -> 1/(s+a^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]
#     ny = 1/2 * a*b*c * quadgk(s -> 1/(s+b^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]
#     nz = 1/2 * a*b*c * quadgk(s -> 1/(s+c^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]

#     return [nx, ny, nz]
# end

# function EllipsoidField(a,b,c,mu,H0)

#     H0x, H0y, H0z = H0
#     nx, ny, nz = ellipsoid_demagnetization_coefficients(a,b,c)

#     Hix = H0x/(1 + (mu-1)*nx)
#     Hiy = H0y/(1 + (mu-1)*ny)
#     Hiz = H0z/(1 + (mu-1)*nz)

#     return [Hix,Hiy,Hiz]
# end

# normals = Array(Float64,size(points)...)
# NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

# for xkey in 1:size(points,2)
#     x,y,z = points[:,xkey]
#     gradf = [x/a^2,y/b^2,z/c^2]
#     nx = gradf/norm(gradf)
#     normals[:,xkey] = nx
# end


# psit = Array(Float64,size(points,2))
# Ht = Array(Float64,3,size(points,2))
# Htheor = EllipsoidField(a,b,c,10,[1,0,0])
# for xkey in 1:size(points,2)
#     nx = normals[:,xkey]
#     psit[xkey] = dot(Htheor,points[:,xkey])
#     Ht[:,xkey] = Htheor
# end

# @time Hn = NormalField(points,faces,mup,[1,0,0],regularize=true,normals=normals)

# normals = Array(Float64,size(points)...)
# NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

# rarr = Float64[]

# Htheor = EllipsoidField(a,b,c,mup,[1,0,0])
# for xkey in 1:maximum(faces)
#     Hnt = dot(Htheor,normals[:,xkey])
#     r = abs((Hn[xkey] - Hnt)/Hnt)
#     term = dot([1,0,0],normals[:,xkey])/mup
#     #importance = abs((Hn[xkey] - term)/term)
#     importance = abs((Hnt - term)/term)
#     #integral error
#     # integrcal = Hn[xkey] - term
#     # integralther = Ht - term
#     bigerr = abs((Hn[xkey] - Hnt)/(Hnt - term))

#     push!(rarr,r)
#     #println("Hn is $(round(Hnt,3)) and computed $(round(Hn[xkey],3)) and r = $(round(r,3)) importance = $(round(importance,3)) ri = $(round(bigerr,3))")
# end

# println("Average relative error is $(mean(rarr)*100) %")

# ########## Some testing

### All field test 

# @time psi = PotentialSimple(points,faces,10,[1,0,0])
# @time rpoints, rfaces = subdivision(points,faces; method=:paraboloid)
# @time Hn = NormalFieldTrapezodial(rpoints,faces,rfaces,10,[1,0,0];NP=3)
# @time H = HField(points,faces,psi,Hn)

# Htheor = EllipsoidField(a,b,c,10,[1,0,0])
# for xkey in 1:25:maximum(faces)
#     r = norm(H[:,xkey] - Htheor)/norm(Htheor)
#     println("r = $r")
# end
