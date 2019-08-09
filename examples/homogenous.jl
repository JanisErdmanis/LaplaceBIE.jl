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


using Jacobi
using Winston

ϵ = 15

### Now I need to do the same thing nummerically

using QuadGK
using LinearAlgebra
using GeometryTypes
using SurfaceTopology
using SoftSurfaceDynamics # used oly to load meshes
using LaplaceBIE

msh = unitsphere(3)
vertices, faces = msh.vertices, msh.faces

n = SoftSurfaceDynamics.normals(vertices,faces) 

H0 = [0,0,1]

freepotential(x) = dot(x,H0)
freefield(x) = H0

ψ = surfacepotential(vertices,n,faces,ϵ,freepotential) # H0 could be a coordinate dpeendant function 
P∇ψ = tangentderivatives(vertices,n,faces,ψ)
n∇ψ = normalderivatives(vertices,n,faces,P∇ψ,ϵ,freefield) # H0 could be a coordinate dpeendant functiOn

# cosθs = [x[3]/1 for x in vertices] ### if one puts point like charge on the z axis

# scatter(cosθs,n∇ψ)

Hin = EllipsoidField(1,1,1,ϵ,H0)

x = [v[3] for v in vertices]
sp = sortperm(x)

for xkey in sp #1:length(vertices)
    Hnt = dot(Hin,vertices[xkey])
    psit = dot(Hin,vertices[xkey])
    Htan = norm(P∇ψ[xkey])
    Htant = sqrt(norm(Hin)^2-Hnt^2)
    #println("$(vertices[xkey][3]) $Hnt")
    #println("$(vertices[xkey][3]) \t $psit \t $(ψ[xkey]) \t $Htant \t $Htant  \t $(n∇ψ[xkey])")
    println("$(vertices[xkey][3]) \t $psit \t $(ψ[xkey]) \t $(n∇ψ[xkey])")
end



#cosθ = range(-1,1,length=100)
#oplot(cosθ, -Hin[3].*cosθ*ϵ)   ### There


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
