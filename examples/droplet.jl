using LinearAlgebra
using GeometryTypes
using SurfaceTopology
using LaplaceBIE
using SoftSurfaceDynamics
using ElTopo

msh = unitsphere(2)
vertices, faces = msh.vertices, msh.faces

########### Original code from MDrop #########

function energy(points,normals,faces,psi,mup,H0)
    vareas = vertexareas(points,faces)
    Area = sum(vareas)

    s = 0
    for xkey in 1:length(points)
        s += psi[xkey]*dot(H0,normals[xkey]) * vareas[xkey]
    end

    Es = gammap * Area
    Em = 1/8/pi * (1 - mup) * s

    return Es+Em
end

H0 = [1.,0.,0.]
etap = 1.
gammap = 1.
μ = 2.

t = 0.
Δt = 0.1
volume0 = surfacevolume(vertices,faces)

for i in 1:10
    n = SoftSurfaceDynamics.normals(vertices,faces) 

    psi = surfacepotential(vertices,n,faces,μ,H0)
    P∇ψ = tangentderivatives(vertices,n,faces,psi)
    Hn = normalderivatives(vertices,n,faces,P∇ψ,μ,H0)

    E = energy(vertices,n,faces,psi,μ,H0)
    rV = surfacevolume(vertices,faces)/volume0
    @show E,rV

    Ht = [norm(j) for j in P∇ψ]
    tensorn = μ*(μ-1)/8/pi * Hn.^2 + (μ-1)/8/pi * Ht.^2

    vn = stokesvelocity(vertices,n,faces,tensorn,etap,gammap)

    vertices .+= n .* vn * Δt

    global t += Δt

    ### ElTopo stabilization
    par = SurfTrack(allow_vertex_movement=true)
    msh = stabilize(HomogenousMesh(vertices,faces),par)
    global vertices, faces = msh.vertices, msh.faces
end

# [norm(i) for i in vertices]
