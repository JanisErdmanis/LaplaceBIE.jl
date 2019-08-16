using GeometryTypes
using SurfaceTopology
using LaplaceBIE
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

hmag = 10.
H0 = [1.,0,0]

vareas = vertexareas(vertices,faces)
psi = surfacepotential(vertices,n,faces,hmag,H0)
P∇ψ = tangentderivatives(vertices,n,faces,psi)
n∇ψ = normalderivatives(vertices,n,faces,P∇ψ,hmag,H0)

∇ψ = surfacefield(vertices,n,faces,psi,hmag,H0)

energy = fieldenergy(vertices,n,faces,psi,hmag,H0)

