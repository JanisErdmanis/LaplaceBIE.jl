# # A sphere in a homogenous field

# For a three axial ellipsoid in a constant field, we have an analytic solution expressed in terms of demagnetization coefficients $n_i$:
# ```math
# H_i = \frac{H_{0i}}{(1 + (\mu-1)n_i)}
# ```
# whre $\vec H$ is the interior field, and $\vec H_{0}$ is the constant field. That makes ellipsoid a perfect subject to test the numerics.

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
end;

# We define that we have a field of strength $H_0 = e_z$. A corresponding potential for that is $\psi_0 = z$, which we need to pass to the solver.

using LaplaceBIE
using LinearAlgebra

H0 = [0,0,1]

ψ0(x) = dot(x,H0)
∇ψ0(x) = H0;

# For simplicity, we consider a sphere with relative permittivity $\epsilon=15$. The mesh for the sphere we generate from subdivisions of icosahedron given in sphere.jl. We also need normal vectors which we get from the vertex positions with high accuracy

ϵ = 15

include("sphere.jl")
msh = unitsphere(3)
vertices, faces = msh.vertices, msh.faces
n = vertices;

# At the moment users need to take the derivative of the potential themselves which can be done in multiple ways. It is also possible that a user might wish to calculate the field due to the surface current for which `normalderivatives` method could be useful.

# Now we can proceed with calculation. To calculate surface potential everywhere on the surface one executes a method
ψ = surfacepotential(vertices,n,faces,ϵ,ψ0);
# which solves a regularized boundary integral equation with BLAS. That then, for example, can be used to calculate the energy of the field.

# Usually, however, one wants to know the field on the surface for force calculations due to $(M \cdot \nabla)H$. The library offers to do a finite differentiation on the calculated potential for calculating tangential field components. That can be easily achieved with method `tangentderivatives`:

P∇ψ = tangentderivatives(vertices,n,faces,ψ);

# And lastly to calculate the normal derivatives, we can use Biot-Savarat law treating the tangential field as raising from a surface current:

n∇ψ = normalderivatives(vertices,n,faces,P∇ψ,ϵ,∇ψ0);

# We can compare numerics with analytics easily by the use of azimuthal symmetry in the z-direction. The normal field is equal to potential since vertex positions coincide with normals which is thus not shown twice. And so the comparison:

Hin = EllipsoidField(1,1,1,ϵ,H0)

x = [v[3] for v in vertices]
sp = sortperm(x)

for xkey in sp[2:43:end] 
    psit = dot(Hin,vertices[xkey])
    println("$(vertices[xkey][3]) \t $psit \t $(ψ[xkey]) \t $(n∇ψ[xkey])")
end

