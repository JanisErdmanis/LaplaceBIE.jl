# Theory and introduction

Since ancient times people had wondered what makes amber rubbed with fur to attract small light objects and what makes lodestones, naturally magnetized pieces of the minereal magnetite, to atract iron. Knowdays we do have an answer that the force comes from magnetization or polarization gradient. However even for a linear materials the computation of the force is a challenge due to secundary field effects. In this article I propose a numerical method for calculating the field at the objects surface and thus also force for a general shapes as long as thoose are smooth.

To start with we introduce a potential $\psi$ whose gradient is either electric or magnetic field. In absence of singularitities (charges, dipoles and etc.) the potential satisfies Laplace equation $\Delta \psi = 0$ which can also be written in boundary integral form:
```math
\int_{\partial \Omega} \frac{\partial \psi}{\partial n_{\boldsymbol x}} \frac{1}{|\boldsymbol x - \boldsymbol y|} dS_{\vec{x}}
-\int_{\partial \Omega} \psi(\boldsymbol x) \frac{\partial}{\partial n_{\boldsymbol x}} \frac{1}{|\boldsymbol x - \boldsymbol y|} dS_{\vec{x}}
=\left\{
    \begin{array}{ll}
      4\pi \psi(\boldsymbol y) ~&if~\vec{y} \in \Omega \\
      2 \pi \psi(\boldsymbol y)~&if~\vec{y} \in \partial \Omega \\
      0 ~&if~\vec{y} \ni \Omega 
    \end{array}
  \right.
```
This integral is flexible in the sense that we can wrap the surface around singularities, infinity and objects and we can apply it to the interior of objects with linear properties.

Let's consider the system shown in the figuere:

[figure]()

We have singularites denoted by $S_n$, object $\Omega$ and infinity $\infty$. Since at infinity the field perturbation of the object is vanishingly small we can set it equal to free field. Similarly the field near the singularity is equal to the field of the singularity and we can neglect effects of the body. Thus we have identity:
```math
  \int_{\partial \infty} + \int_{\partial S} = 4 \pi \psi_{0} (\vec y)
```
where left side is integral over surface of singularities and infinity and the right side is the field in absence of bodies.

To model the bodies inside the field we only need to set a appropriate boundary conditions which are:

This particular package is built with a following boundary conditions
```math
\psi^{in} = \psi^{out} = \psi; \mu \frac{\partial \psi^{in}}{\partial n} = \frac{\partial \psi^{out}}{\partial n} = \frac{\partial \psi}{\partial n}
```
which implies a linear material. Applying the boundary integral to the interior of the object and matching boundary conditions with exterior surface we arrive at useful formula:
```math
\int_{\partial \Omega} \frac{\partial \psi^{out}}{\partial n_x} \frac{1}{|\vec x - \vec y|} dS_{\vec x} = \frac{2 \mu}{\mu + 1} \int_{\partial \Omega} \psi(\vec x) \frac{\partial}{\partial n_{x}} \frac{1}{|\vec x - \vec y|} dS_{\vec x}
```
Which we can apply to obrtain boundary integral equation for the potential
```math
  \psi(\vec y) = 2 \psi_{0}(\vec y)  + \frac{\mu - 1}{\mu + 1} \int_{\partial \Omega} (\psi(\vec x)-\psi(\vec y)) \frac{\partial}{\partial n_x} \frac{1}{|\vec x - \vec y|} dS_{\vec x}
```
where we used identity $\int_{\partial \Omega} \frac{\partial}{\partial n_x} \frac{1}{|\vec x - \vec y|} dS_{\vec x}=0$ for regularization.

Another way of writting Laplace equation in BIE form is as follows:
```math
\int_{\partial \Omega} (\nabla_{\vec x} \psi \times \vec n_{\vec x}) \nabla_x \frac{1}{|\boldsymbol x - \boldsymbol y|} dS_{\vec{x}}
-\int_{\partial \Omega} \frac{\partial \psi}{\partial n_{\vec x}} \nabla_{\vec x}  \frac{1}{|\boldsymbol x - \boldsymbol y|} dS_{\vec{x}}
=\left\{
    \begin{array}{ll}
      4\pi \nabla \psi(\boldsymbol y) ~&if~\vec{y} \in \Omega \\
      2 \pi \nabla \psi(\boldsymbol y)~&if~\vec{y} \in \partial \Omega \\
      0 ~&if~\vec{y} \ni \Omega 
    \end{array}
  \right.
```math

From a similar arguments we obtain
```math
  \int_{\partial \infty} + \int_{\partial S} = 4 \pi \nabla \psi_{0} (\vec y)
```
Applying boundary conditions give us a formula
```math
\vec n_{\vec y} \cdot \int \frac{\partial \psi}{\partial n} \nabla_{\vec x} \frac{1}{|\vec{x} - \vec{y}|} dS_{\vec x} = \frac{\mu + 1}{2} \vec n_{\vec y} \int (\nabla \psi \times \vec n_{\vec x}) \nabla_{\vec x} \frac{1}{|\vec{x} - \vec{y}|} dS_{\vec x}
```
which we substitute back to obtain:
```math
  \frac{\partial \psi^{in}}{\partial n_{\vec y}} = 2 \nabla \psi_0 \cdot \vec n_{\vec y} - \frac{\mu -1}{4 \pi} \vec n_{\vec y} \cdot \int (\nabla \psi \times \vec n_{\vec x}) \times \nabla_{\vec x} \frac{1}{|\vec{x} - \vec{y}|} dS_{\vec x}
```
It is interesting to note that it is a Biot-Svarat integral over the surface current which models the field configuration in presence of object. A regularized version can be obtained (look into article on magnetic liquid droplet)
```math
  \frac{\partial \psi^{in}}{\partial n_{\vec y}} = 2 \nabla \psi_0 \cdot \vec n_{\vec y} - \frac{\mu -1}{4 \pi} \vec n_{\vec y} \cdot \int (P\nabla \psi(\vec x) - P \nabla \psi (\vec y) ) \times \left( \vec n_{\vec x} \times \nabla_{\vec x} \frac{1}{|\vec{x} - \vec{y}|} \right) dS_{\vec x}
```

# Implementation and API

The boundary integral equations are solved with collocation metheod on a triangular mesh. The simplest quadrature $\int_\Delta f dS = (f_1 + f_2 + f_3)/3$ is used for all nonsingular elements. In order to obtain tangential derivatives of the potential we use numerical differentiation of previosly calculated potential and use least squares to solve the overdetermined system. And lastly for Biot-Savarat integral weakly singular quadrature is implemneted to deal with elements near the singularity where for other elements $\int_\Delta f dS = (f_1 + f_2 + f_3)/3$ is used.

```@autodocs
 Modules = [LaplaceBIE]
```
