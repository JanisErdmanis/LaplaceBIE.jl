# LaplaceBIE.jl theory (Magnetic or electric force for a smooth bodies)

Since ancient times people had wondered what makes amber rubbed with fur to attract small light objects and what makes lodestones, naturally magnetized pieces of the minereal magnetite, to atract iron. Knowdays we do have an answer that what makes the force is either magnetization or polarization gradient. However even for a linear materials the computation of the force is a challenge due to secundary field effects. In this article I propose a numerical method for calculating the field and thus also force for a general shapes as long as thoose are smooth.

To sater with we introduce a potential $\psi$ whose gradient is either electric or magnetic field which satisfies Laplace equation $\Delta \psi = 0$. This equation can be written in a boundary integral form:
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

It is particularly usefull to choose a boundary where one of the interfaces is the object and other reaches the infinity plus some interface of some singular regions. Seperatelly for thoose regions we have a following identities which are tro at the ssurface of the object:
```math
  \int_{\partial \infty} = 4 \pi \psi_{\infty} (\vec y)
```

```math
\int_{\partial s} = 4 \pi \psi_s (\vec y)
```

```math
\int_{\partial \Omega} \frac{\partial \psi}{\partial n_x} \frac{1}{|\vec x - \vec y} dS_{\vec x} = \frac{2 \mu}{\mu + 1} \int_{\partial \Omega} \psi(\vec x) \frac{\partial}{\partial n_{x}} \frac{1}{|\vec x - \vec y|} dS_{\vec x}
```

That in turn allows to obtain the singular equation for the potential:
```math
  \psi(\vec y) = 2 (\psi_{\infty}(\vec y) + \psi_s(\vec y)) + \frac{\mu - 1}{\mu + 1} \int_{\partial \Omega} \psi(\vec x) \frac{\partial}{\partial n_x} \frac{1}{|\vec x - \vec y|} dS_{\vec x}
```

Another form of Laplace equation in boundary integral form is as follows:
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
Since:
```math
  \int_{\partial \infty} = 4 \pi \nabla \psi_{\infty} (\vec y) 
```

```math
\int_{\partial s} = 4 \pi \nabla \psi_{s} (\vec y)
```

```math
\vec n_{\vec y} \cdot \int \frac{\partial \psi}{\partial n} \nabla_{\vec x} \frac{1}{|\vec{x} - \vec{y}|} dS_{\vec x} = \frac{\mu + 1}{2} \vec n_{\vec y} \int (\nabla \psi \times \vec n_{\vec x}) \nabla_{\vec x} \frac{1}{|\vec{x} - \vec{y}|} dS_{\vec x}
```

We can evaluate normal derivative on the surface of the object as follows:
```math
  \frac{\partial \psi}{\partial n_{\vec y}} = 2 \left( \frac{\partial \psi_{\infty}}{\partial n_{\vec y}} +  \frac{\partial \psi_{s}}{\partial n_{\vec y}} \right) - \frac{\mu -1}{4 \pi} \vec n_{\vec y} \cdot \int (\nabla \psi \times \vec n_{\vec x}) \times \nabla_{\vec x} \frac{1}{|\vec{x} - \vec{y}|} dS_{\vec x}
```

\end{document}
