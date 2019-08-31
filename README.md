# LaplaceBIE.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://akels.github.io/LaplaceBIE.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://akels.github.io/LaplaceBIE.jl/dev)
[![Build Status](https://travis-ci.org/akels/LaplaceBIE.jl.svg?branch=master)](https://travis-ci.org/akels/LaplaceBIE.jl)

Since ancient times people had wondered what makes amber rubbed with fur to attract small light objects and what makes lodestones, naturally magnetized pieces of the mineral magnetite, to attract iron. Nowadays we do have an answer that the force comes from magnetization or polarization gradient. However, even for linear materials, the computation of the force is challenging as the body itself perturbs the field. Fortunately, boundary integral methods can save our day, which are implemented in this library and can be used to calculate the field and so also a force at an arbitrary object's surface.

> Erdmanis, J. & Kitenbergs, G. & Perzynski, R. & Cebers, A. (2017) Magnetic micro-droplet in rotating field: numerical simulation and comparison with experiment
