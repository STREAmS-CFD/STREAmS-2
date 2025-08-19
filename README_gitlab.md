Welcome to STREAmS-2 repository!

![Pipeline Status](https://codehub.hlrs.de/coes/excellerat-p2/uc-6/streams/badges/main/pipeline.svg)

STREAmS performs Direct Numerical Simulations of compressible turbulent flows solving the unsteady, fully compressible Navier-Stokes equations for a perfect gas. 

Currently, three canonical wall-bounded flows can be simulated in Cartesian geometry:

* compressible turbulent channel flow
* compressible zero-pressure-gradient turbulent boundary layer
* supersonic oblique shock-wave/turbulent boundary-layer interaction

Besides, the canonical flows can be simulated for 2D (x-y) curvilinear grids:

* compressible curved channel flow,
* external compressible flow on C-meshes (e.g., flow over an airfoil),
* spatially developing boundary layer on curved walls (e.g., flow over compression corners).

STREAmS can be used on both local machines and on massively parallel HPC architectures, including those based on Graphical Processing Units (GPUs).
In particular, pure MPI, CUDA Fortran, OpenMP, OpenMP-offload and HIP backends are supported.

To start using it please go to the documentation page <https://STREAmS-CFD.github.io/STREAmS-2>.

# References

Bernardini, M., Modesti, D., Salvadore, F., & Pirozzoli, S. (2021). STREAmS: A high-fidelity accelerated solver for direct numerical simulation of compressible turbulent flows. Computer Physics Communications, 263, 107906. https://doi.org/10.1016/j.cpc.2021.107906

Bernardini, M., Modesti, D., Salvadore, F., Sathyanarayana, S., Della Posta, G., & Pirozzoli, S. (2023). STREAmS-2.0: Supersonic turbulent accelerated Navier-Stokes solver version 2.0. In Computer Physics Communications, 285, 108644. https://doi.org/10.1016/j.cpc.2022.108644
