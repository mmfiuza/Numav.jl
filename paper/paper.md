---
title: 'Numav.jl: an open source tool for acoustic simulation in Julia'
tags:
  - Acoustics
  - Julia
  - FEM
authors:
  - name: Matheus Machado Fiuza
    orcid: 0009-0006-4578-2845
    affiliation: 1
  - name: Paulo Henrique Mareze
    orcid: 0000-0002-5572-1935
    affiliation: 1
  - name: Pedro Paulo Calmon Paes
    orcid: 0009-0002-8127-286X
    affiliation: 1
affiliations:
 - name: Federal University of Santa Maria (UFSM), Brazil
   index: 1
   ror: 01b78mz79
date: 20 November 2026
bibliography: paper.bib

---

# Summary

Sound propagation in rooms, vehicle cabins, ducts and mechanical components can be predicted numerically before anything is built, which makes simulation a routine part of acoustic design. At low frequencies, where wave effects such as room modes dominate, the most widely used approach is the Finite Element Method (FEM). In it, the space is divided into small elements, and the governing wave equation is converted into a linear system whose solution approximates the sound pressure at every mesh node.

`Numav.jl` (*Numerical Methods for Acoustics and Vibrations*) is an open source library for the Julia programming language that performs such simulations. As of the release of this paper (version 0.4.0), `Numav.jl` implements a single module: FEM in the frequency domain, based on the Helmholtz equation, with tetrahedral elements of first or second order. A simulation is described in a short Julia script and installation is as easy as that of any other registered Julia package.

# Statement of need

Commercial packages offer mature implementations of acoustic FEM. However, their licensing cost, closed source code and limited customizability make them hard to use for teaching, for reproducible research and for experimenting with new formulations or workflows. Open alternatives exist, but they tend to sit at one of two extremes: general-purpose finite element frameworks that are not user-friendly for acoustics, or domain-specific tools aimed at particular applications.

`Numav.jl` targets acoustic engineers, researchers and students who want a transparent, scriptable and acoustics-specific FEM tool. The user thinks in terms of acoustic concepts (materials, surface impedances, monopoles, vibrating surfaces) rather than finite element abstractions. The code is open enough to be read, modified and extended.

Implementing FEM simulations in acoustics from scratch for the first time can be very challenging, primarily because the equations and their derivations are complicated, making it difficult for beginners to understand every step that should be coded to perform a complete simulation. Thus, one goal of this project was to provide a didactic guide explaining all the derivations involved in creating FEM itself. At the time of release of this article, Matheus Fiuza's bachelor's thesis is being written with that purpose. It will start from the Helmholtz equation and from there derive the FEM linear system without skipping mathematical steps. This thesis will end up being the summary of years of research done by the members of the Acoustic Engineering program at the [Federal University of Santa Maria (UFSM)](https://www.ufsm.br/).

# State of the field

COMSOL Multiphysics [@COMSOL] and ANSYS [@ANSYS] are the commercial reference implementations of acoustic FEM. Among open source tools, the most traditional alternative is the FEniCSx project [@Baratta:2023; @Scroggs:2022a; @Scroggs:2022b; @Alnaes:2014]. It is a FEM solver with Helmholtz equation support. However, FEniCSx is general-purpose, not prioritizing acoustics. During the development of `Numav.jl`, VIBRA [@vibra_software] was created. Coincidentally, both projects can apply FEM based on the Helmholtz equation. Crucial differences are that VIBRA is developed in Python and has a GUI.

`Numav.jl` was built as a separate project rather than as a contribution to these tools because its goals differ. First, it is acoustics-specific, so the boundary conditions and sources most used in practice are first-class concepts of the interface. Second, it is organized around frequency steps and designed to avoid recomputing anything that does not depend on frequency. Third, it is written for the Julia ecosystem, where scientific users can combine it with existing packages for data handling, optimization and plotting.

# Software design

`Numav.jl` is a Julia package whose computational core is written in C++. The Julia layer provides the user-facing API, while the C++ layer handles the performance-critical work: mesh reading, assembly of the finite element matrices, application of boundary conditions and solution of the linear systems. The C++ code is exposed to Julia through CxxWrap, compiled for multiple platforms with BinaryBuilder and delivered as a Julia binary package. This allows installation on Linux, Windows, macOS and FreeBSD without a local C++ toolchain. In practice, this architecture does not affect installation or usage. The user may not even be aware that most of the library was written in C++.

In a simulation script, the user loads a mesh, assigns acoustic materials and sound sources, chooses the frequencies to simulate and runs the solver. Meshes should be generated externally by the user and loaded in Nastran's `.bdf` format. Gmsh [@Geuzaine:2009] is an open source suggestion for this process. Regions defined in the mesh are used to assign materials, impedances and sources. The same interface supports point volume velocity sources (monopoles), prescribed normal particle velocity on surfaces and prescribed pressure at points or on surfaces. Density, speed of sound, impedance and source amplitudes can be constants or functions of frequency, and may be complex. This allows porous absorbers to be treated as equivalent fluids. Frequencies can be given explicitly or generated automatically. By default, the sampling density grows quadratically with frequency, following the growth of the modal density of rooms.

The formulation starts from the Helmholtz equation derived from the linear wave equation and applies the Galerkin method with tetrahedral elements of first or second order, using the subparametric formulation. Thus, second order elements with curves are not supported. The input mesh is expected to be first order, so that `Numav.jl` automatically generates the edge nodes.

Each frequency step leads to an independent linear system, which is complex, sparse and symmetric once assembled. Three design decisions concentrate on performance:

- **Reuse across frequencies.** The sparsity pattern of the system matrix depends only on the mesh, so it is computed once and used for all frequencies. This reduces the amount of logical operations and memory usage.
- **Sparse storage.** Only the upper triangle of the symmetric matrix is stored in compressed sparse row format.
- **Symmetric direct solver.** Each system is solved with a sparse $LDL^T$ factorization with pivoting. On x86-64 Linux and Windows, `Numav.jl` uses the Intel oneMKL Direct Sparse Solver [@Intel:oneMKL]. On the remaining platforms, it falls back to Eigen [@Guennebaud:Eigen], which performs worse.

Results are stored in [HDF5](https://www.hdfgroup.org/solutions/hdf5/) files containing the complex pressure amplitude at every node and frequency, which can be post-processed in Julia or other tools. For post processing in Julia, [HDF5.jl](https://juliaio.github.io/HDF5.jl/) is recommended. Result visualization is also provided by `Numav.jl` itself, in graphs made with Makie [@Danisch:2021].

Although `Numav.jl` currently implements a single module, it is designed for extension. The modular C++ core and API are intended to accommodate additional numerical methods for acoustics and vibrations, which the name of the project reflects.

The library's visual identity is on Figure \ref{fig:visual}. The logo was created with [`Luxor.jl`](https://juliagraphics.github.io/LuxorManual) and the documentation was created with [`Documenter.jl`](https://documenter.juliadocs.org).

<!-- \begin{figure}
  \centering
  \begin{minipage}{0.3\textwidth}
    \centering
    \includegraphics[height=0.25\paperwidth]{logo.pdf}\\
    \textbf{(a) Official logo.}
  \end{minipage}
  \hfill
  \begin{minipage}{0.69\textwidth}
    \centering
    \includegraphics[height=0.25\paperwidth]{docs.png}\\
    \textbf{(b) Documentation page.}
  \end{minipage}
  \caption{Numav.jl's visual identity.}
  \label{fig:visual}
\end{figure} -->

\begin{figure}
  \centering
  \begin{minipage}{0.3\textwidth}
    \centering
    \begin{minipage}[c][0.25\paperwidth][c]{\linewidth}
      \centering
      \includegraphics[width=0.8\linewidth,height=0.25\paperwidth,keepaspectratio]{logo.pdf}
    \end{minipage}\par\medskip
    \textbf{(a) Official logo.}
  \end{minipage}
  \hfill
  \begin{minipage}{0.69\textwidth}
    \centering
    \begin{minipage}[c][0.25\paperwidth][c]{\linewidth}
      \centering
      \includegraphics[width=\linewidth,height=0.25\paperwidth,keepaspectratio]{docs.png}
    \end{minipage}\par\medskip
    \textbf{(b) Documentation page.}
  \end{minipage}
  \caption{Numav.jl's visual identity.}
  \label{fig:visual}
\end{figure}

# Research impact statement

`Numav.jl` was developed at the Acoustic Engineering program of UFSM and is released under the [AGPL-3.0 license](https://www.gnu.org/licenses/agpl-3.0.en.html), with documentation at [`https://mmfiuza.github.io/Numav.jl`](https://mmfiuza.github.io/Numav.jl) and source code at [https://github.com/mmfiuza/Numav.jl](https://github.com/mmfiuza/Numav.jl). The package is registered in the Julia package ecosystem. The repository contains example scripts and a test suite.

The test suite checks accuracy by comparing results of `Numav.jl` and COMSOL. The same acoustic problem is defined in both programs, using the exact same underlying numerical method. Thus, results should agree to within numerical precision. The meshes used in these tests are extremely unrefined, so that the tests run fast, but the results are not physically accurate.

For that reason, a comparison was made specifically for this article, using a room with realistic dimensions and mesh refinement, as shown in Figure \ref{fig:room}. Multiple conditions were applied simultaneously, including:

- surface velocity;
- specific acoustic impedance;
- volume velocity at point;
- pressure at point;
- absorber defined as an equivalent fluid.

\begin{figure}
  \centering
  \begin{minipage}{0.495\textwidth}
    \centering
    \includegraphics[width=0.8\linewidth]{room1.png}\\
    \textbf{(a) View 1.}
  \end{minipage}
  \hfill
  \begin{minipage}{0.495\textwidth}
    \centering
    \includegraphics[width=0.8\linewidth]{room2.png}\\
    \textbf{(b) View 2.}
  \end{minipage}
  \caption{Simulated room. Available in Numav.jl's repository under the name "example1.bdf". \label{fig:room}}
\end{figure}

The simulation was set according to the following code (valid in version 0.4.0):

```julia
using Numav

# Create the simulation object
s = create_simulation(
    numerical_method = Fem,
    equation = Helmholtz,
    element_shape = Tetrahedron,
    element_order = Linear # or Quadratic
)

# Set simulated frequencies
set_frequency!(s, min=20, max=200, step=1) # (Hz)

# Load mesh
load_mesh!(s, "example1.bdf")

# Define air
rho0 = 1.2
c0= 343
add_volume_material!(s, region=1, density=rho0, speed_of_sound=c0)

# Define material as an equivalent fluid
rho1 = 1 - 1im
c1 = 400 + 400im
add_volume_material!(s, region=2, density=rho1, speed_of_sound=c1)

# Add a surface impedance boundary condition
Zs = 4000 + 4000im
add_surface_material!(s, region=11, specific_acoustic_impedance=Zs)

# Add a monopole volume velocity source at a point
Q(f) = 0.01/f
add_sound_source!(s, coordinates=[0.0, 5.0, 0.0], volume_velocity=Q)

# Add a particle velocity source on a surface
U(f) = 0.015/f
add_sound_source!(s, region=10, particle_velocity=U)

# Add pressure source
P(f) = 0.002f
add_sound_source!(s, coordinates=[4.0, 5.0, 0.0], pressure=P)

# Set output file and run
set_result_export_path!(s, "result.h5")
run!(s)
```

The sound pressure level (SPL) results at a chosen point are shown in Figure \ref{fig:graph}. As expected, the results are practically identical, confirming `Numav.jl`'s accuracy.

\begin{figure}
  \centering
  \begin{minipage}{0.495\textwidth}
    \centering
    \includegraphics[width=\linewidth]{graph-linear.pdf}
    \textbf{(a) Linear.}
  \end{minipage}
  \hfill
  \begin{minipage}{0.495\textwidth}
    \centering
    \includegraphics[width=\linewidth]{graph-quadratic.pdf}
    \textbf{(b) Quadratic.}
  \end{minipage}
  \caption{Simulated sound pressure level at some point of the room. \label{fig:graph}}
\end{figure}

Currently, `Numav.jl` is being used at the Univertity in the course of [Numerical Methods in Acoustics and Vibrations (EAC1061)](https://www.ufsm.br/ementario/disciplinas/EAC1061). Hopefully, more students will also continue to develop the project in the next years.

# AI usage disclosure

Claude Chat with Sonnet 5.5 was the most used AI service for developing the project. DeepSeek was also used for simpler tasks. These services helped with almost all the aspects of the development.

# Acknowledgements

Acknowledgements here

# References