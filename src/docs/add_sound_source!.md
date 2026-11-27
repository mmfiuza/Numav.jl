Adds sound sources either at a specific point in space (via `coordinates`) or over a surface region (via `region`). Three excitation types are supported: `volume_velocity`, `particle velocity`, and `pressure`. Sign conventions are explained [here](@ref "Conventions").

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | The simulation instance. |
| **Keyword arguments** | | |
| `coordinates` | `Vector{Real}`, `Vector{Vector{Real}}` | `[x, y, z]` location of a point source in m. |
| `region` | `Integer`, `Vector{Integer}` | Region ID of a surface or volume. |
| `volume_velocity` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Volume velocity in m³/s. |
| `particle_velocity` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Normal particle velocity in m/s. |
| `pressure` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Acoustic pressure in Pa. |

---
# Examples

> Volume velocity source (monopole):
> ```julia
> Q(f) = 10/f # Volume velocity in m³/s as a function of frequency
> add_sound_source!(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)
> ```
> Suitable for modeling punctual omnidirectional sources.

> Particle velocity source (vibrating surface):
> ```julia
> U(f) = 15/f # Particle velocity in m/s
> add_sound_source!(s, region=2, particle_velocity=U)
> ```
> Prescribes the normal component of particle velocity on all surfaces of a region. Useful for modeling vibrating panels or pistons.

> Pressure source:
> ```julia
> P(f) = 2f # Pressure in Pa as a function of frequency
> 
> # At a point in space
> add_sound_source!(s, coordinates=[2.0, 2.5, 1.0], pressure=P)
> 
> # On a mesh surface
> add_sound_source!(s, region=3, pressure=P)
> ```
> Prescribes acoustic pressure, either at a point or on a surface.

> Multiple points or regions can be assigned at once:
> ```julia
> p1 = [1.0, 3.0, 2.0]
> p2 = [3.0, 1.0, 1.0]
> add_sound_source!(s, coordinates=[p1,p2], volume_velocity=0.01)
> add_sound_source!(s, region=[3,5,9,2], particle_velocity=0.01)
> ```

!!! note
    Each call to `add_sound_source!` should specify either `coordinates` or `region` (not both), and exactly **one** excitation type keyword.
