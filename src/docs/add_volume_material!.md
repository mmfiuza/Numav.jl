Assigns acoustic material properties to a volumetric region of the mesh, identified by its physical group tag.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | The simulation instance. |
| **Keyword arguments** | | |
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID (or vector of IDs) from the mesh. |
| `density` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Density in kg/m³. |
| `speed_of_sound` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Speed of sound in m/s. |

---
# Examples

> ```julia
> rho(f) = 1.20 # air density in kg/m³
> c(f) = 343 # speed of sound in m/s
> add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)
> ```

> Physical quantities can by complex to model sound absorption:
> ```julia
> rho(f) = 1.20 + 0.001*f
> c(f) = 340 + 0.1*sqrt(f)
> add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)
> ```

> Multiple physical groups can be assigned at once:
> ```julia
> add_volume_material!(s, physical_group=[4,6], density=1.2, speed_of_sound=343)
> ```
