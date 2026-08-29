Assigns a specific acoustic impedance boundary condition to a surface in the mesh. This is used to model absorbers, reflecting surfaces or other boundary treatments.

The specific acoustic impedance is the ratio of complex amplitude of acoustic pressure to complex amplitude of normal particle velocity (in Pa·s/m).

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | The simulation instance. |
| **Keyword arguments** | | |
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID of the boundary surface. |
| `specific_acoustic_impedance` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | specific surface acoustic impedance (Pa·s/m). |

---
# Examples

> ```julia
> Z(f) = 1f + 2im
> add_surface_material!(s, physical_group=4, specific_acoustic_impedance=Z)
> ```

> Multiple physical groups can be assigned at once:
> ```julia
> add_surface_material!(s, physical_group=[3,1], specific_acoustic_impedance=1.0)
> ```
