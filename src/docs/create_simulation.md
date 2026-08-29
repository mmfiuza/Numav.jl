Creates and returns an instance of the `Simulation` type.

| Keyword arguments | Type | Supported options | Description |
|:--|:--|:--|:--|
| `numerical_method` | `Numav.Option` | `Fem` | Numerical method used. |
| `equation` | `Numav.Option`  | `Helmholtz` | Differential equation to be solved. |
| `element_shape` | `Numav.Option`  | `Tetrahedron` | Geometrical shape of elements. |
| `element_order` | `Numav.Option`  | `Linear`, `Quadratic` | Polynomial order of the finite elements. |

---
# Examples

> ```julia
> s = create_simulation(
>     numerical_method = Fem,
>     equation = Helmholtz,
>     element_shape = Tetrahedron,
>     element_order = Linear
> )
> ```
