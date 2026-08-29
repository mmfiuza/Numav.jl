Plots a 3D graph of the pressure field for all space and frequncies.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `file_path` | `String` | Result file path (`.h5`). |
| **Keyword arguments** | | |
| `db` | `Bool` | Controls if the plot uses the decibel scale. Defaults to `false`. |
| `colormap` | `Symbol` | Colors used for the scale. Available options [here](https://docs.makie.org/dev/explanations/colors). Defaults to `:rainbow1`. |

---
# Examples

> ```julia
> plot_pressure_field("result.h5")
> ```

!!! warning
    This function is not accurate for second order elements. The graph shown does not account for the non-vertex nodes.
