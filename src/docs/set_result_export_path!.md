Specifies the file path where simulation results will be written when calling [`run!`](@ref).

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | The simulation instance. |
| `file_path` | `String` | Result file path (`.h5`). |

---
# Examples

> ```julia
> set_result_export_path!(s, "result.h5")
> ```
