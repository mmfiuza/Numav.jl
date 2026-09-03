Assembles and solves the system of equations for all frequencies in the defined range while writing the results to the specified export path during the solution.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | The simulation instance. |

---
# Examples

> ```julia
> run!(s)
> ```

---
# Output format

Results are exported as files in the [HDF5 format](https://www.hdfgroup.org/solutions/hdf5/) (`.h5`). To read it, it is recommended to use [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/) or [HDFView](https://www.hdfgroup.org/download-hdfview/).

!!! tip

    With [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/), you can post-process results in Julia like:
    ```julia
    using HDF5
    r = h5open("result.h5", "r")
    # inspect contents
    ```
