Assembles and solves the system of equations for all frequencies in the defined range while writing the results to the specified export path during the solution. The exported file contains the [complex pressure amplitude](@ref "Complex amplitudes") of all nodes for all frequencies.

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

Results are exported as files in the [HDF5 format](https://www.hdfgroup.org/solutions/hdf5/) (`.h5`).

!!! tip

    With [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/), you can post-process results in Julia like:
    ```julia
    using HDF5
    r = h5open("result.h5", "r")
    # inspect contents
    ```
    Alternatively, [HDFView](https://www.hdfgroup.org/download-hdfview/) can be used to check the results in a graphical interface.
