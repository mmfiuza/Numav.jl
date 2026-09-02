Loads the mesh file. It defines the geometry of the domain and must contain the physical groups that are referenced when assigning boundary conditions. To generate meshes, a good open-source option is [Gmsh](https://gmsh.info/).

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | The simulation instance. |
| `path_to_mesh` | `String` | Path to the mesh file (`.bdf`). |

---
# Examples

> ```julia
> load_mesh!(s, "some_mesh_file.bdf")
> ```

!!! warning
    Currently, only meshes in the **small field** BDF format are supported.

!!! warning
    The coordinates in the mesh file must be given in **meters**.
