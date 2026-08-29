# How to start the development container
The development container can be started in two ways:
- directly in a linux terminal environment;
- through the Dev Containers VS Code extension.

## Linux terminal
To start the container directly in a Linux terminal environment, go to the root directory of this repository and run:
```
./docker/start-container.sh <path-to-your-github-ssh-keys>
```
`<path-to-your-github-ssh-keys>` is an optional argument. Its default value is `~/.ssh`.

## VS Code
To start the container through the Dev Containers VS Code extension, install it and open the root folder of this repository. Then, click the remote indicator in the bottom-left corner and select **Reopen in Container**.

# How to build libnumav (static)
Run:
```
rm -rf build &&
cmake c++ -B build -D CMAKE_BUILD_TYPE=Release \
-D CMAKE_PREFIX_PATH="\
/opt/intel/oneapi/mkl/2025.2/lib/cmake;\
/HDF_Group/HDF5/2.1.1" &&
cmake --build build --parallel ${nproc}
```

# How to build libnumav (dynamic)
Run:
```
rm -rf build &&
cmake c++ -B build -D CMAKE_BUILD_TYPE=Release \
-D BUILD_SHARED_LIBS=TRUE \
-D CMAKE_PREFIX_PATH="\
/opt/intel/oneapi/mkl/2025.2/lib/cmake;\
/HDF_Group/HDF5/2.1.1" &&
cmake --build build --parallel ${nproc}
```

# How to build the JLL for local testing
```
rm -rf build && rm -rf install && rm -rf products && julia +1.12.4 \
/workspace/Yggdrasil/N/numav_julia/build_tarballs.jl \
--deploy-jll=local x86_64-linux-gnu-julia_version+1.12.0
```

# Dev the Julia packages
```
dev /usr/local/share/julia/dev/numav_julia_jll
dev /workspace
```

# Generate override for numav_julia_jll
```
julia c++/generate_override.jl
```

# How to build numav_julia for local testing
```
rm -rf build &&
rm -rf c++/override/lib &&
cmake c++ -B build \
-D CMAKE_BUILD_TYPE=Release \
-D SOLVER=ONEMKL \
-D BIND_JULIA=TRUE \
-D CMAKE_INSTALL_PREFIX=c++/override \
-D CMAKE_PREFIX_PATH="\
/opt/intel/oneapi/mkl/2025.2/lib/cmake;\
/usr/local/share/julia/artifacts/6d260b2393efd5030d726537e5efe3573d0fbd28;\
/HDF_Group/HDF5/2.1.1" &&
cmake --build build --parallel ${nproc} &&
cmake --install ./build
```

# How to build docs
```
julia --project=docs docs/make.jl
```
