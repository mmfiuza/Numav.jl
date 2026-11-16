# Contributing

This is a guide to expose the overall structure of the Numav project. It covers the most important aspects worth knowing if you want to contribute to the development.

## Docker development container

To contribute to Numav, the first step is to use the Docker development container. This ensures a development environment that is consistent across multiple machines. It can be started in two ways:
- directly in a Linux terminal environment;
- through the Dev Containers VS Code extension.

### Via Linux terminal
To start the container directly in a Linux terminal environment, go to the root directory of this repository and run:
```
./docker/start-container.sh <path-to-your-github-ssh-keys>
```
`<path-to-your-github-ssh-keys>` is an optional argument. Its default value is `~/.ssh`.

### Via VS Code
To start the container through the Dev Containers VS Code extension, install it and open the root folder of this repository. Then, click the remote indicator in the bottom-left corner and select **Reopen in Container**.

## Interaction between Julia and C++

Numav is a Julia library mostly written in C++. The communication between the Julia and C++ parts is done via CxxWrap, which is a pair of two libraries: `libcxxwrap_julia` and `CxxWrap.jl`. The first one is a C++ library that is able to convert elements of C++ libraries to symbols in a shared library that can be accessed through Julia. On the other hand, `CxxWrap.jl` is a Julia library that can read the symbols in a shared library created with `libcxxwrap_julia`. Thus, `libcxxwrap_julia` and `CxxWrap.jl` are libraries from different languages that are able to exchange information. This process makes it possible to map C++ functions to Julia functions and vice versa, making a two-way bind.

Thus, the complete Numav ecosystem is split into four libraries:
- `libnumav`;
- `libnumav_julia`;
- `numav_julia_jll`;
- `Numav.jl`.

### 1: `libnumav`

`libnumav` is a pure C++ library, having the core calculation engine of Numav ecosystem. It is located in the `c++` directory of this repository. It is not recommended to be used as a standalone C++ library. Its API was never intended to be user-friendly or well documented. The only purpose of `libnumav` is to have a tiny amount of functions that can be called from Julia later.

### 2: `libnumav_julia`

`libnumav_julia` is a C++ library that provides the functions of `libnumav` in a way that is accessible from Julia. More specifically, `libcxxwrap_julia` is used with `libnumav`, resulting in `libnumav_julia`.

### 3: `numav_julia_jll`

`numav_julia_jll` is a Julia library that contains `libnumav_julia` compiled as a dynamic library. More specifically, it is considered a [JLL package](https://docs.binarybuilder.org/stable/jll/), registered as an official Julia package downloadable through the package manager ([Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/)). The main objective of this package is to create an easy way for the Julia users to download a precompiled copy of `libnumav_julia`.

### 4: `Numav.jl`

`Numav.jl` is the top-level Julia library, that the user will manually install. It has `numav_julia_jll` as a dependency that is automatically installed, bringing Numav's C++ component. The API of Numav is meant to be very user-friendly and well documented. It exploits the flexibility of Julia while converting instructions to the C++ component, which is not flexible, but has great performance.

## Building `libnumav`

If needed, `libnumav` can be built as a static or dynamic standalone C++ library by the following commands:

- **Static build**
```
cd c++ &&
rm -rf build &&
cmake . -B build -D CMAKE_BUILD_TYPE=Release \
-D CMAKE_PREFIX_PATH="\
/opt/intel/oneapi/mkl/2025.2/lib/cmake" &&
cmake --build build --parallel $(nproc) &&
cd ..
```

- **Dynamic build**
```
cd c++ &&
rm -rf build &&
cmake . -B build -D CMAKE_BUILD_TYPE=Release \
-D BUILD_SHARED_LIBS=TRUE \
-D CMAKE_PREFIX_PATH="\
/opt/intel/oneapi/mkl/2025.2/lib/cmake" &&
cmake --build build --parallel $(nproc) &&
cd ..
```

## Developing the package locally

In Julia's package manager ([Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/)), `dev` (short for `develop`) tells Julia to use a package directly from a folder on your machine, instead of using the registered version that Pkg downloads and stores as a read-only copy. For example, `dev some_path` in the Pkg mode of the Julia REPL makes your current environment load the package found in `some_path`. Thus, any edit you make to the source code is used the next time the package is loaded, without having to reinstall anything. Also, `dev` can be reverted with the command `free`, like `dev some_package`.

`Numav.jl` is meant to be developed inside the Docker development container described above. What you need to set up for local development depends on what you intend to change:

- **Option 1 (for Julia part only):** you only need to `dev` `Numav.jl` itself. The C++ component is taken from the precompiled `numav_julia_jll` that is already registered, so nothing has to be compiled.
- **Option 2 (for C++ part):** you need to build your own local copy of `numav_julia_jll` and `libnumav_julia`, and make `Numav.jl` use them instead of the registered ones.

### Versions of `numav_julia_jll` and `Numav.jl`

`numav_julia_jll` and `Numav.jl` are two different packages, and **they have different version numbers**. `Numav.jl` follows the versioning of the Julia API that users interact with, while `numav_julia_jll` follows the versioning of the compiled binaries. A new version of one of them does not imply a new version of the other: a change that only touches the Julia code does not require a new JLL, and a change in the C++ code requires a new JLL and a matching update of the version required by `Numav.jl`.

### Yggdrasil and the JLL recipe

`numav_julia_jll` cannot be built by hand. It has to be built by its recipe in [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil), the community repository of build recipes for JLL packages. The recipe of Numav is the file `N/numav_julia/build_tarballs.jl` of [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil).

For that reason, if you plan to change the C++ part, it is **highly recommended to clone Yggdrasil to the root of this repository**:
```
git clone https://github.com/JuliaPackaging/Yggdrasil.git
```
The commands in this guide assume that the recipe is located at `Yggdrasil/N/numav_julia/build_tarballs.jl` relative to the root of the repository (that is, `/workspace/Yggdrasil/...` inside the container). Yggdrasil is not part of Numav, so do not commit it to this repository.

### Option 1: changing only the Julia part

Open Julia in the root of the repository and, in the package manager (press `]` in the REPL), run:
```
dev /workspace
```
This replaces the registered version of `Numav.jl` by the local copy, so every change you make to the Julia code is used directly.

### Option 2: changing the C++ part too

Follow these steps in order.

**1. Build the JLL for local testing.** From the root of the repository, run the recipe in Yggdrasil with the `--deploy-jll=local` option:
```
rm -rf build && rm -rf install && rm -rf products && julia +1.12.4 \
/workspace/Yggdrasil/N/numav_julia/build_tarballs.jl \
--deploy-jll=local x86_64-linux-gnu-julia_version+1.12.0
```
This builds `numav_julia_jll` and deploys it locally, without registering or publishing anything. The package is placed in `/usr/local/share/julia/dev/numav_julia_jll`.

**2. Dev the Julia packages.** In the package manager, `dev` both the local JLL and `Numav.jl`:
```
dev /usr/local/share/julia/dev/numav_julia_jll
dev /workspace
```

**3. Generate the override for `numav_julia_jll`.** The override makes the JLL load the `libnumav_julia` library that you compile locally, instead of the precompiled one:
```
julia c++/generate_override.jl
```

**4. Build `libnumav_julia` for local testing.** This compiles the C++ code with the Julia bindings and installs it into the override directory:
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
/usr/local/share/julia/artifacts/6d260b2393efd5030d726537e5efe3573d0fbd28" &&
cmake --build build --parallel $(nproc) &&
cmake --install ./build
```
The second entry of `CMAKE_PREFIX_PATH` points to a Julia artifact directory. If the hash in that path does not exist on your machine, adjust it to the artifact directory that was installed in your environment.

From now on, whenever you change the C++ code, repeat step 4 to rebuild and reinstall `libnumav_julia`. Steps 1 to 3 only need to be repeated if the recipe or the JLL itself changes.

## Documentation

The documentation is built with the following command:
```
julia --project=docs docs/make.jl
```

To preview it locally, run:
```
julia docs/preview.jl
```
and open http://localhost:8000/ in you browser.

## Testing

In the `test` directory, Numav has a test suite. It is meant to ensure that every functionality works as intended, providing accurate results.

To run the tests, type in the Pkg mode of the Julia REPL (press `]`):
```julia
test Numav
```
This will test the instance of Numav that is currently loaded by Pkg, whether it is a `dev` copy or the registered version.

New features and bug fixes should come with tests.

## Releasing new versions

`Numav.jl` and `numav_julia_jll` have different version numbers and are released separately, each through its own process. If a change touches the C++ part, release the JLL first, because `Numav.jl` needs to depend on the new JLL version.

### Releasing `Numav.jl`

1. Update the `version` field in the `Project.toml` of `Numav.jl`. If the release depends on a new version of `numav_julia_jll`, also update its entry in the `[compat]` section.
2. Make sure the commit you want to release is merged into `main`.
3. On that commit in GitHub, write the following comment:
```
@JuliaRegistrator register branch=main
```
The comment will trigger [@JuliaRegistrator](https://github.com/JuliaRegistrator) to open a pull request to the Julia General registry with the new version. Once it is merged, the new version becomes available to users through the package manager. In this process the documentation of the new version will be automatically built and deployed to the GitGub Pages. Also, a new release and tag will be created automatically on Numav's repository.

### Releasing `numav_julia_jll`

`numav_julia_jll` is built and registered from its recipe in [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil), so a new version is released by opening a pull request to that repository:

1. Edit the recipe at `N/numav_julia/build_tarballs.jl`. Update the `version` to the new JLL version, and update the source code URL so that it points to the code you want to build.
2. Open a pull request with the changes to Yggdrasil.
3. When the pull request is merged, Yggdrasil's automated build process builds the binaries for all supported platforms and registers the new `numav_julia_jll` version in the General registry.

Remember that the version number of the JLL is independent of the version number of `Numav.jl`, so do not try to keep them equal.
