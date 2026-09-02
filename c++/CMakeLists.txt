# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

cmake_minimum_required(VERSION 3.10)
project(numav LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_SHARED_LIBS "Build shared library" OFF)

# find dependencies
find_package(Eigen3 REQUIRED)

# solver-specific configuration
set(SOLVER "ONEMKL" CACHE STRING "The linear system solver")
set_property(CACHE SOLVER PROPERTY STRINGS "INTERNAL" "EIGEN" "ONEMKL" "MUMPS")
if(SOLVER STREQUAL "INTERNAL")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_INTERNAL_SOLVER)
    set(SOLVER_SPECIFIC_LINK_LIBRARIES ldlt_solver)
    add_subdirectory(third-party/ldlt-solver)
elseif(SOLVER STREQUAL "EIGEN")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_EIGEN_SOLVER)
    # set(SOLVER_SPECIFIC_LINK_LIBRARIES -lblas)
    # set(SOLVER_SPECIFIC_DEFINITIONS EIGEN_USE_BLAS)
elseif(SOLVER STREQUAL "ONEMKL")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_ONEMKL_SOLVER)
    set(MKL_LINK sdl) 
    find_package(MKL REQUIRED)
    set(SOLVER_SPECIFIC_LINK_LIBRARIES MKL::MKL)
    set(SOLVER_SPECIFIC_DEFINITIONS MKL_ILP64)
elseif(SOLVER STREQUAL "MUMPS")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_MUMPS_SOLVER)
    set(SOLVER_SPECIFIC_INCLUDE_DIRECTORIES /usr/include)
    set(SOLVER_SPECIFIC_LINK_LIBRARIES -lzmumps_seq)
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/lib")

# add numav
file(GLOB_RECURSE NUMAV_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE NUMAV_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_library(numav ${NUMAV_SOURCES} ${NUMAV_HEADERS})

target_include_directories(numav 
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/third-party/safe-ptr/include
        ${SOLVER_SPECIFIC_INCLUDE_DIRECTORIES}
)

target_compile_definitions(numav PUBLIC 
    FMT_HEADER_ONLY
    $<$<CONFIG:Debug>:SAFE_PTR_DEBUG> 
    EIGEN_DEFAULT_DENSE_INDEX_TYPE=int64_t
    $<$<NOT:$<CONFIG:Debug>>:EIGEN_NO_DEBUG>
    EIGEN_NO_AUTOMATIC_RESIZING
    EIGEN_FAST_MATH
    ${SELECTED_SOLVER_DEFINITION}
    ${SOLVER_SPECIFIC_DEFINITIONS}
    ${ADDITIONAL_COMPILE_DEFINITIONS}
)

target_link_libraries(numav PUBLIC
    Eigen3::Eigen
    ${SOLVER_SPECIFIC_LINK_LIBRARIES}
)

# set numav's compile and link flags
target_compile_options(numav PRIVATE
    -Wall -Wextra
    $<$<CONFIG:Debug>:-g -O0 -fno-inline>
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:RelWithDebInfo>:-O2 -g> # Release with symbols
    $<$<CONFIG:MinSizeRel>:-Os> # Minimum size
    $<$<BOOL:${BIND_JULIA}>:-fPIC>
)
target_link_options(numav PRIVATE
    $<$<NOT:$<CONFIG:Debug>>:-flto> # LTO for everything except Debug
)

if(BIND_JULIA)
    find_package(JlCxx REQUIRED)
    get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
    get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
    
    add_library(numav_julia SHARED
        ${CMAKE_CURRENT_SOURCE_DIR}/numav_julia.cpp
    )
    target_link_libraries(numav_julia
        numav
        JlCxx::cxxwrap_julia
        JlCxx::cxxwrap_julia_stl
    )
    install(TARGETS numav_julia LIBRARY DESTINATION lib)
endif()
