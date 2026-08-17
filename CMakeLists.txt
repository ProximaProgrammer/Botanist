cmake_minimum_required(VERSION 3.12)
project(absorptivity_module LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Suppress the legacy finder policy warning
cmake_policy(SET CMP0148 OLD)

# Force CMake to look directly inside your active Conda environment paths
include_directories("$ENV{CONDA_PREFIX}/include")
include_directories("$ENV{CONDA_PREFIX}/include/python3.13")
link_directories("$ENV{CONDA_PREFIX}/lib")

find_package(pybind11 REQUIRED)

pybind11_add_module(absorptivity molar_absorptivity.cpp)

set_target_properties(absorptivity PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
