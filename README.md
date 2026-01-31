Project for the course "Implementation of the Finite Element Method"

This repository contains a Finite Element Method (FEM) implementation for solving Navier-Stokes equations with Coriolis force considerations.

Follow these steps to build the project on your local machine.
Prerequisites

    CMake (3.10 or higher)

    C++ Compiler (GCC, Clang, or MSVC)

    Make

Compilation Steps

From the project root directory, execute the following commands in your terminal:

    Create and enter the build directory:
    Bash

    mkdir build && cd build

    Configure the project with CMake:
    Bash

    cmake ..

    Prepare assets: Copy the shaders into the build folder so the executable can locate them at runtime.
    Bash

    cp -r ../shaders .

    Build the executable:
    Bash

    make

    ./test_coriolis sphere/cube -n
