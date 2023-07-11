# Boyle Library: The Fundamental Math Basis for Autonomous Driving Vehicles and Robotics

Boyle Library aims to provide the most fundamental numeric implementation for autonomous driving vehicles and robotics. It contains piecewise polynomial functions/curves, sparse matrix, osqp wrapper, convex optimization modules for trajectory generations.

## Environment Setup

Prerequisite:

* Modern Linux OS (recommend Arch Linux and Debian 12, but could work on Ubuntu 20.04/22.04),
* gcc >= 12.0.0 or clang > 15.0.0,
* cmake >= 3.25,
* boost >= 1.71.0,
* clangd >= 15.0.0
* clang-format >= 15.0.0.
* clang-tidy >= 15.0.0
* gnuplot (for graphic unit test showing)

Necessary Visual Studio Code Plugins:

* C/C++ Extension Pack,
* CMake Tools,
* clangd,
* Doxygen Document Generator,
* GitLens.

The recommended vscode settings and git configurations are included under directory:

* developers/houchen_li/personal_configs.

## Run Graphic Unit Test

1. cd to "out/build/gcc-linux_x64/\<the relative path of the unit test source code\>".
2. execute the binary file with an argument "--plot-graph".
