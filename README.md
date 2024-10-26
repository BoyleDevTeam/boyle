# Boyle Library: The Fundamental Math Basis for Autonomous Driving Vehicles and Robotics

Boyle Library aims to provide the most fundamental numeric implementation for autonomous driving vehicles and robotics. It contains piecewise polynomial functions/curves, sparse matrix, convex optimization modules for trajectory generations.

## Environment Setup

Prerequisite:

* Modern Linux OS (recommend Arch Linux and Debian 13, but could work on Ubuntu 22.04/24.04),
* gcc >= 14.0.0 or clang > 19.0.0,
* cmake >= 4.0.0,
* clangd >= 19.0.0,
* clang-format >= 19.0.0,
* clang-tidy >= 19.0.0,
* ccache >= 4.10.0,
* gnuplot (for graphic unit test showing)

Necessary Visual Studio Code Plugins:

* C/C++ Extension Pack,
* CMake Tools,
* Better C++ Syntax,
* clangd,
* Doxygen Document Generator,
* Better Comments,
* GitLens.

The recommended vscode settings and git configurations are included under directory:

* developers/houchen_li/personal_configs.

## Run Graphic Unit Test

1. cd to "out/build/linux-gcc-x64/\<the relative path of the unit test source code\>".
2. execute the binary file with an argument "--plot-graph".
