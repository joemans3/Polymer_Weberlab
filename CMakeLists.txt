#CMake Entry Point
cmake_minimum_required(VERSION 2.6)
project(Polymer)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
subdirs(src)