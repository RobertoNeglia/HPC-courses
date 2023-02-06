cmake_minimum_required(VERSION 3.12.0)
project(lab01 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(diffusion-advection-exam src/diffusion-advection-exam.cpp src/Poisson2D.cpp)
deal_ii_setup_target(diffusion-advection-exam)

add_executable(stokes_exam src/stokes_exam.cpp src/Stokes.cpp)
deal_ii_setup_target(stokes_exam)