cmake_minimum_required(VERSION 3.12.0)
project(lab01 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(diffusion-reaction-exam src/diffusion-reaction-exam.cpp src/Poisson2D.cpp)
deal_ii_setup_target(diffusion-reaction-exam)

add_executable(heat-equation-exam src/heat-equation-exam.cpp src/Heat.cpp)
deal_ii_setup_target(heat-equation-exam)