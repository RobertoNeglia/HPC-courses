cmake_minimum_required(VERSION 3.12.0)
project(lab01 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(poisson_exam src/poisson_exam.cpp src/Poisson2D.cpp)
deal_ii_setup_target(poisson_exam)

# add_executable(stokes_exam_sim src/stokes_exam_sim.cpp src/Stokes.cpp)
# deal_ii_setup_target(stokes_exam_sim)