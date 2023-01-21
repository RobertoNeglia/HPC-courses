cmake_minimum_required(VERSION 3.12.0)
project(lab01 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(exam_sim src/exam_sim.cpp src/Poisson2D.cpp)
deal_ii_setup_target(exam_sim)