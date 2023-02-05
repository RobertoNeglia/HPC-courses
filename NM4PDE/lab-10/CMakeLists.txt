cmake_minimum_required(VERSION 3.12.0)
project(lab01 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(lab-10 src/lab-10.cpp src/Poisson2D.cpp)
deal_ii_setup_target(lab-10)