cmake_minimum_required(VERSION 3.12.0)
project(lab01 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(lab-06-exercise1 src/lab-06-exercise1.cpp src/NonLinearDiffusion.cpp)
deal_ii_setup_target(lab-06-exercise1)

add_executable(lab-06-exercise2 src/lab-06-exercise2.cpp src/LinearElasticity.cpp)
deal_ii_setup_target(lab-06-exercise2)
