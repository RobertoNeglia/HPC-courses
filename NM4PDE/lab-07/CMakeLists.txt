cmake_minimum_required(VERSION 3.12.0)
project(lab01 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(lab-07 src/lab-07.cpp src/Heat.cpp)
deal_ii_setup_target(lab-07)
