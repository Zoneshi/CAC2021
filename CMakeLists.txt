cmake_minimum_required(VERSION 3.20)

project(CAC_2021)

set(CMAKE_CXX_STANDARD 17)

include_directories("/Users/zhongjiaoshi/VSCode/include")

find_package(Eigen3 REQUIRED)

add_executable(cac main.cpp)

target_link_libraries(cac Eigen3::Eigen)