project(aco-tsp)

cmake_minimum_required(VERSION 3.8)

add_executable(main src/main.cpp)

target_link_libraries(main PRIVATE matplotlib_cpp)

include_directories("lib/matplotlib-cpp")
add_subdirectory("lib/matplotlib-cpp")

