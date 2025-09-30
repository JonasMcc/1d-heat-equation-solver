cmake_minimum_required(VERSION 3.16)
project(heat-equation-solver VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Collect all source files
file(GLOB_RECURSE SOURCES
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/main.cpp
)

# Create executable with a valid target name
add_executable(HeatEquationSolver ${SOURCES})

# Include directories scoped to the executable
target_include_directories(HeatEquationSolver PRIVATE ${PROJECT_SOURCE_DIR}/include)
