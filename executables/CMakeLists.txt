cmake_minimum_required (VERSION 3.10)

add_executable (LinearProgramming main.cpp)

target_link_libraries(LinearProgramming PRIVATE
    core)

target_include_directories(LinearProgramming PRIVATE
    ${ROOT_DIRECTORY}/include
    ${ROOT_DIRECTORY}/utils
    ${CLP_INCLUDE})

file(COPY ${ROOT_DIRECTORY}/data/mps/
    DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/data/mps)

file(COPY ${ROOT_DIRECTORY}/data/feasible_solutions/
    DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/data/feasible_solutions)