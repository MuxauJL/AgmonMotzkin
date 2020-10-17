cmake_minimum_required (VERSION 3.10)

if (NOT DEFINED CLP_LIBS)
    message(FATAL_ERROR "Please set the variable CLP_LIBS")
endif()

if (NOT DEFINED CLP_INCLUDE)
    message(FATAL_ERROR "Please set the variable CLP_INCLUDE")
endif()

file(GLOB SRC_FILES 
    ${ROOT_DIRECTORY}/src/*.cpp
    ${ROOT_DIRECTORY}/utils/*.cpp)

add_executable (LinearProgramming ${SRC_FILES})

target_link_libraries(LinearProgramming PRIVATE
    "${CLP_LIBS}/Clp.lib"
    "${CLP_LIBS}/CoinUtils.lib"
    "${CLP_LIBS}/Osi.lib")

target_include_directories(LinearProgramming PRIVATE
    ${CLP_INCLUDE}
    ${ROOT_DIRECTORY}/utils)

file(COPY ${ROOT_DIRECTORY}/data/mps/
    DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/data/mps)

file(COPY ${ROOT_DIRECTORY}/data/feasible_solutions/
    DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/data/feasible_solutions)