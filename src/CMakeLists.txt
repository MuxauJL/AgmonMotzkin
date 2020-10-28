cmake_minimum_required (VERSION 3.10)

if (NOT DEFINED CLP_LIBS)
    message(FATAL_ERROR "Please set the variable CLP_LIBS")
endif()

if (NOT DEFINED CLP_INCLUDE)
    message(FATAL_ERROR "Please set the variable CLP_INCLUDE")
endif()

file(GLOB SRC_FILES 
    ${ROOT_DIRECTORY}/src/*.cpp)

add_library(core STATIC ${SRC_FILES})

target_link_libraries(core PRIVATE
    "${CLP_LIBS}/Clp.lib"
    "${CLP_LIBS}/CoinUtils.lib"
    "${CLP_LIBS}/Osi.lib")

target_include_directories(core PRIVATE
    ${ROOT_DIRECTORY}/include
    ${ROOT_DIRECTORY}/utils
    ${CLP_INCLUDE})