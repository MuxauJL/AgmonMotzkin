﻿cmake_minimum_required (VERSION 3.10)

project (LinearProgramming LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(ROOT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS YES)
endif()

add_subdirectory(src)
add_subdirectory(executables)
