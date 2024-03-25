#[[
Copyright (C) 2024 Geoffrey Daniels. https://gpdaniels.com/

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License only.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
]]

# Define the minimum version of CMake that is required.
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1)
CMAKE_POLICY(VERSION 3.5.1)

# Project name.
SET(PROJECT_NAME "pnd3d")
PROJECT(${PROJECT_NAME})

# Print status.
MESSAGE(STATUS "Configuring ${PROJECT_NAME}...")

# Require c11.
SET(CMAKE_C_STANDARD 11)
SET(CMAKE_C_STANDARD_REQUIRED ON)
SET(CMAKE_C_EXTENSIONS OFF)

# Require c++17.
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

# Set compiler build flags.
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++98-compat")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++98-compat-pedantic")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-newline-eof")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reserved-id-macro")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-padded")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-disabled-macro-expansion")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-global-constructors")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-exit-time-destructors")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-weak-vtables")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-warning-option")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unsafe-buffer-usage")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-macros")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++98-compat")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++98-compat-pedantic")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-newline-eof")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reserved-id-macro")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-padded")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-disabled-macro-expansion")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-global-constructors")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-exit-time-destructors")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-weak-vtables")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-warning-option")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unsafe-buffer-usage")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-macros")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:strict")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_HAS_EXCEPTIONS=0")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GR-")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ITERATOR_DEBUG_LEVEL=0")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_USE_MATH_DEFINES=1")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOMINMAX=1")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ENABLE_ATOMIC_ALIGNMENT_FIX=1")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CRT_NONSTDC_NO_WARNINGS=1")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CRT_SECURE_NO_WARNINGS=1")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
ELSE()
    MESSAGE(FATAL_ERROR "Unsupported compiler.")
ENDIF()

# Set compiler linker flags.
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -latomic")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    # None
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # None
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # None
ELSE("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    MESSAGE(WARNING "Unsupported compiler.")
ENDIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

# Enable threading.
SET(THREADS_PREFER_PTHREAD_FLAG ON)
FIND_PACKAGE(Threads REQUIRED)

# Create the executable.
FILE(GLOB_RECURSE SOURCE_FILES "${PROJECT_SOURCE_DIR}/source/*")
ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCE_FILES})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/source/")
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE Threads::Threads)

IF(UNIX)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE GL)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE X11)
ENDIF()

# Print status.
MESSAGE(STATUS "Finished configuring ${PROJECT_NAME}.")
