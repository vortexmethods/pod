#/*--------------------------------*- POD -*------------------*---------------*\
#|    #####    ####   #####      |                            | Version 1.0    |
#|    ##  ##  ##  ##  ##  ##     |  POD: Proper Orthogonal    | 2022/08/01     |
#|    #####   ##  ##  ##  ##     |  Decomposition method      *----------------*
#|    ##      ##  ##  ##  ##     |  Open Source Code                           |
#|    ##       ####   #####      |  https://www.github.com/vortexmethods/pod   |
#|                                                                             |
#| Copyright (C) 2020-2022 Ilia Marchevsky, Irina Soldatova, Kseniia Sokol     |
#*-----------------------------------------------------------------------------*
#| File name: CMakeLists.txt                                                   |
#| Info: Source code of POD                                                    |
#|                                                                             |
#| This file is part of POD.                                                   |
#| POD is free software: you can redistribute it and/or modify it              |
#| under the terms of the GNU General Public License as published by           |
#| the Free Software Foundation, either version 3 of the License, or           |
#| (at your option) any later version.                                         |
#|                                                                             |
#| POD is distributed in the hope that it will be useful, but WITHOUT          |
#| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
#| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
#| for more details.                                                           |
#|                                                                             |
#| You should have received a copy of the GNU General Public License           |
#| along with POD.  If not, see <http://www.gnu.org/licenses/>.                |
#\*---------------------------------------------------------------------------*/

cmake_minimum_required(VERSION 2.8)
project(POD LANGUAGES CXX)

set(CMAKE_CONFIGURATION_TYPE "Debug;Release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(SRC_LIST src/pod.cpp)

FILE(GLOB ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cuh
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cuh
  )

add_executable(${PROJECT_NAME} ${ALLFILES})

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()


SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${dirlist})
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
include_directories(${dirlist})


#OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})