# PROJECT CONFIGURATION
cmake_minimum_required(VERSION 3.1)
# set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")

project(RotPriorRPp LANGUAGES C CXX VERSION 1.0.0)
set(LIBRARY_TARGET_NAME "RotPriorRPp")
set(LIBRARY_TARGET_NAME_EXPORT "${LIBRARY_TARGET_NAME}Export")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # We require C++ 14


# build the examples
set(BUILD_${LIBRARY_TARGET_NAME}_EXAMPLE ON)


# Build type

# Directory for built libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place libraries built by this project")
# Directory for built executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")

# BUILD CONFIGURATIONS
option(CMAKE_VERBOSE_MAKEFILE "Generate verbose makefiles?" OFF)

set(CODE_PROFILING OFF CACHE BOOL "Turn on code profiling?")
if(${CODE_PROFILING})
message(STATUS "Turning on code profiling for Essential Matrix Estimation")
endif()

# Add the .cmake files that ship with Eigen3 to the CMake module path (useful for finding other stuff)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" CACHE STRING "The CMake module path used for this project")

# FIND EIGEN3
set( ENV{EIGEN3_ROOT_DIR} ${CMAKE_SOURCE_DIR}/eigen)
find_package(Eigen3 3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION_STRING})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()


ADD_DEFINITIONS ( -DPQXX_HIDE_EXP_OPTIONAL )
#  do not include this if add_subdirectory
#  for essential is used
add_subdirectory(dependences/Rosen_optimization)

# add_subdirectory(dependences/iterativeCertifier)
find_package(IterCertAlg REQUIRED)

# IMPORT OpenGV library
# add_definitions(-march=native)



set(${LIBRARY_TARGET_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(${LIBRARY_TARGET_NAME}_SOURCE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(${LIBRARY_TARGET_NAME}_EXAMPLES_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/examples)


# Expose the include directories for this project
set(${LIBRARY_TARGET_NAME}_ADD_INCLUDES ${ESSENTIAL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} )
set(${LIBRARY_TARGET_NAME}_CERT_INCLUDES ${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS} ${${LIBRARY_TARGET_NAME}_ADD_INCLUDES})

message(STATUS "Essential headers in ${${LIBRARY_TARGET_NAME}_CERT_INCLUDES} \n")


# Get the set of Essential header and source files
set(${LIBRARY_TARGET_NAME}_HDRS
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/RotPrior.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/RotPriorManifold.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/RotPriorProblem.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/RotPriorTypes.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/RotPriorUtils.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/RotPriorConstraints.h
)


set(${LIBRARY_TARGET_NAME}_SRCS
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/RotPrior.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/RotPriorManifold.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/RotPriorProblem.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/RotPriorUtils.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/RotPriorConstraints.cpp
)




message(STATUS "Building in ${CMAKE_BUILD_TYPE} mode\n")

# Build the Essential library
add_library(${LIBRARY_TARGET_NAME} ${${LIBRARY_TARGET_NAME}_HDRS} ${${LIBRARY_TARGET_NAME}_SRCS}  )

# Set two minimum target properties for the library.
# See https://cmake.org/cmake/help/latest/command/set_target_properties.html
# Properties are: 1) Library version number 2) list of corresponding public headers
#set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES VERSION       ${${PROJECT_NAME}_VERSION}
#                                                        PUBLIC_HEADER ${${LIBRARY_TARGET_NAME}_HDRS}
#)


# Option 1:
# target_include_directories(${LIBRARY_TARGET_NAME} PUBLIC ${${LIBRARY_TARGET_NAME}_CERT_INCLUDES})
# Option 2:
target_include_directories(${LIBRARY_TARGET_NAME} PUBLIC
        # only when building from the source tree
        $<BUILD_INTERFACE:${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}>
        # ---------------------$<INSTALL_INTERFACE:${INSTALL_CMAKE_DIR}"/${LIBRARY_TARGET_NAME}EstimationRosen/include">) or:
        # $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}>)
        # only when using the lib from the install path
        $<INSTALL_INTERFACE:include>
        ${${LIBRARY_TARGET_NAME}_ADD_INCLUDES}
)


target_link_libraries(${LIBRARY_TARGET_NAME}
        PUBLIC  ${BLAS_LIBRARIES} ${M} ${LAPACK} 
)




if(${CODE_PROFILING})
  set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES COMPILE_FLAGS "-pg -g" LINK_FLAGS "-pg -g")
endif()


set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS ${LIBRARY_TARGET_NAME}
  EXPORT  ${LIBRARY_TARGET_NAME_EXPORT}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX}
  INCLUDES DESTINATION "include"
  PUBLIC_HEADER DESTINATION "include/${LIBRARY_TARGET_NAME}"
 )



# Build the example executable
IF(BUILD_${LIBRARY_TARGET_NAME}_EXAMPLE)
message(STATUS "Adding examples to build")
add_subdirectory(examples)
endif()

include(cmake/setup_installation.cmake)

