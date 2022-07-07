project(RotPriorRPp-Examples CXX)


# Find Eigen library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")
find_package(Eigen3 3.3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()

# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})


add_library(exp_gen ../utils/generatePointCloud.h ../utils/generatePointCloud.cpp)


# Fast Essential Matrix Estimation
add_executable(example_basic ${CMAKE_CURRENT_SOURCE_DIR}/example_basic.cpp)
target_link_libraries(example_basic RotPriorRPp exp_gen)

# Fast Iterative certification
add_executable(example_certifier ${CMAKE_CURRENT_SOURCE_DIR}/example_certifier.cpp)
target_link_libraries(example_certifier RotPriorRPp IterCertAlg exp_gen)


# staircase-alike certification
add_executable(example_rpp_staircase ${CMAKE_CURRENT_SOURCE_DIR}/example_rpp_staircase.cpp)
target_link_libraries(example_rpp_staircase RotPriorRPp IterCertAlg exp_gen)

