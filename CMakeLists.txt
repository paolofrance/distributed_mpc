cmake_minimum_required(VERSION 3.1)

project(distributed_mpc)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(catkin REQUIRED COMPONENTS
roscpp
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)


catkin_package(
 INCLUDE_DIRS include
#  LIBRARIES gt_traj_deformation
 CATKIN_DEPENDS 
 DEPENDS EIGEN3
)

include_directories( include ${catkin_INCLUDE_DIRS} )



add_library(
d_mpc src/distributed_mpc/d_mpc.cpp
)

add_executable(mpc_node src/distributed_mpc/mpc_node.cpp)
target_link_libraries(mpc_node d_mpc ${catkin_LIBRARIES})


