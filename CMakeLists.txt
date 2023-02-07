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
 LIBRARIES ${PROJECT_NAME}
 CATKIN_DEPENDS 
 DEPENDS EIGEN3
)

include_directories( include ${catkin_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})



add_library(
${PROJECT_NAME} src/distributed_mpc/d_mpc.cpp
)
target_link_libraries(${PROJECT_NAME}  ${catkin_LIBRARIES} Eigen3::Eigen)

install(
  TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)
install(
  DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

add_executable(mpc_node src/distributed_mpc/mpc_node.cpp)
target_link_libraries(mpc_node ${PROJECT_NAME} ${catkin_LIBRARIES})


