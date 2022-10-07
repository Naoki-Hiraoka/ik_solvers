add_compile_options(-std=c++14)
set(CMAKE_BUILD_TYPE RelWithDebInfo)

find_package(catkin REQUIRED COMPONENTS
  ik_constraint
  prioritized_qp_osqp
  roslib
  )
find_package(choreonoid REQUIRED)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${CHOREONOID_INCLUDE_DIRS}
)

# 相対パスを絶対パスに直す
set(CHOREONOID_BODY_LIBRARIES_ABS)
foreach(lib ${CHOREONOID_BODY_LIBRARIES})
  find_library(${lib}_abs NAMES ${lib} PATHS ${CHOREONOID_LIBRARY_DIRS})
  set(CHOREONOID_BODY_LIBRARIES_ABS ${CHOREONOID_BODY_LIBRARIES_ABS} ${${lib}_abs})
endforeach(lib)

add_executable(sampleSR1
  SampleSR1.cpp
)

target_link_libraries(sampleSR1
  ${catkin_LIBRARIES}
  ${CHOREONOID_BODY_LIBRARIES_ABS}
  ${PROJECT_NAME}
  )

install(TARGETS sampleSR1
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  )
