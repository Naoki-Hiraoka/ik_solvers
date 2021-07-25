find_package(catkin REQUIRED COMPONENTS
  vclip #VclipLinkPair.cpp, VclipCollisionConstraint.cpp
  cddeigen #CddSCFRConstraint.cpp #CddSCFRConstraint2.cpp
  qhulleigen #CddSCFRConstraint.cpp #CddSCFRConstraint2.cpp
  clpeigen #LPSCFRConstraint.cpp
  )

#main
set(target IKConstraints)
set(sources
  IKConstraint.cpp
  PositionConstraint.cpp
  MinMaxJointConstraint.cpp
  CollisionConstraint.cpp
  AISTCollisionConstraint.cpp
  VclipLinkPair.cpp
  VClipCollisionConstraint.cpp
  SCFRConstraint.cpp
  CddSCFRConstraint.cpp
  CddSCFRConstraint2.cpp
  LPSCFRConstraint.cpp
  )

add_library(${target} SHARED ${sources})
target_include_directories(${target} PUBLIC ${catkin_INCLUDE_DIRS})

target_link_libraries(${target} PUBLIC
  ${catkin_LIBRARIES}
  RobotConfig
  )

target_compile_definitions(${target} PUBLIC ${CHOREONOID_DEFINITIONS})
