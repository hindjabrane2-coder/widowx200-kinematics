cmake_minimum_required(VERSION 3.8)
project(widowx200_interfaces)

find_package(ament_cmake REQUIRED)
find_package(rosidl_default_generators REQUIRED)
find_package(geometry_msgs REQUIRED)

rosidl_generate_interfaces(${PROJECT_NAME}
  "msg/TrajectoryPoint.msg"
  "srv/ComputeFK.srv"
  "srv/ComputeIK.srv"
  "srv/PlanTrajectory.srv"
  DEPENDENCIES geometry_msgs
)

ament_package()
