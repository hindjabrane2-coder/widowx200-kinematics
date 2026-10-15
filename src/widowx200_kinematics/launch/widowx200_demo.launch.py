import os
from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch_ros.actions import Node


def generate_launch_description():
    config = os.path.join(
        get_package_share_directory('widowx200_kinematics'),
        'config',
        'robot_params.yaml',
    )

    return LaunchDescription([
        Node(
            package='widowx200_kinematics',
            executable='kinematics_server',
            name='kinematics_server',
            parameters=[config],
            output='screen',
        ),
        Node(
            package='widowx200_kinematics',
            executable='trajectory_planner',
            name='trajectory_planner',
            parameters=[config],
            output='screen',
        ),
    ])
