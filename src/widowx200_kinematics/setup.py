from setuptools import find_packages, setup
import os
from glob import glob

package_name = 'widowx200_kinematics'

setup(
    name=package_name,
    version='1.0.0',
    packages=find_packages(exclude=['test']),
    data_files=[
        ('share/ament_index/resource_index/packages', ['resource/' + package_name]),
        ('share/' + package_name, ['package.xml']),
        (os.path.join('share', package_name, 'launch'), glob('launch/*.py')),
        (os.path.join('share', package_name, 'config'), glob('config/*.yaml')),
    ],
    install_requires=['setuptools', 'numpy'],
    zip_safe=True,
    maintainer='Hind Jabrane',
    maintainer_email='hindhindo067@gmail.com',
    description='Cinematique et trajectoire du WidowX-200',
    license='MIT',
    entry_points={
        'console_scripts': [
            'kinematics_server = widowx200_kinematics.kinematics_server:main',
            'trajectory_planner = widowx200_kinematics.trajectory_planner:main',
            'workspace_analyzer = widowx200_kinematics.workspace_analyzer:main',
        ],
    },
)
