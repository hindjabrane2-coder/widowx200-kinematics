import numpy as np
import rclpy
from rclpy.node import Node
from sensor_msgs.msg import JointState

from widowx200_interfaces.srv import PlanTrajectory
from widowx200_interfaces.msg import TrajectoryPoint

from . import robot_config as cfg
from .cinematique import mgd, mgi_geometrique, choisir_meilleure_solution
from .trajectoire import move_j, move_l, move_c


class TrajectoryPlanner(Node):

    def __init__(self):
        super().__init__('trajectory_planner')

        self.declare_parameter('default_joint_velocity', 1.0)
        self.declare_parameter('default_joint_acceleration', 2.0)
        self.declare_parameter('default_cartesian_velocity', 50.0)
        self.declare_parameter('default_cartesian_acceleration', 200.0)
        self.declare_parameter('default_num_points', 200)

        self.srv_plan = self.create_service(
            PlanTrajectory, 'plan_trajectory', self.callback_plan
        )

        # Publication de l'etat articulaire courant
        self.pub_joint_state = self.create_publisher(
            JointState, 'joint_states', 10
        )

        self.q_current = cfg.HOME_POSITION.copy()

        self.get_logger().info('Planificateur de trajectoire demarre')

    def callback_plan(self, request, response):
        mode = request.mode.lower()
        n_pts = request.num_points if request.num_points > 0 else (
            self.get_parameter('default_num_points').value
        )

        q_start = np.array(request.start_joints) if len(request.start_joints) == cfg.NB_AXES else self.q_current

        if mode == 'joint':
            response = self._plan_joint(request, response, q_start, n_pts)
        elif mode == 'cartesian' or mode == 'linear':
            response = self._plan_linear(request, response, q_start, n_pts)
        else:
            response.success = False
            response.message = f'Mode inconnu : {mode}. Utiliser "joint" ou "cartesian".'

        return response

    def _plan_joint(self, request, response, q_start, n_pts):
        if len(request.goal_joints) != cfg.NB_AXES:
            response.success = False
            response.message = f'goal_joints doit contenir {cfg.NB_AXES} valeurs'
            return response

        q_goal = np.array(request.goal_joints)

        v_max = request.max_velocity if request.max_velocity > 0 else (
            self.get_parameter('default_joint_velocity').value
        )
        a_max = request.max_acceleration if request.max_acceleration > 0 else (
            self.get_parameter('default_joint_acceleration').value
        )

        profil = 'scurve' if request.profile_type.lower() in ('scurve', 's-curve', 's_curve') else 'trapeze'

        trajectoire, temps = move_j(q_start, q_goal, v_max, a_max, profil, n_pts)

        response.points = self._convertir_points(trajectoire, temps)
        response.duration = float(temps[-1]) if len(temps) > 0 else 0.0
        response.success = True
        response.message = f'Trajectoire articulaire ({profil}), {len(temps)} points'

        self.q_current = q_goal.copy()
        self._publier_joint_state(q_goal)

        self.get_logger().info(
            f'MoveJ : {np.round(np.degrees(q_start), 1)} -> '
            f'{np.round(np.degrees(q_goal), 1)}, '
            f'duree={response.duration:.2f}s'
        )

        return response

    def _plan_linear(self, request, response, q_start, n_pts):
        pos_goal = np.array([
            request.goal_pose.position.x,
            request.goal_pose.position.y,
            request.goal_pose.position.z,
        ])

        v_cart = request.max_velocity if request.max_velocity > 0 else (
            self.get_parameter('default_cartesian_velocity').value
        )
        a_cart = request.max_acceleration if request.max_acceleration > 0 else (
            self.get_parameter('default_cartesian_acceleration').value
        )

        trajectoire, temps, succes = move_l(q_start, pos_goal, v_cart, a_cart, n_pts)

        response.points = self._convertir_points(trajectoire, temps)
        response.duration = float(temps[-1]) if len(temps) > 0 else 0.0
        response.success = succes

        if succes:
            response.message = f'Trajectoire lineaire, {len(temps)} points'
            self.q_current = trajectoire[-1].copy()
            self._publier_joint_state(self.q_current)
        else:
            response.message = 'Trajectoire lineaire partiellement resolue (singularites rencontrees)'

        self.get_logger().info(
            f'MoveL : vers ({pos_goal[0]:.1f}, {pos_goal[1]:.1f}, {pos_goal[2]:.1f}), '
            f'duree={response.duration:.2f}s, succes={succes}'
        )

        return response

    def _convertir_points(self, trajectoire, temps):
        """Convertit les tableaux numpy en messages TrajectoryPoint."""
        points = []
        for i in range(len(temps)):
            pt = TrajectoryPoint()
            pt.time_from_start = float(temps[i])
            pt.positions = trajectoire[i].tolist()

            # Vitesses par differences finies
            if i == 0:
                vel = np.zeros(cfg.NB_AXES)
            elif i == len(temps) - 1:
                vel = np.zeros(cfg.NB_AXES)
            else:
                dt = temps[i + 1] - temps[i - 1]
                if dt > 0:
                    vel = (trajectoire[i + 1] - trajectoire[i - 1]) / dt
                else:
                    vel = np.zeros(cfg.NB_AXES)
            pt.velocities = vel.tolist()
            pt.accelerations = []
            points.append(pt)

        return points

    def _publier_joint_state(self, q):
        msg = JointState()
        msg.header.stamp = self.get_clock().now().to_msg()
        msg.name = [f'joint_{i+1}' for i in range(cfg.NB_AXES)]
        msg.position = q.tolist()
        self.pub_joint_state.publish(msg)


def main(args=None):
    rclpy.init(args=args)
    node = TrajectoryPlanner()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
