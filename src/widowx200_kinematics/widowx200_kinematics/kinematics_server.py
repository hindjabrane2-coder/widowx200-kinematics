import numpy as np
import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Pose, Point, Quaternion

from widowx200_interfaces.srv import ComputeFK, ComputeIK

from . import robot_config as cfg
from .cinematique import (
    mgd, manipulabilite, mgi_geometrique, mgi_numerique,
    choisir_meilleure_solution,
)


def rotation_vers_quaternion(R):
    """Convertit une matrice de rotation 3x3 en quaternion (x, y, z, w)."""
    trace = R[0, 0] + R[1, 1] + R[2, 2]
    if trace > 0:
        s = 0.5 / np.sqrt(trace + 1.0)
        w = 0.25 / s
        x = (R[2, 1] - R[1, 2]) * s
        y = (R[0, 2] - R[2, 0]) * s
        z = (R[1, 0] - R[0, 1]) * s
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = 2.0 * np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
        w = (R[2, 1] - R[1, 2]) / s
        x = 0.25 * s
        y = (R[0, 1] + R[1, 0]) / s
        z = (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = 2.0 * np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
        w = (R[0, 2] - R[2, 0]) / s
        x = (R[0, 1] + R[1, 0]) / s
        y = 0.25 * s
        z = (R[1, 2] + R[2, 1]) / s
    else:
        s = 2.0 * np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
        w = (R[1, 0] - R[0, 1]) / s
        x = (R[0, 2] + R[2, 0]) / s
        y = (R[1, 2] + R[2, 1]) / s
        z = 0.25 * s
    return x, y, z, w


def quaternion_vers_rotation(qx, qy, qz, qw):
    """Convertit un quaternion en matrice de rotation 3x3."""
    R = np.zeros((3, 3))
    R[0, 0] = 1 - 2 * (qy ** 2 + qz ** 2)
    R[0, 1] = 2 * (qx * qy - qz * qw)
    R[0, 2] = 2 * (qx * qz + qy * qw)
    R[1, 0] = 2 * (qx * qy + qz * qw)
    R[1, 1] = 1 - 2 * (qx ** 2 + qz ** 2)
    R[1, 2] = 2 * (qy * qz - qx * qw)
    R[2, 0] = 2 * (qx * qz - qy * qw)
    R[2, 1] = 2 * (qy * qz + qx * qw)
    R[2, 2] = 1 - 2 * (qx ** 2 + qy ** 2)
    return R


class KinematicsServer(Node):

    def __init__(self):
        super().__init__('kinematics_server')

        self.declare_parameter('ik_max_iterations', 200)
        self.declare_parameter('ik_tolerance', 1e-3)
        self.declare_parameter('ik_lambda_init', 0.01)

        self.srv_fk = self.create_service(
            ComputeFK, 'compute_fk', self.callback_fk
        )
        self.srv_ik = self.create_service(
            ComputeIK, 'compute_ik', self.callback_ik
        )

        self.get_logger().info(
            'Serveur de cinematique demarre (FK et IK sur WidowX-200)'
        )

    def callback_fk(self, request, response):
        q = np.array(request.joint_positions)

        if len(q) != cfg.NB_AXES:
            self.get_logger().warn(
                f'FK : attendu {cfg.NB_AXES} articulations, recu {len(q)}'
            )
            response.success = False
            return response

        if not cfg.dans_limites(q):
            self.get_logger().warn('FK : configuration hors limites articulaires')

        T = mgd(q)
        pos = T[:3, 3]
        R = T[:3, :3]
        qx, qy, qz, qw = rotation_vers_quaternion(R)

        response.pose = Pose()
        response.pose.position = Point(
            x=float(pos[0]), y=float(pos[1]), z=float(pos[2])
        )
        response.pose.orientation = Quaternion(
            x=float(qx), y=float(qy), z=float(qz), w=float(qw)
        )
        response.manipulability = float(manipulabilite(q))
        response.success = True

        self.get_logger().info(
            f'FK : q={np.round(np.degrees(q), 1)} -> '
            f'pos=({pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f}) mm, '
            f'manip={response.manipulability:.4f}'
        )

        return response

    def callback_ik(self, request, response):
        cible = np.array([
            request.target_pose.position.x,
            request.target_pose.position.y,
            request.target_pose.position.z,
        ])

        tol = request.tolerance if request.tolerance > 0 else (
            self.get_parameter('ik_tolerance').value
        )

        # Solution geometrique en priorite
        solutions = mgi_geometrique(cible[0], cible[1], cible[2])

        q_init = None
        if len(request.seed_joints) == cfg.NB_AXES:
            q_init = np.array(request.seed_joints)

        q_sol = choisir_meilleure_solution(solutions, q_init)

        if q_sol is not None:
            T_verif = mgd(q_sol)
            erreur = np.linalg.norm(cible - T_verif[:3, 3])
            response.joint_positions = q_sol.tolist()
            response.position_error = float(erreur)
            response.iterations = 0
            response.success = True
        else:
            # Repli sur la methode numerique
            max_iter = self.get_parameter('ik_max_iterations').value
            lam = self.get_parameter('ik_lambda_init').value

            q_num, erreur, n_iter, ok = mgi_numerique(
                cible, q_init=q_init, max_iter=max_iter,
                tol=tol, lambda_init=lam,
            )
            response.joint_positions = q_num.tolist()
            response.position_error = float(erreur)
            response.iterations = int(n_iter)
            response.success = bool(ok)

        if response.success:
            self.get_logger().info(
                f'IK : cible=({cible[0]:.1f}, {cible[1]:.1f}, {cible[2]:.1f}) '
                f'-> q={np.round(np.degrees(response.joint_positions), 1)}, '
                f'erreur={response.position_error:.3f} mm'
            )
        else:
            self.get_logger().warn(
                f'IK : echec pour cible ({cible[0]:.1f}, {cible[1]:.1f}, '
                f'{cible[2]:.1f}), erreur={response.position_error:.3f} mm'
            )

        return response


def main(args=None):
    rclpy.init(args=args)
    node = KinematicsServer()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
