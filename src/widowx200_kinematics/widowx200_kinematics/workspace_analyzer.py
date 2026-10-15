import numpy as np
import rclpy
from rclpy.node import Node
from std_msgs.msg import Float64MultiArray

from . import robot_config as cfg
from .espace_travail import (
    analyse_monte_carlo, volume_espace_travail,
    cartographie_dexterite, portee_maximale, portee_minimale,
    hauteur_maximale, zones_singuliere,
)


class WorkspaceAnalyzer(Node):

    def __init__(self):
        super().__init__('workspace_analyzer')

        self.declare_parameter('monte_carlo_samples', 50000)
        self.declare_parameter('dexterity_resolution', 100)
        self.declare_parameter('dexterity_z_height', 200.0)
        self.declare_parameter('output_directory', '/tmp/widowx200_workspace')

        self.pub_workspace = self.create_publisher(
            Float64MultiArray, 'workspace_points', 10
        )

        self.get_logger().info('Analyseur d\'espace de travail demarre')
        self.get_logger().info(
            f'Portee : {portee_minimale():.0f} - {portee_maximale():.0f} mm, '
            f'hauteur max : {hauteur_maximale():.0f} mm'
        )

    def analyser(self):
        """Lance l'analyse complete de l'espace de travail."""
        n_samples = self.get_parameter('monte_carlo_samples').value
        self.get_logger().info(
            f'Analyse Monte Carlo avec {n_samples} echantillons...'
        )

        resultats = analyse_monte_carlo(n_samples)
        positions = resultats['positions']
        manip = resultats['manipulabilite']

        vol = volume_espace_travail(positions)
        self.get_logger().info(f'Volume de l\'espace de travail : {vol:.2f} litres')
        self.get_logger().info(f'Points valides : {resultats["n_valides"]}')

        # Statistiques de manipulabilite
        manip_valide = manip[manip > 0]
        if len(manip_valide) > 0:
            self.get_logger().info(
                f'Manipulabilite : moy={np.mean(manip_valide):.4f}, '
                f'max={np.max(manip_valide):.4f}, '
                f'min={np.min(manip_valide):.6f}'
            )

        # Cartographie de dexterite
        z = self.get_parameter('dexterity_z_height').value
        res = self.get_parameter('dexterity_resolution').value
        self.get_logger().info(
            f'Cartographie de dexterite a z={z:.0f} mm...'
        )
        x, y, carte = cartographie_dexterite(z, res)

        n_atteignables = np.sum(~np.isnan(carte))
        self.get_logger().info(
            f'Points atteignables dans le plan z={z:.0f} : '
            f'{n_atteignables}/{res*res}'
        )

        # Zones singulieres
        self.get_logger().info('Recherche des zones singulieres...')
        singularites = zones_singuliere(n_samples // 5)
        self.get_logger().info(
            f'{len(singularites)} configurations singulieres detectees'
        )

        # Publication des points pour visualisation
        self._publier_points(positions)

        # Export CSV
        output_dir = self.get_parameter('output_directory').value
        self._exporter_resultats(output_dir, positions, manip, x, y, carte)

        return resultats

    def _publier_points(self, positions):
        """Publie les points de l'espace de travail."""
        msg = Float64MultiArray()
        msg.data = positions.flatten().tolist()
        self.pub_workspace.publish(msg)

    def _exporter_resultats(self, output_dir, positions, manip, x, y, carte):
        """Exporte les resultats dans des fichiers CSV."""
        import os
        os.makedirs(output_dir, exist_ok=True)

        # Points de l'espace de travail
        fichier_points = os.path.join(output_dir, 'workspace_points.csv')
        with open(fichier_points, 'w') as f:
            f.write('x_mm,y_mm,z_mm,manipulabilite\n')
            for i in range(len(positions)):
                f.write(
                    f'{positions[i, 0]:.2f},{positions[i, 1]:.2f},'
                    f'{positions[i, 2]:.2f},{manip[i]:.6f}\n'
                )

        # Carte de dexterite
        fichier_carte = os.path.join(output_dir, 'dexterity_map.csv')
        with open(fichier_carte, 'w') as f:
            f.write('x_mm,y_mm,manipulabilite\n')
            for i, xi in enumerate(x):
                for j, yj in enumerate(y):
                    if not np.isnan(carte[j, i]):
                        f.write(f'{xi:.2f},{yj:.2f},{carte[j, i]:.6f}\n')

        self.get_logger().info(f'Resultats exportes dans {output_dir}')


def main(args=None):
    rclpy.init(args=args)
    node = WorkspaceAnalyzer()
    node.analyser()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
