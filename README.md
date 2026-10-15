# WidowX-200 Kinematics (ROS2)

Package ROS2 pour la cinematique, la planification de trajectoire et l'analyse de l'espace de travail du robot WidowX-200 (Interbotix). Le robot utilise 4 servomoteurs Dynamixel XL430 et possede 4 degres de liberte.


## Architecture

Le workspace contient deux packages :

```
src/
    widowx200_interfaces/       Messages et services ROS2 personnalises
        msg/TrajectoryPoint.msg
        srv/ComputeFK.srv
        srv/ComputeIK.srv
        srv/PlanTrajectory.srv

    widowx200_kinematics/       Noeuds ROS2 et bibliotheque de calcul
        widowx200_kinematics/
            robot_config.py         Parametres DH, limites articulaires
            cinematique.py          MGD, MGI geometrique et numerique, jacobienne
            profil_vitesse.py       Profils trapeze et S-curve, spline cubique
            trajectoire.py          MoveJ, MoveL, MoveC
            espace_travail.py       Monte Carlo, dexterite, singularites
            kinematics_server.py    Noeud serveur FK/IK
            trajectory_planner.py   Noeud planificateur de trajectoire
            workspace_analyzer.py   Noeud analyse de l'espace de travail
        launch/
            widowx200_demo.launch.py
        config/
            robot_params.yaml
```


## Parametres DH

| Axe | a (mm) | alpha (rad) | d (mm) | offset (rad) |
|-----|--------|-------------|--------|--------------|
| 1   | 0      | pi/2        | 110    | 0            |
| 2   | 200    | 0           | 0      | -pi/2        |
| 3   | 200    | 0           | 0      | 0            |
| 4   | 130    | 0           | 0      | 0            |

Longueurs : L1 = 110, L2 = 200, L3 = 200, L4 = 130 mm


## Noeuds ROS2

### kinematics_server

Expose deux services :

`compute_fk` (ComputeFK) : calcule la position et l'orientation de l'effecteur a partir des positions articulaires. Retourne aussi l'indice de manipulabilite de Yoshikawa.

`compute_ik` (ComputeIK) : resout la cinematique inverse. Essaie d'abord la solution geometrique analytique (coude haut et coude bas), puis bascule sur la methode numerique de Newton-Raphson avec amortissement de Levenberg-Marquardt si la solution geometrique n'existe pas.

### trajectory_planner

Service `plan_trajectory` (PlanTrajectory) avec deux modes :

Mode "joint" : interpolation dans l'espace articulaire avec profil de vitesse trapezoldal ou en S (jerk limite). L'axe le plus lent cadence le mouvement.

Mode "cartesian" : interpolation lineaire dans l'espace cartesien. L'effecteur suit une ligne droite. La cinematique inverse est resolue a chaque point intermediaire.

Publie l'etat articulaire sur le topic `joint_states` (sensor_msgs/JointState).

### workspace_analyzer

Analyse Monte Carlo de l'espace de travail : volume, cartographie de dexterite dans un plan horizontal, detection des zones singulieres. Exporte les resultats en CSV.


## Cinematique inverse

Deux methodes implementees :

Solution geometrique : resolution analytique pour 4 DDL. Le poignet est calcule pour maintenir l'effecteur horizontal. Deux solutions possibles (coude haut, coude bas), la plus proche de la configuration courante est retenue.

Solution numerique : Newton-Raphson avec amortissement de Levenberg-Marquardt. La jacobienne est calculee par differences finies. Le facteur d'amortissement lambda est ajuste dynamiquement : reduit quand l'erreur diminue, augmente quand elle augmente. Saturation automatique aux butees articulaires a chaque iteration.


## Profils de vitesse

Trapezoldal : acceleration constante, phase de croisiere a vitesse max, deceleration constante. Passe en profil triangulaire si la distance est trop courte pour atteindre la vitesse maximale.

S-curve : 7 phases avec jerk borne. L'acceleration monte et descend progressivement, ce qui reduit les vibrations mecaniques et les a-coups sur les servomoteurs.

Spline cubique naturelle : interpolation douce par resolution du systeme tridiagonal. Parametre cumulatif base sur la distance entre les points.


## Trajectoires

MoveJ : mouvement articulaire point a point. Chaque axe parcourt sa distance angulaire proportionnellement au meme profil temporel.

MoveL : mouvement lineaire cartesien. Discretisation de la ligne droite, resolution de la MGI a chaque point. Si la solution geometrique echoue, bascule sur la methode numerique.

MoveC : mouvement circulaire defini par trois points (depart, intermediaire, arrivee). Calcul du centre et du rayon du cercle par intersection des mediatrices, puis interpolation le long de l'arc.


## Compilation et lancement

```bash
cd widowx200-kinematics
colcon build
source install/setup.bash
ros2 launch widowx200_kinematics widowx200_demo.launch.py
```

Tester le service FK :
```bash
ros2 service call /compute_fk widowx200_interfaces/srv/ComputeFK \
  "{joint_positions: [0.0, 0.5, -0.3, 0.0]}"
```

Tester le service IK :
```bash
ros2 service call /compute_ik widowx200_interfaces/srv/ComputeIK \
  "{target_pose: {position: {x: 300.0, y: 0.0, z: 200.0}}, tolerance: 0.01}"
```

Planifier une trajectoire articulaire :
```bash
ros2 service call /plan_trajectory widowx200_interfaces/srv/PlanTrajectory \
  "{mode: 'joint', start_joints: [0,0,0,0], goal_joints: [1.0,0.5,-0.3,0.2], \
    max_velocity: 1.0, max_acceleration: 2.0, profile_type: 'scurve', num_points: 100}"
```

Lancer l'analyse de l'espace de travail :
```bash
ros2 run widowx200_kinematics workspace_analyzer
```


## Dependances

ROS2 Humble ou superieur, Python 3, numpy.
