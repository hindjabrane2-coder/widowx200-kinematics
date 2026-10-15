import numpy as np


# Parametres DH du WidowX-200 (4 DDL)
# Convention DH modifiee : a, alpha, d, theta_offset

NB_AXES = 4

# Longueurs des segments (mm)
L1 = 110.0    # base vers epaule
L2 = 200.0    # bras superieur
L3 = 200.0    # avant-bras
L4 = 130.0    # poignet vers effecteur

DH_TABLE = np.array([
    [0.0,    np.pi / 2, L1,  0.0],       # axe 1 : rotation base
    [L2,     0.0,       0.0, -np.pi / 2], # axe 2 : epaule
    [L3,     0.0,       0.0, 0.0],        # axe 3 : coude
    [L4,     0.0,       0.0, 0.0],        # axe 4 : poignet
])

# Limites articulaires (rad)
JOINT_LIMITS_MIN = np.array([-np.pi, -1.85, -2.15, -np.pi])
JOINT_LIMITS_MAX = np.array([np.pi, 1.60, 1.53, np.pi])

# Vitesses maximales par axe (rad/s) - Dynamixel XL430
JOINT_VEL_MAX = np.array([3.5, 3.5, 3.5, 3.5])

# Couples maximaux (N.m)
JOINT_TORQUE_MAX = np.array([1.5, 1.5, 1.5, 1.5])

# Position de repos
HOME_POSITION = np.array([0.0, 0.0, 0.0, 0.0])

# Tolerances
POSITION_TOLERANCE = 1e-3    # mm
ORIENTATION_TOLERANCE = 1e-4 # rad
SINGULARITY_THRESHOLD = 1e-4


def dans_limites(q):
    """Verifie si une configuration articulaire respecte les butees."""
    return np.all(q >= JOINT_LIMITS_MIN) and np.all(q <= JOINT_LIMITS_MAX)


def saturer_articulations(q):
    """Sature les positions articulaires aux butees mecaniques."""
    return np.clip(q, JOINT_LIMITS_MIN, JOINT_LIMITS_MAX)
