import numpy as np
from . import robot_config as cfg
from .cinematique import mgd, manipulabilite


def analyse_monte_carlo(n_echantillons=50000):
    """Analyse de l'espace de travail par methode de Monte Carlo.

    Tire des configurations articulaires aleatoires, calcule la position
    de l'effecteur, et retourne les points atteignables avec leur
    indice de manipulabilite.
    """
    positions = np.zeros((n_echantillons, 3))
    manip = np.zeros(n_echantillons)
    valides = np.zeros(n_echantillons, dtype=bool)

    for i in range(n_echantillons):
        q = cfg.JOINT_LIMITS_MIN + np.random.random(cfg.NB_AXES) * (
            cfg.JOINT_LIMITS_MAX - cfg.JOINT_LIMITS_MIN
        )
        T = mgd(q)
        positions[i] = T[:3, 3]
        manip[i] = manipulabilite(q)
        valides[i] = True

    return {
        'positions': positions[valides],
        'manipulabilite': manip[valides],
        'n_valides': int(np.sum(valides)),
    }


def volume_espace_travail(positions):
    """Estime le volume de l'espace de travail a partir des points Monte Carlo.

    Utilise une grille 3D pour compter les cellules occupees.
    """
    resolution = 10.0  # mm

    mins = np.min(positions, axis=0)
    maxs = np.max(positions, axis=0)

    nx = int((maxs[0] - mins[0]) / resolution) + 1
    ny = int((maxs[1] - mins[1]) / resolution) + 1
    nz = int((maxs[2] - mins[2]) / resolution) + 1

    grille = np.zeros((nx, ny, nz), dtype=bool)

    for p in positions:
        ix = int((p[0] - mins[0]) / resolution)
        iy = int((p[1] - mins[1]) / resolution)
        iz = int((p[2] - mins[2]) / resolution)
        ix = min(ix, nx - 1)
        iy = min(iy, ny - 1)
        iz = min(iz, nz - 1)
        grille[ix, iy, iz] = True

    n_cellules = np.sum(grille)
    volume = n_cellules * resolution ** 3
    return volume / 1e9  # litres


def cartographie_dexterite(z_coupe, n_points=100):
    """Cartographie de dexterite dans un plan horizontal a la hauteur z_coupe.

    Pour chaque point (x, y) d'une grille, on calcule l'indice de
    manipulabilite si le point est atteignable.
    """
    r_max = cfg.L2 + cfg.L3 + cfg.L4
    x = np.linspace(-r_max, r_max, n_points)
    y = np.linspace(-r_max, r_max, n_points)

    carte = np.full((n_points, n_points), np.nan)

    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            r = np.sqrt(xi ** 2 + yj ** 2)
            if r < 10 or r > r_max:
                continue

            q1 = np.arctan2(yj, xi)
            rw = r - cfg.L4
            zw = z_coupe - cfg.L1
            d = np.sqrt(rw ** 2 + zw ** 2)

            if d > cfg.L2 + cfg.L3 or d < abs(cfg.L2 - cfg.L3):
                continue

            cos_q3 = (d ** 2 - cfg.L2 ** 2 - cfg.L3 ** 2) / (2 * cfg.L2 * cfg.L3)
            if abs(cos_q3) > 1:
                continue

            q3 = np.arccos(np.clip(cos_q3, -1, 1))
            beta = np.arctan2(zw, rw)
            phi = np.arctan2(cfg.L3 * np.sin(q3), cfg.L2 + cfg.L3 * np.cos(q3))
            q2 = beta - phi
            q4 = -(q2 + q3)

            q = np.array([q1, q2, q3, q4])
            if cfg.dans_limites(q):
                carte[j, i] = manipulabilite(q)

    return x, y, carte


def portee_maximale():
    """Calcule la portee maximale du robot (bras completement tendu)."""
    return cfg.L2 + cfg.L3 + cfg.L4


def portee_minimale():
    """Calcule la portee minimale (bras replie au maximum)."""
    return abs(cfg.L2 - cfg.L3) + cfg.L4


def hauteur_maximale():
    """Hauteur maximale atteignable."""
    return cfg.L1 + cfg.L2 + cfg.L3 + cfg.L4


def zones_singuliere(n_echantillons=10000):
    """Identifie les configurations proches des singularites."""
    configs_singulieres = []

    for _ in range(n_echantillons):
        q = cfg.JOINT_LIMITS_MIN + np.random.random(cfg.NB_AXES) * (
            cfg.JOINT_LIMITS_MAX - cfg.JOINT_LIMITS_MIN
        )
        m = manipulabilite(q)
        if m < cfg.SINGULARITY_THRESHOLD * 10:
            T = mgd(q)
            configs_singulieres.append({
                'q': q.copy(),
                'position': T[:3, 3].copy(),
                'manipulabilite': m,
            })

    return configs_singulieres
