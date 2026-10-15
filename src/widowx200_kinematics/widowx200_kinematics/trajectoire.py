import numpy as np
from . import robot_config as cfg
from .cinematique import mgd, mgi_geometrique, mgi_numerique, choisir_meilleure_solution
from .profil_vitesse import profil_trapeze, profil_s_curve


def move_j(q_start, q_goal, v_max=1.0, a_max=2.0, profil='trapeze', n_points=200):
    """Mouvement articulaire (interpolation dans l'espace des joints).

    Chaque axe suit le meme profil temporel, cadence par l'axe le plus lent.
    """
    dq = q_goal - q_start
    dq_max = np.max(np.abs(dq))

    if dq_max < 1e-10:
        return np.tile(q_start, (1, 1)), np.zeros(1)

    if profil == 'scurve':
        t, s, _, _ = profil_s_curve(1.0, v_max / dq_max, a_max / dq_max,
                                     a_max * 10 / dq_max, n_points)
    else:
        t, s, _, _ = profil_trapeze(1.0, v_max / dq_max, a_max / dq_max, n_points)

    trajectoire = np.zeros((n_points, cfg.NB_AXES))
    for i in range(n_points):
        trajectoire[i] = q_start + s[i] * dq

    return trajectoire, t


def move_l(q_start, pos_goal, v_cart=50.0, a_cart=200.0, n_points=200):
    """Mouvement lineaire dans l'espace cartesien.

    L'effecteur suit une ligne droite de sa position actuelle vers pos_goal.
    La cinematique inverse est resolue a chaque point intermediaire.
    """
    T_start = mgd(q_start)
    pos_start = T_start[:3, 3]

    direction = pos_goal - pos_start
    distance = np.linalg.norm(direction)

    if distance < cfg.POSITION_TOLERANCE:
        return np.tile(q_start, (1, 1)), np.zeros(1), True

    t, s, _, _ = profil_trapeze(distance, v_cart, a_cart, n_points)

    trajectoire = np.zeros((n_points, cfg.NB_AXES))
    trajectoire[0] = q_start
    q_prev = q_start.copy()
    succes = True

    for i in range(1, n_points):
        ratio = s[i] / distance if distance > 0 else 0
        pos_inter = pos_start + ratio * direction

        solutions = mgi_geometrique(pos_inter[0], pos_inter[1], pos_inter[2])
        q_sol = choisir_meilleure_solution(solutions, q_prev)

        if q_sol is None:
            q_sol, err, _, ok = mgi_numerique(pos_inter, q_init=q_prev)
            if not ok:
                succes = False
                q_sol = q_prev

        trajectoire[i] = q_sol
        q_prev = q_sol

    return trajectoire, t, succes


def move_c(q_start, pos_via, pos_goal, v_cart=50.0, a_cart=200.0, n_points=200):
    """Mouvement circulaire passant par trois points.

    Le cercle est defini par la position actuelle, un point intermediaire
    et le point final. L'interpolation suit l'arc de cercle.
    """
    T_start = mgd(q_start)
    pos_start = T_start[:3, 3]

    p1 = pos_start
    p2 = pos_via
    p3 = pos_goal

    # Calcul du centre et du rayon du cercle
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    norm_n = np.linalg.norm(normal)

    if norm_n < 1e-10:
        # Points colineaires, on fait un MoveL
        return move_l(q_start, pos_goal, v_cart, a_cart, n_points)

    normal = normal / norm_n

    # Centre du cercle par intersection des mediateurs
    d1 = np.dot(v1, v1)
    d2 = np.dot(v2, v2)
    d12 = np.dot(v1, v2)
    denom = 2 * (d1 * d2 - d12 ** 2)

    if abs(denom) < 1e-10:
        return move_l(q_start, pos_goal, v_cart, a_cart, n_points)

    alpha = d2 * (d1 - d12) / denom
    beta = d1 * (d2 - d12) / denom
    centre = p1 + alpha * v1 + beta * v2
    rayon = np.linalg.norm(p1 - centre)

    # Repere local du cercle
    u = (p1 - centre) / rayon
    w = normal
    v = np.cross(w, u)

    # Angles des trois points dans le plan du cercle
    angle_start = 0.0
    dp2 = p2 - centre
    angle_via = np.arctan2(np.dot(dp2, v), np.dot(dp2, u))
    dp3 = p3 - centre
    angle_end = np.arctan2(np.dot(dp3, v), np.dot(dp3, u))

    # Assurer la continuite de l'arc
    if angle_via < 0:
        angle_via += 2 * np.pi
    if angle_end < 0:
        angle_end += 2 * np.pi
    if angle_end < angle_via:
        angle_end += 2 * np.pi

    arc_length = rayon * angle_end
    t, s, _, _ = profil_trapeze(arc_length, v_cart, a_cart, n_points)

    trajectoire = np.zeros((n_points, cfg.NB_AXES))
    trajectoire[0] = q_start
    q_prev = q_start.copy()
    succes = True

    for i in range(1, n_points):
        ratio = s[i] / arc_length if arc_length > 0 else 0
        angle = angle_start + ratio * angle_end
        pos_inter = centre + rayon * (np.cos(angle) * u + np.sin(angle) * v)

        solutions = mgi_geometrique(pos_inter[0], pos_inter[1], pos_inter[2])
        q_sol = choisir_meilleure_solution(solutions, q_prev)

        if q_sol is None:
            q_sol, err, _, ok = mgi_numerique(pos_inter, q_init=q_prev)
            if not ok:
                succes = False
                q_sol = q_prev

        trajectoire[i] = q_sol
        q_prev = q_sol

    return trajectoire, t, succes


def comparer_trajectoires(q_start, q_goal, v_max=1.0, a_max=2.0, n_points=200):
    """Compare les trajectoires MoveJ (trapeze vs S-curve) et MoveL.

    Retourne un dictionnaire avec les metriques de chaque trajectoire.
    """
    pos_goal = position_effecteur_from_q(q_goal)

    traj_trap, t_trap = move_j(q_start, q_goal, v_max, a_max, 'trapeze', n_points)
    traj_scurve, t_scurve = move_j(q_start, q_goal, v_max, a_max, 'scurve', n_points)
    traj_lin, t_lin, ok_lin = move_l(q_start, pos_goal, 50.0, 200.0, n_points)

    def metriques_cartesiennes(traj):
        positions = np.array([mgd(q)[:3, 3] for q in traj])
        distances = np.linalg.norm(np.diff(positions, axis=0), axis=1)
        return {
            'longueur_totale': np.sum(distances),
            'ecart_max_lineaire': _ecart_lineaire_max(positions),
        }

    return {
        'trapeze': {
            'duree': t_trap[-1] if len(t_trap) > 0 else 0,
            'cartesien': metriques_cartesiennes(traj_trap),
        },
        'scurve': {
            'duree': t_scurve[-1] if len(t_scurve) > 0 else 0,
            'cartesien': metriques_cartesiennes(traj_scurve),
        },
        'lineaire': {
            'duree': t_lin[-1] if len(t_lin) > 0 else 0,
            'cartesien': metriques_cartesiennes(traj_lin),
            'succes_ik': ok_lin,
        },
    }


def position_effecteur_from_q(q):
    T = mgd(q)
    return T[:3, 3].copy()


def _ecart_lineaire_max(positions):
    """Ecart maximal de la trajectoire par rapport a la ligne droite."""
    if len(positions) < 2:
        return 0.0
    direction = positions[-1] - positions[0]
    longueur = np.linalg.norm(direction)
    if longueur < 1e-10:
        return 0.0
    direction = direction / longueur

    ecarts = []
    for p in positions:
        v = p - positions[0]
        projection = np.dot(v, direction)
        point_proj = positions[0] + projection * direction
        ecarts.append(np.linalg.norm(p - point_proj))

    return max(ecarts)
