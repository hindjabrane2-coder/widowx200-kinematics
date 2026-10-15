import numpy as np
from . import robot_config as cfg


def matrice_dh(a, alpha, d, theta):
    """Matrice de transformation homogene 4x4 selon la convention DH."""
    ct = np.cos(theta)
    st = np.sin(theta)
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    return np.array([
        [ct, -st * ca,  st * sa, a * ct],
        [st,  ct * ca, -ct * sa, a * st],
        [0.0, sa,       ca,      d],
        [0.0, 0.0,      0.0,     1.0],
    ])


def mgd(q):
    """Modele geometrique direct. Retourne la matrice 4x4 de l'effecteur."""
    T = np.eye(4)
    for i in range(cfg.NB_AXES):
        a, alpha, d, offset = cfg.DH_TABLE[i]
        T = T @ matrice_dh(a, alpha, d, q[i] + offset)
    return T


def position_effecteur(q):
    """Position cartesienne de l'effecteur (x, y, z) en mm."""
    T = mgd(q)
    return T[:3, 3].copy()


def orientation_effecteur(q):
    """Matrice de rotation 3x3 de l'effecteur."""
    T = mgd(q)
    return T[:3, :3].copy()


def matrices_intermediaires(q):
    """Retourne les matrices T0i pour chaque articulation (i=0..NB_AXES)."""
    matrices = [np.eye(4)]
    T = np.eye(4)
    for i in range(cfg.NB_AXES):
        a, alpha, d, offset = cfg.DH_TABLE[i]
        T = T @ matrice_dh(a, alpha, d, q[i] + offset)
        matrices.append(T.copy())
    return matrices


def jacobienne(q, delta=1e-6):
    """Jacobienne geometrique 6xN par differences finies."""
    p0 = mgd(q)
    pos0 = p0[:3, 3]
    R0 = p0[:3, :3]

    J = np.zeros((6, cfg.NB_AXES))
    for i in range(cfg.NB_AXES):
        qp = q.copy()
        qp[i] += delta
        pp = mgd(qp)

        J[:3, i] = (pp[:3, 3] - pos0) / delta

        dR = pp[:3, :3] @ R0.T
        J[3, i] = (dR[2, 1] - dR[1, 2]) / (2 * delta)
        J[4, i] = (dR[0, 2] - dR[2, 0]) / (2 * delta)
        J[5, i] = (dR[1, 0] - dR[0, 1]) / (2 * delta)

    return J


def manipulabilite(q):
    """Indice de manipulabilite (Yoshikawa) : sqrt(det(J * J^T))."""
    J = jacobienne(q)
    Jp = J[:3, :]
    JJt = Jp @ Jp.T
    det = np.linalg.det(JJt)
    return np.sqrt(max(det, 0.0))


def est_singulier(q):
    """Detecte si la configuration est proche d'une singularite."""
    return manipulabilite(q) < cfg.SINGULARITY_THRESHOLD


# Cinematique inverse geometrique (solution analytique pour 4 DDL)

def mgi_geometrique(px, py, pz):
    """Cinematique inverse geometrique.
    Retourne une liste de solutions (q1, q2, q3, q4) atteignant (px, py, pz).
    L'orientation du poignet (q4) est calculee pour garder l'effecteur horizontal.
    """
    solutions = []

    # Axe 1 : rotation de base
    q1 = np.arctan2(py, px)

    r = np.sqrt(px ** 2 + py ** 2)
    z = pz - cfg.L1

    # Position du poignet (on recule de L4 le long de la direction horizontale)
    rw = r - cfg.L4
    zw = z

    d = np.sqrt(rw ** 2 + zw ** 2)

    # Verification de l'accessibilite
    if d > cfg.L2 + cfg.L3 or d < abs(cfg.L2 - cfg.L3):
        return solutions

    # Axe 3 : coude (deux solutions, coude haut et coude bas)
    cos_q3 = (d ** 2 - cfg.L2 ** 2 - cfg.L3 ** 2) / (2 * cfg.L2 * cfg.L3)
    cos_q3 = np.clip(cos_q3, -1.0, 1.0)

    for signe in [1, -1]:
        q3 = signe * np.arccos(cos_q3)

        # Axe 2 : epaule
        beta = np.arctan2(zw, rw)
        phi = np.arctan2(cfg.L3 * np.sin(q3), cfg.L2 + cfg.L3 * np.cos(q3))
        q2 = beta - phi

        # Axe 4 : poignet (effecteur horizontal)
        q4 = -(q2 + q3)

        q = np.array([q1, q2, q3, q4])

        if cfg.dans_limites(q):
            solutions.append(q)

    return solutions


# Cinematique inverse numerique (Newton-Raphson avec amortissement de Levenberg-Marquardt)

def mgi_numerique(cible_pos, cible_orient=None, q_init=None,
                  max_iter=200, tol=1e-3, lambda_init=0.01):
    """Cinematique inverse iterative par methode de Newton-Raphson amortie.

    cible_pos : position cible (x, y, z) en mm
    cible_orient : matrice de rotation 3x3 cible (optionnel)
    q_init : configuration initiale (si None, utilise la position de repos)
    """
    if q_init is None:
        q_init = cfg.HOME_POSITION.copy()

    q = q_init.copy()
    lam = lambda_init

    for it in range(max_iter):
        T = mgd(q)
        pos = T[:3, 3]

        erreur_pos = cible_pos - pos

        if cible_orient is not None:
            R = T[:3, :3]
            Re = cible_orient @ R.T
            trace = np.clip((np.trace(Re) - 1) / 2, -1, 1)
            angle = np.arccos(trace)
            if abs(angle) < 1e-10:
                erreur_orient = np.zeros(3)
            else:
                axe = np.array([
                    Re[2, 1] - Re[1, 2],
                    Re[0, 2] - Re[2, 0],
                    Re[1, 0] - Re[0, 1],
                ]) / (2 * np.sin(angle))
                erreur_orient = angle * axe
            erreur = np.concatenate([erreur_pos, erreur_orient])
        else:
            erreur = erreur_pos

        norme = np.linalg.norm(erreur[:3])
        if norme < tol:
            return q, norme, it + 1, True

        J = jacobienne(q)
        if cible_orient is None:
            J = J[:3, :]

        # Levenberg-Marquardt : dq = (J^T J + lambda I)^-1 J^T e
        JtJ = J.T @ J
        n = JtJ.shape[0]
        A = JtJ + lam * np.eye(n)
        b = J.T @ erreur

        try:
            dq = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            lam *= 10
            continue

        q_new = q + dq
        q_new = cfg.saturer_articulations(q_new)

        T_new = mgd(q_new)
        erreur_new = np.linalg.norm(cible_pos - T_new[:3, 3])

        if erreur_new < norme:
            q = q_new
            lam = max(lam / 2, 1e-6)
        else:
            lam *= 5
            q = q + 0.1 * dq
            q = cfg.saturer_articulations(q)

    T_final = mgd(q)
    erreur_finale = np.linalg.norm(cible_pos - T_final[:3, 3])
    return q, erreur_finale, max_iter, erreur_finale < tol


def choisir_meilleure_solution(solutions, q_actuelle=None):
    """Parmi les solutions de la MGI, choisit la plus proche de q_actuelle."""
    if not solutions:
        return None
    if q_actuelle is None:
        q_actuelle = cfg.HOME_POSITION

    distances = [np.linalg.norm(sol - q_actuelle) for sol in solutions]
    return solutions[np.argmin(distances)]
