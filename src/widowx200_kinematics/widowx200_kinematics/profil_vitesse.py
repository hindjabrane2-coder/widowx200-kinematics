import numpy as np


def profil_trapeze(distance, v_max, a_max, n_points=200):
    """Profil de vitesse trapezoldal (acceleration constante).

    Retourne les tableaux temps, position, vitesse, acceleration.
    La distance peut etre negative (le profil est symetrise).
    """
    signe = np.sign(distance) if distance != 0 else 1.0
    d = abs(distance)

    # Temps d'acceleration pour atteindre v_max
    t_acc = v_max / a_max

    # Distance parcourue pendant l'acceleration
    d_acc = 0.5 * a_max * t_acc ** 2

    if 2 * d_acc >= d:
        # Profil triangulaire (v_max jamais atteinte)
        t_acc = np.sqrt(d / a_max)
        t_cruise = 0.0
        v_pic = a_max * t_acc
    else:
        d_cruise = d - 2 * d_acc
        t_cruise = d_cruise / v_max
        v_pic = v_max

    t_total = 2 * t_acc + t_cruise

    t = np.linspace(0, t_total, n_points)
    pos = np.zeros(n_points)
    vel = np.zeros(n_points)
    acc = np.zeros(n_points)

    for i, ti in enumerate(t):
        if ti <= t_acc:
            # Phase d'acceleration
            acc[i] = a_max
            vel[i] = a_max * ti
            pos[i] = 0.5 * a_max * ti ** 2
        elif ti <= t_acc + t_cruise:
            # Phase a vitesse constante
            dt = ti - t_acc
            acc[i] = 0.0
            vel[i] = v_pic
            pos[i] = d_acc + v_pic * dt
        else:
            # Phase de deceleration
            dt = ti - t_acc - t_cruise
            acc[i] = -a_max
            vel[i] = v_pic - a_max * dt
            pos[i] = d_acc + v_pic * t_cruise + v_pic * dt - 0.5 * a_max * dt ** 2

    return t, signe * pos, signe * vel, signe * acc


def profil_s_curve(distance, v_max, a_max, j_max, n_points=200):
    """Profil de vitesse en S (jerk limite, 7 phases).

    Le jerk est borne, ce qui produit une acceleration progressive
    et reduit les vibrations mecaniques.
    """
    signe = np.sign(distance) if distance != 0 else 1.0
    d = abs(distance)

    # Phase 1 : jerk positif (montee d'acceleration)
    t_j = a_max / j_max
    a_pic = j_max * t_j

    # Vitesse atteinte a la fin de la rampe d'acceleration
    v_fin_rampe = 0.5 * j_max * t_j ** 2

    if 2 * v_fin_rampe >= v_max:
        # Pas de palier d'acceleration, profil reduit
        t_j = np.sqrt(v_max / j_max)
        t_a = 0.0
        v_pic = j_max * t_j ** 2
    else:
        # Palier d'acceleration (jerk = 0)
        t_a = (v_max - 2 * v_fin_rampe) / a_pic
        v_pic = v_max

    # Distances des phases d'acceleration et deceleration
    d_phase_acc = _distance_phase_acc(t_j, t_a, j_max)
    d_phase_dec = d_phase_acc

    if 2 * d_phase_acc >= d:
        # Pas de phase a vitesse constante, reduire v_pic
        facteur = d / (2 * d_phase_acc) if d_phase_acc > 0 else 0
        t_j *= np.sqrt(facteur)
        t_a *= facteur
        t_cruise = 0.0
        d_phase_acc = d / 2
        d_phase_dec = d / 2
    else:
        d_cruise = d - 2 * d_phase_acc
        t_cruise = d_cruise / v_pic

    # Duree totale
    t_total = 2 * (2 * t_j + t_a) + t_cruise

    t = np.linspace(0, t_total, n_points)
    pos = np.zeros(n_points)
    vel = np.zeros(n_points)
    acc = np.zeros(n_points)

    # Calcul par integration numerique des 7 phases
    dt_sim = t_total / (n_points - 1) if n_points > 1 else t_total
    v_curr = 0.0
    p_curr = 0.0

    t1 = t_j                       # fin jerk+
    t2 = t_j + t_a                 # fin acc constante
    t3 = 2 * t_j + t_a            # fin jerk- (vitesse max atteinte)
    t4 = t3 + t_cruise             # fin vitesse constante
    t5 = t4 + t_j                  # debut deceleration constante
    t6 = t4 + t_j + t_a           # fin deceleration constante
    t7 = t_total                   # arret

    for i, ti in enumerate(t):
        if ti <= t1:
            jerk = j_max
        elif ti <= t2:
            jerk = 0.0
        elif ti <= t3:
            jerk = -j_max
        elif ti <= t4:
            jerk = 0.0
        elif ti <= t5:
            jerk = -j_max
        elif ti <= t6:
            jerk = 0.0
        else:
            jerk = j_max

        if i > 0:
            dt_i = t[i] - t[i - 1]
            acc[i] = acc[i - 1] + jerk * dt_i
            vel[i] = vel[i - 1] + acc[i - 1] * dt_i + 0.5 * jerk * dt_i ** 2
            pos[i] = pos[i - 1] + vel[i - 1] * dt_i + 0.5 * acc[i - 1] * dt_i ** 2

    return t, signe * pos, signe * vel, signe * acc


def _distance_phase_acc(t_j, t_a, j_max):
    """Distance parcourue pendant les phases d'acceleration (jerk+, palier, jerk-)."""
    a_pic = j_max * t_j
    d1 = (1.0 / 6) * j_max * t_j ** 3
    v1 = 0.5 * j_max * t_j ** 2
    d2 = v1 * t_a + 0.5 * a_pic * t_a ** 2
    v2 = v1 + a_pic * t_a
    d3 = v2 * t_j + 0.5 * a_pic * t_j ** 2 - (1.0 / 6) * j_max * t_j ** 3
    return d1 + d2 + d3


def interpolation_cubique(points, n_inter=50):
    """Spline cubique naturelle passant par les points donnes.

    points : array (N, M) ou N est le nombre de points et M la dimension
    Retourne un array (N_total, M) de points interpoles.
    """
    n = len(points)
    if n < 2:
        return points.copy()

    points = np.atleast_2d(points)
    dim = points.shape[1]

    # Parametre t (cumulatif, base sur la distance)
    t = np.zeros(n)
    for i in range(1, n):
        t[i] = t[i - 1] + np.linalg.norm(points[i] - points[i - 1])
    if t[-1] == 0:
        return points.copy()
    t /= t[-1]

    # Resolution du systeme tridiagonal pour chaque dimension
    resultat = []
    t_out = np.linspace(0, 1, (n - 1) * n_inter + 1)

    for d in range(dim):
        y = points[:, d]
        coeffs = _spline_tridiagonal(t, y)
        vals = _evaluer_spline(t, y, coeffs, t_out)
        resultat.append(vals)

    return np.column_stack(resultat)


def _spline_tridiagonal(t, y):
    """Resolution du systeme tridiagonal pour les coefficients de la spline."""
    n = len(t)
    h = np.diff(t)
    b = np.diff(y) / h

    # Systeme tridiagonal pour les moments M
    A = np.zeros((n, n))
    rhs = np.zeros(n)
    A[0, 0] = 1.0
    A[-1, -1] = 1.0

    for i in range(1, n - 1):
        A[i, i - 1] = h[i - 1]
        A[i, i] = 2 * (h[i - 1] + h[i])
        A[i, i + 1] = h[i]
        rhs[i] = 6 * (b[i] - b[i - 1])

    M = np.linalg.solve(A, rhs)
    return M


def _evaluer_spline(t, y, M, t_eval):
    """Evalue la spline cubique aux points t_eval."""
    n = len(t)
    h = np.diff(t)
    result = np.zeros(len(t_eval))

    for k, te in enumerate(t_eval):
        # Trouver l'intervalle
        i = np.searchsorted(t, te, side='right') - 1
        i = max(0, min(i, n - 2))

        dt = te - t[i]
        hi = h[i]
        result[k] = (y[i]
                     + dt * ((y[i + 1] - y[i]) / hi
                             - hi * (2 * M[i] + M[i + 1]) / 6)
                     + dt ** 2 * M[i] / 2
                     + dt ** 3 * (M[i + 1] - M[i]) / (6 * hi))

    return result
