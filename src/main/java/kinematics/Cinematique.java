package kinematics;

import static kinematics.RobotConfig.*;

// Cinematique directe et inverse, jacobienne, singularites, Newton-Raphson

public class Cinematique {

    // Matrice de transformation DH 4x4

    public static double[][] matriceDH(double theta, double d, double a, double alpha) {
        double ct = Math.cos(theta), st = Math.sin(theta);
        double ca = Math.cos(alpha), sa = Math.sin(alpha);
        return new double[][] {
            { ct, -st * ca,  st * sa, a * ct },
            { st,  ct * ca, -ct * sa, a * st },
            { 0,   sa,       ca,      d      },
            { 0,   0,        0,       1      }
        };
    }

    // Cinematique directe complete par multiplication de matrices

    public static double[][] cinematiqueDirecteMatrice(double[] q) {
        double[][] t01 = matriceDH(q[0], L1, 0, Math.PI / 2);
        double[][] t12 = matriceDH(q[1], 0, L2, 0);
        double[][] t23 = matriceDH(q[2], 0, L3, 0);
        double[][] t34 = matriceDH(q[3], 0, L4, 0);
        double[][] t04 = Mat.mult4x4(Mat.mult4x4(Mat.mult4x4(t01, t12), t23), t34);
        return t04;
    }

    // Cinematique directe (x, y, z)

    public static double[] cinematiqueDirecte(double[] q) {
        double q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];
        double r = L2 * Math.cos(q2) + L3 * Math.cos(q2 + q3) + L4 * Math.cos(q2 + q3 + q4);
        double z = L1 + L2 * Math.sin(q2) + L3 * Math.sin(q2 + q3) + L4 * Math.sin(q2 + q3 + q4);
        double x = r * Math.cos(q1);
        double y = r * Math.sin(q1);
        return new double[] { x, y, z };
    }

    // Orientation de la pince dans le plan vertical

    public static double orientationEffecteur(double[] q) {
        return q[1] + q[2] + q[3];
    }

    // Jacobienne analytique 3x4

    public static double[][] jacobienneAnalytique(double[] q) {
        double q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];
        double s1 = Math.sin(q1), c1 = Math.cos(q1);
        double s23  = Math.sin(q2 + q3);
        double c23  = Math.cos(q2 + q3);
        double s234 = Math.sin(q2 + q3 + q4);
        double c234 = Math.cos(q2 + q3 + q4);

        double r = L2 * Math.cos(q2) + L3 * c23 + L4 * c234;

        double drDq2 = -L2 * Math.sin(q2) - L3 * s23 - L4 * s234;
        double drDq3 = -L3 * s23 - L4 * s234;
        double drDq4 = -L4 * s234;

        double dzDq2 = L2 * Math.cos(q2) + L3 * c23 + L4 * c234;
        double dzDq3 = L3 * c23 + L4 * c234;
        double dzDq4 = L4 * c234;

        return new double[][] {
            { -r * s1, c1 * drDq2, c1 * drDq3, c1 * drDq4 },
            {  r * c1, s1 * drDq2, s1 * drDq3, s1 * drDq4 },
            {  0,      dzDq2,      dzDq3,      dzDq4       }
        };
    }

    // Jacobienne numerique par differences finies

    public static double[][] jacobienneNumerique(double[] q) {
        double eps = 1e-6;
        double[] ref = cinematiqueDirecte(q);
        double[][] jac = new double[3][4];
        for (int j = 0; j < 4; j++) {
            double[] qp = q.clone();
            qp[j] += eps;
            double[] pp = cinematiqueDirecte(qp);
            for (int i = 0; i < 3; i++) {
                jac[i][j] = (pp[i] - ref[i]) / eps;
            }
        }
        return jac;
    }

    // Manipulabilite de Yoshikawa : sqrt(det(J*Jt))

    public static double manipulabilite(double[] q) {
        double[][] j = jacobienneAnalytique(q);
        double[][] jjt = Mat.multABt(j, j, 3, 4);
        double det = Mat.det3x3(jjt);
        return Math.sqrt(Math.max(det, 0.0));
    }

    // Detection de singularite

    public static boolean estSingulier(double[] q, double seuil) {
        double[][] j = jacobienneAnalytique(q);
        double[][] jjt = Mat.multABt(j, j, 3, 4);
        double det = Mat.det3x3(jjt);
        return Math.abs(det) < seuil;
    }

    public static boolean estSingulier(double[] q) {
        return estSingulier(q, 1e-3);
    }

    // Verification espace de travail

    public static boolean dansEspaceTravail(double x, double y, double z) {
        double rh = Math.sqrt(x * x + y * y);
        double d3d = Math.sqrt(rh * rh + (z - L1) * (z - L1));
        if (d3d > L2 + L3 + L4) return false;
        if (d3d < Math.abs(L2 - L3)) return false;
        return z >= 0;
    }

    // Verification des limites articulaires

    public static boolean limiteOk(double[] q) {
        for (int i = 0; i < 4; i++) {
            if (q[i] < Q_MIN[i] - 1e-6 || q[i] > Q_MAX[i] + 1e-6)
                return false;
        }
        return true;
    }

    // Cinematique inverse geometrique

    public static double[] cinematiqueInverse(double x, double y, double z,
                                              double phi, boolean coudeHaut) {
        if (!dansEspaceTravail(x, y, z)) {
            throw new IllegalArgumentException(
                String.format("Position (%.1f, %.1f, %.1f) hors espace", x, y, z));
        }

        double q1 = Math.atan2(y, x);
        double rTotal = Math.sqrt(x * x + y * y);
        double xw = rTotal - L4 * Math.cos(phi);
        double zw = z - L1 - L4 * Math.sin(phi);
        double d = Math.sqrt(xw * xw + zw * zw);

        if (d > L2 + L3 || d < Math.abs(L2 - L3)) {
            throw new IllegalArgumentException(
                String.format("Position (%.1f, %.1f, %.1f) hors de l'espace atteignable "
                    + "(d_poignet=%.1f, bornes=[%.1f, %.1f])",
                    x, y, z, d, Math.abs(L2 - L3), L2 + L3));
        }

        double cosQ3 = (d * d - L2 * L2 - L3 * L3) / (2 * L2 * L3);
        cosQ3 = Math.max(-1.0, Math.min(1.0, cosQ3));

        double q3 = coudeHaut ? -Math.acos(cosQ3) : Math.acos(cosQ3);
        double alpha = Math.atan2(zw, xw);
        double beta  = Math.atan2(L3 * Math.sin(q3), L2 + L3 * Math.cos(q3));
        double q2 = alpha - beta;
        double q4 = phi - q2 - q3;

        double[] q = { q1, q2, q3, q4 };
        if (!limiteOk(q)) {
            throw new IllegalArgumentException("Articulation hors limites");
        }
        return q;
    }

    public static double[] cinematiqueInverse(double x, double y, double z) {
        return cinematiqueInverse(x, y, z, 0.0, true);
    }

    // Deux solutions (coude haut et bas)

    public static double[][] deuxSolutions(double x, double y, double z, double phi) {
        java.util.List<double[]> solutions = new java.util.ArrayList<>();
        for (boolean coudeH : new boolean[] { true, false }) {
            try {
                solutions.add(cinematiqueInverse(x, y, z, phi, coudeH));
            } catch (IllegalArgumentException ignored) {}
        }
        return solutions.toArray(new double[0][]);
    }

    // Cinematique inverse iterative (Newton-Raphson)

    public static double[] cinematiqueInverseIterative(double[] cible, double phi,
                                                       double[] qInit, int maxIter,
                                                       double tol, double amortissement) {
        double[] q;
        if (qInit != null) {
            q = qInit.clone();
        } else {
            q = new double[] { Math.atan2(cible[1], cible[0]), 0.5, -0.5, 0.0 };
        }

        for (int iter = 0; iter < maxIter; iter++) {
            double[] posActuelle = cinematiqueDirecte(q);
            double[] erreur = { cible[0] - posActuelle[0],
                                cible[1] - posActuelle[1],
                                cible[2] - posActuelle[2] };
            double norme = Math.sqrt(erreur[0] * erreur[0] + erreur[1] * erreur[1]
                                   + erreur[2] * erreur[2]);
            if (norme < tol) {
                if (!limiteOk(q)) {
                    for (int i = 0; i < 4; i++) {
                        q[i] = Math.max(Q_MIN[i], Math.min(Q_MAX[i], q[i]));
                    }
                    continue;
                }
                return q;
            }

            double[][] jac = jacobienneAnalytique(q);

            // pseudo-inverse avec amortissement de Levenberg-Marquardt
            double lambda = 0.01;
            double[][] jjt = Mat.multABt(jac, jac, 3, 4);
            for (int i = 0; i < 3; i++) jjt[i][i] += lambda;

            double[][] jjtInv = Mat.inv3x3(jjt);
            // J_pinv = Jt * (J*Jt + lambda*I)^-1
            double[][] jPinv = new double[4][3];
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 3; j++) {
                    double s = 0;
                    for (int k = 0; k < 3; k++) s += jac[k][i] * jjtInv[k][j];
                    jPinv[i][j] = s;
                }
            }

            double[] dq = new double[4];
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 3; j++) dq[i] += jPinv[i][j] * erreur[j];
            }

            for (int i = 0; i < 4; i++) q[i] += amortissement * dq[i];

            // corriger orientation
            double phiActuel = q[1] + q[2] + q[3];
            q[3] += (phi - phiActuel) * 0.3;

            for (int i = 0; i < 4; i++) {
                q[i] = Math.max(Q_MIN[i], Math.min(Q_MAX[i], q[i]));
            }
        }

        throw new IllegalArgumentException("Pas de convergence apres " + maxIter + " iterations");
    }

    public static double[] cinematiqueInverseIterative(double[] cible, double phi) {
        return cinematiqueInverseIterative(cible, phi, null, 200, 0.01, 0.5);
    }

    // Validation : erreur entre la solution et la cible

    public static double validerSolution(double[] q, double[] cible) {
        double[] pos = cinematiqueDirecte(q);
        double dx = pos[0] - cible[0], dy = pos[1] - cible[1], dz = pos[2] - cible[2];
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
}
