package kinematics;

import static kinematics.RobotConfig.*;

// Deplacement articulaire synchronise, multi-points, splines cubiques

public class MoveJ {

    // Resultat d'une trajectoire multi-axes
    public static class Trajectoire {
        public double[]   temps;
        public double[][] positions;      // [echantillon][axe]
        public double[][] vitesses;
        public double[][] accelerations;
        public int n;
    }

    // MoveJ synchronise sur tous les axes

    public static Trajectoire planifier(double[] qStart, double[] qEnd,
                                        double vmax, double amax,
                                        double dt, String profil) {
        int nAxes = qStart.length;
        double[] deltas = new double[nAxes];
        for (int j = 0; j < nAxes; j++) deltas[j] = Math.abs(qEnd[j] - qStart[j]);

        double[] durees = new double[nAxes];
        double tSync = 0;
        for (int j = 0; j < nAxes; j++) {
            durees[j] = ProfilVitesse.calculerDureeAxe(deltas[j], vmax, amax);
            if (durees[j] > tSync) tSync = durees[j];
        }

        if (tSync < 1e-6) {
            Trajectoire t = new Trajectoire();
            t.temps = new double[] { 0 };
            t.positions = new double[][] { qStart.clone() };
            t.vitesses = new double[][] { new double[nAxes] };
            t.accelerations = new double[][] { new double[nAxes] };
            t.n = 1;
            return t;
        }

        int nPts = (int) (tSync / dt) + 1;
        Trajectoire traj = new Trajectoire();
        traj.temps = new double[nPts];
        traj.positions = new double[nPts][nAxes];
        traj.vitesses = new double[nPts][nAxes];
        traj.accelerations = new double[nPts][nAxes];
        traj.n = nPts;

        for (int i = 0; i < nPts; i++) traj.temps[i] = i * dt;

        for (int j = 0; j < nAxes; j++) {
            if (deltas[j] < 1e-6) {
                for (int i = 0; i < nPts; i++) traj.positions[i][j] = qStart[j];
                continue;
            }

            // adapter la vitesse de cet axe pour respecter tSync
            double discr = tSync * tSync * amax * amax - 4 * deltas[j] * amax;
            double vj;
            if (discr > 0) {
                vj = (tSync * amax - Math.sqrt(discr)) / 2;
            } else {
                vj = Math.sqrt(deltas[j] * amax);
            }
            vj = Math.min(vj, vmax);
            vj = Math.max(vj, 0.01);

            ProfilVitesse pv;
            if ("s".equals(profil)) {
                pv = ProfilVitesse.sCurve(qStart[j], qEnd[j], vj, amax, JMAX, dt);
            } else {
                pv = ProfilVitesse.trapeze(qStart[j], qEnd[j], vj, amax, dt);
            }

            int len = Math.min(nPts, pv.longueur);
            for (int i = 0; i < len; i++) {
                traj.positions[i][j] = pv.position[i];
                traj.vitesses[i][j] = pv.vitesse[i];
                traj.accelerations[i][j] = pv.acceleration[i];
            }
            for (int i = len; i < nPts; i++) {
                traj.positions[i][j] = qEnd[j];
            }
        }

        return traj;
    }

    public static Trajectoire planifier(double[] qStart, double[] qEnd) {
        return planifier(qStart, qEnd, VMAX, AMAX, DT, "trapeze");
    }

    // Trajectoire multi-points avec concatenation

    public static Trajectoire multiPoints(double[][] viaPoints,
                                          double vmax, double amax, double dt) {
        int nSeg = viaPoints.length - 1;
        if (nSeg < 1) {
            Trajectoire t = new Trajectoire();
            t.temps = new double[] { 0 };
            t.positions = new double[][] { viaPoints[0].clone() };
            t.vitesses = new double[][] { new double[viaPoints[0].length] };
            t.n = 1;
            return t;
        }

        java.util.List<double[]> tList   = new java.util.ArrayList<>();
        java.util.List<double[]> posList = new java.util.ArrayList<>();
        java.util.List<double[]> velList = new java.util.ArrayList<>();
        double tOffset = 0;

        for (int seg = 0; seg < nSeg; seg++) {
            Trajectoire segTraj = planifier(viaPoints[seg], viaPoints[seg + 1],
                                            vmax, amax, dt, "trapeze");
            for (int i = 0; i < segTraj.n; i++) {
                tList.add(new double[] { segTraj.temps[i] + tOffset });
                posList.add(segTraj.positions[i]);
                velList.add(segTraj.vitesses[i]);
            }
            tOffset += segTraj.temps[segTraj.n - 1];
        }

        Trajectoire result = new Trajectoire();
        result.n = tList.size();
        result.temps = new double[result.n];
        result.positions = new double[result.n][];
        result.vitesses = new double[result.n][];
        for (int i = 0; i < result.n; i++) {
            result.temps[i] = tList.get(i)[0];
            result.positions[i] = posList.get(i);
            result.vitesses[i] = velList.get(i);
        }
        return result;
    }

    // Spline cubique naturelle dans l'espace articulaire

    public static Trajectoire splineCubique(double[][] viaPoints, double[] durees, double dt) {
        int nPts = viaPoints.length;
        int nAxes = viaPoints[0].length;
        int nSeg = nPts - 1;

        double[] tNoeuds = new double[nPts];
        for (int i = 0; i < nSeg; i++) {
            tNoeuds[i + 1] = tNoeuds[i] + durees[i];
        }

        int totalPts = (int) (tNoeuds[nPts - 1] / dt) + 1;
        Trajectoire traj = new Trajectoire();
        traj.temps = new double[totalPts];
        traj.positions = new double[totalPts][nAxes];
        traj.vitesses = new double[totalPts][nAxes];
        traj.n = totalPts;

        for (int i = 0; i < totalPts; i++) traj.temps[i] = i * dt;

        for (int ax = 0; ax < nAxes; ax++) {
            double[] y = new double[nPts];
            for (int k = 0; k < nPts; k++) y[k] = viaPoints[k][ax];

            double[] h = durees.clone();

            // systeme tridiagonal pour les derivees secondes (spline naturelle)
            double[] lower = new double[nPts];
            double[] diag  = new double[nPts];
            double[] upper = new double[nPts];
            double[] rhs   = new double[nPts];

            diag[0] = 1;
            diag[nPts - 1] = 1;
            for (int i = 1; i < nPts - 1; i++) {
                lower[i] = h[i - 1];
                diag[i]  = 2 * (h[i - 1] + h[i]);
                upper[i] = h[i];
                rhs[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
            }

            double[] c = Mat.resoudreTridiagonal(lower, diag, upper, rhs);

            for (int seg = 0; seg < nSeg; seg++) {
                double aCoeff = y[seg];
                double bCoeff = (y[seg + 1] - y[seg]) / h[seg]
                              - h[seg] * (2 * c[seg] + c[seg + 1]) / 3;
                double cCoeff = c[seg];
                double dCoeff = (c[seg + 1] - c[seg]) / (3 * h[seg]);

                for (int i = 0; i < totalPts; i++) {
                    double t = traj.temps[i];
                    if (t >= tNoeuds[seg] - 1e-9 && t <= tNoeuds[seg + 1] + 1e-9) {
                        double tl = t - tNoeuds[seg];
                        traj.positions[i][ax] = aCoeff + bCoeff * tl
                                              + cCoeff * tl * tl + dCoeff * tl * tl * tl;
                        traj.vitesses[i][ax] = bCoeff + 2 * cCoeff * tl
                                             + 3 * dCoeff * tl * tl;
                    }
                }
            }
        }
        return traj;
    }

    // Estimation des couples articulaires (modele simplifie)

    public static double[][] estimerCouples(double[][] accelerations,
                                            double[][] vitesses) {
        int n = accelerations.length;
        int nAxes = accelerations[0].length;
        double frot = 0.01;
        double[][] couples = new double[n][nAxes];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < nAxes; j++) {
                couples[i][j] = INERTIE_AXES[j] * accelerations[i][j]
                              + frot * vitesses[i][j];
            }
        }
        return couples;
    }

    // Verification des limites articulaires sur toute la trajectoire

    public static int verifierTrajectoire(double[][] positions) {
        int depassements = 0;
        for (double[] pos : positions) {
            for (int j = 0; j < pos.length; j++) {
                if (pos[j] < Q_MIN[j] - 0.01 || pos[j] > Q_MAX[j] + 0.01) {
                    depassements++;
                }
            }
            double[] cart = Cinematique.cinematiqueDirecte(pos);
            if (cart[2] < -10) depassements++;
        }
        return depassements;
    }
}
