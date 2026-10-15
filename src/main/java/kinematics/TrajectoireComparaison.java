package kinematics;

import static kinematics.RobotConfig.*;

// Comparaison MoveJ / MoveL : metriques de performance, etude parametrique, echantillonnage

public class TrajectoireComparaison {

    // Metriques d'une trajectoire

    public static class Metriques {
        public double duree;
        public double longueurChemin;
        public double distanceDirecte;
        public double ratioChemin;
        public double deviationMax;
        public double vitesseCartMax;
        public double vitesseCartMoy;
        public double accelCartMax;
        public double vitesseArtMax;
        public double accelArtMax;
        public double lissage;
        public double energie;
        public int nPoints;
        public double[] deviations;
        public double[] vCart;
        public double[] aCart;
    }

    // Resultat d'une comparaison

    public static class ResultatComparaison {
        public Metriques moveJ;
        public Metriques moveL;
    }

    // Resultat d'une etude parametrique en vitesse

    public static class ResultatVitesse {
        public double vitesse;
        public double dureeJ, dureeL;
        public double devMaxJ, devMaxL;
        public double energieJ, energieL;
    }

    // Vitesse cartesienne instantanee

    public static double[] vitesseCartesienne(double[] t, double[][] cart) {
        double[] v = new double[t.length];
        for (int i = 1; i < t.length; i++) {
            double dti = t[i] - t[i - 1];
            if (dti > 0) {
                double dx = cart[i][0] - cart[i - 1][0];
                double dy = cart[i][1] - cart[i - 1][1];
                double dz = cart[i][2] - cart[i - 1][2];
                v[i] = Math.sqrt(dx * dx + dy * dy + dz * dz) / dti;
            }
        }
        return v;
    }

    // Acceleration cartesienne

    public static double[] accelerationCartesienne(double[] t, double[] vitesse) {
        double[] a = new double[t.length];
        for (int i = 1; i < t.length; i++) {
            double dti = t[i] - t[i - 1];
            if (dti > 0) {
                a[i] = (vitesse[i] - vitesse[i - 1]) / dti;
            }
        }
        return a;
    }

    // Longueur totale du chemin cartesien

    public static double longueurChemin(double[][] cart) {
        double l = 0;
        for (int i = 1; i < cart.length; i++) {
            double dx = cart[i][0] - cart[i - 1][0];
            double dy = cart[i][1] - cart[i - 1][1];
            double dz = cart[i][2] - cart[i - 1][2];
            l += Math.sqrt(dx * dx + dy * dy + dz * dz);
        }
        return l;
    }

    // Deviation par rapport a la droite ideale

    public static double[] deviationDroite(double[][] cart, double[] pStart, double[] pEnd) {
        double[] dir = {
            pEnd[0] - pStart[0], pEnd[1] - pStart[1], pEnd[2] - pStart[2]
        };
        double longueur = Mat.norme(dir);
        double[] devs = new double[cart.length];

        if (longueur < 1e-6) return devs;
        for (int i = 0; i < 3; i++) dir[i] /= longueur;

        for (int i = 0; i < cart.length; i++) {
            double vx = cart[i][0] - pStart[0];
            double vy = cart[i][1] - pStart[1];
            double vz = cart[i][2] - pStart[2];
            double dot = vx * dir[0] + vy * dir[1] + vz * dir[2];
            double ex = vx - dot * dir[0];
            double ey = vy - dot * dir[1];
            double ez = vz - dot * dir[2];
            devs[i] = Math.sqrt(ex * ex + ey * ey + ez * ez);
        }
        return devs;
    }

    // Vitesse articulaire maximale

    public static double[] vitesseArtMax(double[] t, double[][] q) {
        int nAxes = q[0].length;
        double[] vMax = new double[nAxes];
        for (int i = 1; i < t.length; i++) {
            double dti = t[i] - t[i - 1];
            if (dti > 0) {
                for (int j = 0; j < nAxes; j++) {
                    double v = Math.abs(q[i][j] - q[i - 1][j]) / dti;
                    if (v > vMax[j]) vMax[j] = v;
                }
            }
        }
        return vMax;
    }

    // Acceleration articulaire maximale

    public static double[] accelArtMax(double[] t, double[][] q) {
        int nAxes = q[0].length;
        double[] aMax = new double[nAxes];
        double[] vPrec = new double[nAxes];

        for (int i = 1; i < t.length; i++) {
            double dti = t[i] - t[i - 1];
            if (dti > 0) {
                double[] v = new double[nAxes];
                for (int j = 0; j < nAxes; j++) {
                    v[j] = (q[i][j] - q[i - 1][j]) / dti;
                }
                if (i > 1) {
                    for (int j = 0; j < nAxes; j++) {
                        double a = Math.abs(v[j] - vPrec[j]) / dti;
                        if (a > aMax[j]) aMax[j] = a;
                    }
                }
                vPrec = v;
            }
        }
        return aMax;
    }

    // Indice de lissage : variation moyenne de direction

    public static double indiceLissage(double[][] cart) {
        if (cart.length < 3) return 0;

        double variations = 0;
        int nValid = 0;

        for (int i = 1; i < cart.length - 1; i++) {
            double v1x = cart[i][0] - cart[i - 1][0];
            double v1y = cart[i][1] - cart[i - 1][1];
            double v1z = cart[i][2] - cart[i - 1][2];
            double v2x = cart[i + 1][0] - cart[i][0];
            double v2y = cart[i + 1][1] - cart[i][1];
            double v2z = cart[i + 1][2] - cart[i][2];

            double n1 = Math.sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
            double n2 = Math.sqrt(v2x * v2x + v2y * v2y + v2z * v2z);

            if (n1 > 1e-6 && n2 > 1e-6) {
                double cosAngle = (v1x * v2x + v1y * v2y + v1z * v2z) / (n1 * n2);
                cosAngle = Math.max(-1, Math.min(1, cosAngle));
                variations += Math.abs(Math.acos(cosAngle));
                nValid++;
            }
        }
        return nValid == 0 ? 0 : variations / nValid;
    }

    // Estimation de l'energie consommee (modele simplifie)

    public static double estimationEnergie(double[] t, double[][] q) {
        int nAxes = q[0].length;
        double[][] vitesses = new double[t.length][nAxes];
        double[][] accelerations = new double[t.length][nAxes];

        for (int i = 1; i < t.length; i++) {
            double dti = t[i] - t[i - 1];
            if (dti > 0) {
                for (int j = 0; j < nAxes; j++) {
                    vitesses[i][j] = (q[i][j] - q[i - 1][j]) / dti;
                    if (i > 1) {
                        accelerations[i][j] = (vitesses[i][j] - vitesses[i - 1][j]) / dti;
                    }
                }
            }
        }

        double energie = 0;
        for (int i = 1; i < t.length; i++) {
            double dti = t[i] - t[i - 1];
            double puissance = 0;
            for (int j = 0; j < nAxes; j++) {
                double couple = INERTIE_AXES[j] * accelerations[i][j];
                puissance += Math.abs(couple * vitesses[i][j]);
            }
            energie += puissance * dti;
        }
        return energie;
    }

    // Erreur d'echantillonnage

    public static class ResultatEchantillonnage {
        public double erreurMax;
        public double erreurMoyenne;
    }

    public static ResultatEchantillonnage erreurEchantillonnage(double[] tRef,
                                                                 double[][] cartRef,
                                                                 double freqEch) {
        double dtEch = 1.0 / freqEch;
        int nEch = (int) (tRef[tRef.length - 1] / dtEch) + 1;
        double[][] cartEch = new double[nEch][3];
        double[] tEch = new double[nEch];

        for (int i = 0; i < nEch; i++) {
            tEch[i] = i * dtEch;
            int idx = java.util.Arrays.binarySearch(tRef, tEch[i]);
            if (idx < 0) idx = -idx - 1;
            idx = Math.min(idx, tRef.length - 1);
            cartEch[i] = cartRef[idx].clone();
        }

        double errMax = 0, errSomme = 0;
        int count = 0;

        for (int i = 0; i < tRef.length; i++) {
            int idx = java.util.Arrays.binarySearch(tEch, tRef[i]);
            if (idx < 0) idx = -idx - 2;
            idx = Math.max(0, Math.min(idx, nEch - 2));

            double alpha = 0;
            if (tEch[idx + 1] > tEch[idx]) {
                alpha = (tRef[i] - tEch[idx]) / (tEch[idx + 1] - tEch[idx]);
            }
            alpha = Math.max(0, Math.min(1, alpha));

            double px = cartEch[idx][0] * (1 - alpha) + cartEch[idx + 1][0] * alpha;
            double py = cartEch[idx][1] * (1 - alpha) + cartEch[idx + 1][1] * alpha;
            double pz = cartEch[idx][2] * (1 - alpha) + cartEch[idx + 1][2] * alpha;

            double dx = px - cartRef[i][0];
            double dy = py - cartRef[i][1];
            double dz = pz - cartRef[i][2];
            double err = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (err > errMax) errMax = err;
            errSomme += err;
            count++;
        }

        ResultatEchantillonnage re = new ResultatEchantillonnage();
        re.erreurMax = errMax;
        re.erreurMoyenne = count > 0 ? errSomme / count : 0;
        return re;
    }

    // Calcul complet des metriques pour une trajectoire

    public static Metriques calculerMetriques(double[] t, double[][] cart,
                                               double[][] q, double[] pStart,
                                               double[] pEnd) {
        Metriques m = new Metriques();
        m.duree = t[t.length - 1];
        m.nPoints = t.length;
        m.distanceDirecte = Mat.distance(pStart, pEnd);
        m.longueurChemin = longueurChemin(cart);
        m.ratioChemin = m.distanceDirecte > 0 ? m.longueurChemin / m.distanceDirecte : 0;

        m.deviations = deviationDroite(cart, pStart, pEnd);
        m.deviationMax = 0;
        for (double d : m.deviations) if (d > m.deviationMax) m.deviationMax = d;

        m.vCart = vitesseCartesienne(t, cart);
        m.aCart = accelerationCartesienne(t, m.vCart);

        m.vitesseCartMax = 0;
        double sommeV = 0;
        int countV = 0;
        for (double v : m.vCart) {
            if (v > m.vitesseCartMax) m.vitesseCartMax = v;
            if (v > 0) { sommeV += v; countV++; }
        }
        m.vitesseCartMoy = countV > 0 ? sommeV / countV : 0;

        m.accelCartMax = 0;
        for (double a : m.aCart) {
            double aa = Math.abs(a);
            if (aa > m.accelCartMax) m.accelCartMax = aa;
        }

        double[] vArt = vitesseArtMax(t, q);
        m.vitesseArtMax = 0;
        for (double v : vArt) if (v > m.vitesseArtMax) m.vitesseArtMax = v;

        double[] aArt = accelArtMax(t, q);
        m.accelArtMax = 0;
        for (double a : aArt) if (a > m.accelArtMax) m.accelArtMax = a;

        m.lissage = indiceLissage(cart);
        m.energie = estimationEnergie(t, q);
        return m;
    }

    // Comparaison MoveJ / MoveL entre deux points

    public static ResultatComparaison comparer(double[] pStart, double[] pEnd, double phi) {
        double[] qStart = Cinematique.cinematiqueInverse(pStart[0], pStart[1], pStart[2],
                                                          phi, true);
        double[] qEnd = Cinematique.cinematiqueInverse(pEnd[0], pEnd[1], pEnd[2],
                                                        phi, true);

        // MoveJ
        MoveJ.Trajectoire trajJ = MoveJ.planifier(qStart, qEnd);
        double[][] cartJ = new double[trajJ.n][3];
        for (int i = 0; i < trajJ.n; i++) {
            cartJ[i] = Cinematique.cinematiqueDirecte(trajJ.positions[i]);
        }

        // MoveL
        MoveL.ResultatMoveL trajL = MoveL.movel(pStart, pEnd, phi, phi);

        ResultatComparaison rc = new ResultatComparaison();
        rc.moveJ = calculerMetriques(trajJ.temps, cartJ, trajJ.positions, pStart, pEnd);
        rc.moveL = calculerMetriques(trajL.temps, trajL.positionsCart,
                                      trajL.positionsArt, pStart, pEnd);
        return rc;
    }

    // Etude parametrique en vitesse

    public static ResultatVitesse[] etudeVitesse(double[] pStart, double[] pEnd,
                                                  double[] vitesses, double phi) {
        ResultatVitesse[] resultats = new ResultatVitesse[vitesses.length];

        for (int v = 0; v < vitesses.length; v++) {
            double vi = vitesses[v];
            double scale = vi / VITESSE_LINEAIRE;
            double vmaxJ = Math.min(VMAX * scale, VMAX * 2);
            double amaxJ = Math.min(AMAX * scale, AMAX * 2);

            double[] qStart = Cinematique.cinematiqueInverse(pStart[0], pStart[1], pStart[2],
                                                              phi, true);
            double[] qEnd = Cinematique.cinematiqueInverse(pEnd[0], pEnd[1], pEnd[2],
                                                            phi, true);

            MoveJ.Trajectoire trajJ = MoveJ.planifier(qStart, qEnd, vmaxJ, amaxJ, DT, "trapeze");
            double[][] cartJ = new double[trajJ.n][3];
            for (int i = 0; i < trajJ.n; i++) {
                cartJ[i] = Cinematique.cinematiqueDirecte(trajJ.positions[i]);
            }

            MoveL.ResultatMoveL trajL = MoveL.movel(pStart, pEnd, vi, DT, phi, phi,
                                                     ACCEL_LINEAIRE, null, true);

            double[] devJ = deviationDroite(cartJ, pStart, pEnd);
            double[] devL = deviationDroite(trajL.positionsCart, pStart, pEnd);

            ResultatVitesse rv = new ResultatVitesse();
            rv.vitesse = vi;
            rv.dureeJ = trajJ.temps[trajJ.n - 1];
            rv.dureeL = trajL.temps[trajL.n - 1];
            rv.devMaxJ = 0;
            for (double d : devJ) if (d > rv.devMaxJ) rv.devMaxJ = d;
            rv.devMaxL = 0;
            for (double d : devL) if (d > rv.devMaxL) rv.devMaxL = d;
            rv.energieJ = estimationEnergie(trajJ.temps, trajJ.positions);
            rv.energieL = estimationEnergie(trajL.temps, trajL.positionsArt);

            resultats[v] = rv;
        }
        return resultats;
    }

    // Scenarios de test

    public static double[][] scenarioCourt() {
        return new double[][] { { 350, 0, 150 }, { 370, 20, 140 } };
    }

    public static double[][] scenarioLong() {
        return new double[][] { { 400, -50, 110 }, { 300, 100, 200 } };
    }

    public static double[][] scenarioSingularite() {
        return new double[][] { { 450, 0, 110 }, { 350, 50, 150 } };
    }

    public static double[][] scenarioFrontiere() {
        return new double[][] { { 400, -50, 110 }, { 400, 50, 110 } };
    }
}
