package kinematics;

import static kinematics.RobotConfig.*;

// Deplacement lineaire cartesien, mouvement circulaire, analyse de deviation

public class MoveL {

    // Resultat d'une trajectoire cartesienne

    public static class ResultatMoveL {
        public double[]   temps;
        public double[][] positionsCart;    // [echantillon][xyz]
        public double[][] positionsArt;     // [echantillon][axe]
        public double[]   vitessesCart;
        public int n;
    }

    // Interpolation SLERP simplifiee pour un angle (gere le passage par pi)

    public static double interpolerSlerp(double phiStart, double phiEnd, double s) {
        double diff = phiEnd - phiStart;
        diff = ((diff + Math.PI) % (2 * Math.PI)) - Math.PI;
        if (diff < -Math.PI) diff += 2 * Math.PI;
        return phiStart + s * diff;
    }

    // Profil de vitesse trapezoidal cartesien (retourne le parametre s entre 0 et 1)

    public static double[][] profilVitesseLineaire(double distance, double vitesseMax,
                                                    double accelMax, double dt) {
        if (distance < 1e-3) {
            return new double[][] { { 0.0 }, { 0.0 } };
        }

        double tAcc = vitesseMax / accelMax;
        double dAcc = 0.5 * accelMax * tAcc * tAcc;
        double tCruise, vPeak;

        if (2 * dAcc > distance) {
            tAcc = Math.sqrt(distance / accelMax);
            tCruise = 0;
            vPeak = accelMax * tAcc;
        } else {
            tCruise = (distance - 2 * dAcc) / vitesseMax;
            vPeak = vitesseMax;
        }

        double tTotal = 2 * tAcc + tCruise;
        int nPts = (int) (tTotal / dt) + 1;
        double[] t = new double[nPts];
        double[] s = new double[nPts];

        double dAccReal = 0.5 * accelMax * tAcc * tAcc;

        for (int i = 0; i < nPts; i++) {
            double ti = i * dt;
            t[i] = ti;
            double d;
            if (ti <= tAcc) {
                d = 0.5 * accelMax * ti * ti;
            } else if (ti <= tAcc + tCruise) {
                d = dAccReal + vPeak * (ti - tAcc);
            } else {
                double td = ti - tAcc - tCruise;
                d = dAccReal + vPeak * tCruise + vPeak * td - 0.5 * accelMax * td * td;
            }
            s[i] = d / distance;
        }
        s[nPts - 1] = 1.0;
        return new double[][] { t, s };
    }

    // Decalage TCP (Tool Center Point)

    public static double[] appliquerTcp(double[] position, double phi, double[] tcpOffset) {
        double dx = tcpOffset[0], dz = tcpOffset[1];
        double q1 = Math.atan2(position[1], position[0]);
        double rOffset = dx * Math.cos(phi) - dz * Math.sin(phi);
        double zOffset = dx * Math.sin(phi) + dz * Math.cos(phi);

        double rTotal = Math.sqrt(position[0] * position[0]
                                + position[1] * position[1]) + rOffset;
        double xNew = rTotal * Math.cos(q1);
        double yNew = rTotal * Math.sin(q1);
        double zNew = position[2] + zOffset;
        return new double[] { xNew, yNew, zNew };
    }

    // MoveL : deplacement lineaire cartesien

    public static ResultatMoveL movel(double[] pStart, double[] pEnd,
                                       double vitesse, double dt,
                                       double phiStart, double phiEnd,
                                       double accel, double[] tcpOffset,
                                       boolean eviterSingularites) {
        double dx = pEnd[0] - pStart[0];
        double dy = pEnd[1] - pStart[1];
        double dz = pEnd[2] - pStart[2];
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

        if (distance < 1e-3) {
            ResultatMoveL r = new ResultatMoveL();
            r.temps = new double[] { 0 };
            r.positionsCart = new double[][] { pStart.clone() };
            double[] q = Cinematique.cinematiqueInverse(pStart[0], pStart[1], pStart[2],
                                                        phiStart, true);
            r.positionsArt = new double[][] { q };
            r.vitessesCart = new double[] { 0 };
            r.n = 1;
            return r;
        }

        double[][] profil = profilVitesseLineaire(distance, vitesse, accel, dt);
        double[] t = profil[0];
        double[] s = profil[1];
        int nPts = t.length;

        ResultatMoveL res = new ResultatMoveL();
        res.temps = t;
        res.positionsCart = new double[nPts][3];
        res.positionsArt = new double[nPts][4];
        res.vitessesCart = new double[nPts];
        res.n = nPts;

        double[] qPrec = null;

        for (int i = 0; i < nPts; i++) {
            double si = s[i];
            double[] p = {
                pStart[0] + si * (pEnd[0] - pStart[0]),
                pStart[1] + si * (pEnd[1] - pStart[1]),
                pStart[2] + si * (pEnd[2] - pStart[2])
            };
            double phiI = interpolerSlerp(phiStart, phiEnd, si);

            if (tcpOffset != null) {
                p = appliquerTcp(p, phiI, tcpOffset);
            }

            res.positionsCart[i] = p;

            if (i > 0 && t[i] > t[i - 1]) {
                double ddx = res.positionsCart[i][0] - res.positionsCart[i - 1][0];
                double ddy = res.positionsCart[i][1] - res.positionsCart[i - 1][1];
                double ddz = res.positionsCart[i][2] - res.positionsCart[i - 1][2];
                res.vitessesCart[i] = Math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
                                    / (t[i] - t[i - 1]);
            }

            try {
                double[] q = Cinematique.cinematiqueInverse(p[0], p[1], p[2], phiI, true);

                if (eviterSingularites && qPrec != null) {
                    double w = Cinematique.manipulabilite(q);
                    // on peut ajuster le profil en amont si w < SEUIL_SINGULARITE
                }

                qPrec = q;
                res.positionsArt[i] = q;
            } catch (IllegalArgumentException e) {
                if (qPrec != null) {
                    res.positionsArt[i] = qPrec.clone();
                } else {
                    throw new IllegalArgumentException("Point de depart hors espace : " + e.getMessage());
                }
            }
        }

        return res;
    }

    public static ResultatMoveL movel(double[] pStart, double[] pEnd) {
        return movel(pStart, pEnd, VITESSE_LINEAIRE, DT, 0.0, 0.0,
                     ACCEL_LINEAIRE, null, true);
    }

    public static ResultatMoveL movel(double[] pStart, double[] pEnd,
                                       double phiStart, double phiEnd) {
        return movel(pStart, pEnd, VITESSE_LINEAIRE, DT, phiStart, phiEnd,
                     ACCEL_LINEAIRE, null, true);
    }

    // MoveC : deplacement circulaire passant par trois points

    public static ResultatMoveL movec(double[] pStart, double[] pVia, double[] pEnd,
                                       double vitesse, double dt, double phi) {
        double[] v1 = {
            pVia[0] - pStart[0], pVia[1] - pStart[1], pVia[2] - pStart[2]
        };
        double[] v2 = {
            pEnd[0] - pStart[0], pEnd[1] - pStart[1], pEnd[2] - pStart[2]
        };

        // normale au plan (produit vectoriel)
        double[] normal = {
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        };
        double normN = Mat.norme(normal);

        if (normN < 1e-6) {
            // points alignes : fallback en MoveL
            return movel(pStart, pEnd, vitesse, dt, phi, phi, ACCEL_LINEAIRE, null, true);
        }

        for (int i = 0; i < 3; i++) normal[i] /= normN;

        // centre du cercle par intersection des mediatrices
        double[] mid1 = {
            (pStart[0] + pVia[0]) / 2,
            (pStart[1] + pVia[1]) / 2,
            (pStart[2] + pVia[2]) / 2
        };
        double[] mid2 = {
            (pVia[0] + pEnd[0]) / 2,
            (pVia[1] + pEnd[1]) / 2,
            (pVia[2] + pEnd[2]) / 2
        };

        double[] dir1 = produitVectoriel(v1, normal);
        double[] v3 = { pEnd[0] - pVia[0], pEnd[1] - pVia[1], pEnd[2] - pVia[2] };
        double[] dir2 = produitVectoriel(v3, normal);

        double nDir1 = Mat.norme(dir1);
        double nDir2 = Mat.norme(dir2);
        for (int i = 0; i < 3; i++) { dir1[i] /= nDir1; dir2[i] /= nDir2; }

        // systeme surdetermine : mid1 + t1*dir1 = mid2 + t2*dir2
        // resolution par moindres carres (3 equations, 2 inconnues)
        double[] b = { mid2[0] - mid1[0], mid2[1] - mid1[1], mid2[2] - mid1[2] };
        double t1 = resoudreMoindresCarres(dir1, dir2, b);
        double[] centre = {
            mid1[0] + t1 * dir1[0],
            mid1[1] + t1 * dir1[1],
            mid1[2] + t1 * dir1[2]
        };

        double rayon = Mat.distance(pStart, centre);

        // repere local du cercle
        double[] u = new double[3];
        for (int i = 0; i < 3; i++) u[i] = (pStart[i] - centre[i]) / rayon;
        double[] w = produitVectoriel(normal, u);

        // angles sur le cercle
        double thetaStart = angleSurCercle(pStart, centre, u, w);
        double thetaVia = angleSurCercle(pVia, centre, u, w);
        double thetaEnd = angleSurCercle(pEnd, centre, u, w);

        // normaliser les angles
        thetaVia = normaliserAngle(thetaVia, thetaStart);
        thetaEnd = normaliserAngle(thetaEnd, thetaStart);

        if (thetaEnd > thetaStart) {
            if (thetaVia < thetaStart || thetaVia > thetaEnd) {
                thetaEnd -= 2 * Math.PI;
            }
        } else {
            if (thetaVia > thetaStart || thetaVia < thetaEnd) {
                thetaEnd += 2 * Math.PI;
            }
        }

        double arcLength = Math.abs(thetaEnd - thetaStart) * rayon;
        double duree = arcLength / vitesse;
        int nPts = Math.max((int) (duree / dt), 2);

        ResultatMoveL res = new ResultatMoveL();
        res.temps = new double[nPts];
        res.positionsCart = new double[nPts][3];
        res.positionsArt = new double[nPts][4];
        res.vitessesCart = new double[nPts];
        res.n = nPts;

        for (int i = 0; i < nPts; i++) {
            double frac = (double) i / (nPts - 1);
            double theta = thetaStart + frac * (thetaEnd - thetaStart);
            res.temps[i] = frac * duree;

            double cosT = Math.cos(theta), sinT = Math.sin(theta);
            res.positionsCart[i][0] = centre[0] + rayon * (cosT * u[0] + sinT * w[0]);
            res.positionsCart[i][1] = centre[1] + rayon * (cosT * u[1] + sinT * w[1]);
            res.positionsCart[i][2] = centre[2] + rayon * (cosT * u[2] + sinT * w[2]);

            if (i > 0 && res.temps[i] > res.temps[i - 1]) {
                double ddx = res.positionsCart[i][0] - res.positionsCart[i - 1][0];
                double ddy = res.positionsCart[i][1] - res.positionsCart[i - 1][1];
                double ddz = res.positionsCart[i][2] - res.positionsCart[i - 1][2];
                res.vitessesCart[i] = Math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
                                    / (res.temps[i] - res.temps[i - 1]);
            }

            try {
                double[] q = Cinematique.cinematiqueInverse(
                    res.positionsCart[i][0], res.positionsCart[i][1],
                    res.positionsCart[i][2], phi, true);
                res.positionsArt[i] = q;
            } catch (IllegalArgumentException e) {
                if (i > 0) {
                    res.positionsArt[i] = res.positionsArt[i - 1].clone();
                } else {
                    throw e;
                }
            }
        }

        return res;
    }

    public static ResultatMoveL movec(double[] pStart, double[] pVia, double[] pEnd) {
        return movec(pStart, pVia, pEnd, VITESSE_LINEAIRE, DT, 0.0);
    }

    // Analyse de deviation par rapport a la droite ideale

    public static class StatsDeviation {
        public double[] deviations;
        public double max;
        public double moyenne;
        public double rms;
        public double ecartType;
    }

    public static StatsDeviation analyserDeviation(double[][] cart,
                                                    double[] pStart, double[] pEnd) {
        double[] direction = {
            pEnd[0] - pStart[0], pEnd[1] - pStart[1], pEnd[2] - pStart[2]
        };
        double longueur = Mat.norme(direction);

        StatsDeviation sd = new StatsDeviation();
        sd.deviations = new double[cart.length];

        if (longueur < 1e-6) {
            sd.max = 0; sd.moyenne = 0; sd.rms = 0; sd.ecartType = 0;
            return sd;
        }

        for (int i = 0; i < 3; i++) direction[i] /= longueur;

        double somme = 0, sommeCarres = 0;
        sd.max = 0;

        for (int i = 0; i < cart.length; i++) {
            double vx = cart[i][0] - pStart[0];
            double vy = cart[i][1] - pStart[1];
            double vz = cart[i][2] - pStart[2];
            double dot = vx * direction[0] + vy * direction[1] + vz * direction[2];
            double ex = vx - dot * direction[0];
            double ey = vy - dot * direction[1];
            double ez = vz - dot * direction[2];
            double dev = Math.sqrt(ex * ex + ey * ey + ez * ez);
            sd.deviations[i] = dev;
            somme += dev;
            sommeCarres += dev * dev;
            if (dev > sd.max) sd.max = dev;
        }

        int n = cart.length;
        sd.moyenne = somme / n;
        sd.rms = Math.sqrt(sommeCarres / n);
        double variance = sommeCarres / n - sd.moyenne * sd.moyenne;
        sd.ecartType = Math.sqrt(Math.max(0, variance));
        return sd;
    }

    // Calcul de vitesse cartesienne instantanee

    public static double[] calculerVitesseCartesienne(double[] t, double[][] cart) {
        double[] vitesse = new double[t.length];
        for (int i = 1; i < t.length; i++) {
            double dti = t[i] - t[i - 1];
            if (dti > 0) {
                double ddx = cart[i][0] - cart[i - 1][0];
                double ddy = cart[i][1] - cart[i - 1][1];
                double ddz = cart[i][2] - cart[i - 1][2];
                vitesse[i] = Math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz) / dti;
            }
        }
        return vitesse;
    }

    // Calcul d'acceleration cartesienne

    public static double[] calculerAccelerationCartesienne(double[] t, double[] vitesse) {
        double[] accel = new double[t.length];
        for (int i = 1; i < t.length; i++) {
            double dti = t[i] - t[i - 1];
            if (dti > 0) {
                accel[i] = (vitesse[i] - vitesse[i - 1]) / dti;
            }
        }
        return accel;
    }

    // Transformation repere base vers monde

    public static double[] repereBaseVersMonde(double[] position, double[] offsetMonde,
                                                double rotationZ) {
        double cz = Math.cos(rotationZ), sz = Math.sin(rotationZ);
        double x = cz * position[0] - sz * position[1] + offsetMonde[0];
        double y = sz * position[0] + cz * position[1] + offsetMonde[1];
        double z = position[2] + offsetMonde[2];
        return new double[] { x, y, z };
    }

    // Transformation repere monde vers base

    public static double[] repereMondeVersBase(double[] posMonde, double[] offsetMonde,
                                                double rotationZ) {
        double cz = Math.cos(rotationZ), sz = Math.sin(rotationZ);
        double px = posMonde[0] - offsetMonde[0];
        double py = posMonde[1] - offsetMonde[1];
        double pz = posMonde[2] - offsetMonde[2];
        double x = cz * px + sz * py;
        double y = -sz * px + cz * py;
        return new double[] { x, y, pz };
    }

    // Utilitaires internes

    private static double[] produitVectoriel(double[] a, double[] b) {
        return new double[] {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        };
    }

    private static double angleSurCercle(double[] p, double[] centre,
                                          double[] u, double[] w) {
        double vx = p[0] - centre[0];
        double vy = p[1] - centre[1];
        double vz = p[2] - centre[2];
        double dotU = vx * u[0] + vy * u[1] + vz * u[2];
        double dotW = vx * w[0] + vy * w[1] + vz * w[2];
        return Math.atan2(dotW, dotU);
    }

    private static double normaliserAngle(double a, double ref) {
        while (a - ref < -Math.PI) a += 2 * Math.PI;
        while (a - ref > Math.PI) a -= 2 * Math.PI;
        return a;
    }

    // Resolution moindres carres pour le systeme mid1 + t1*dir1 = mid2 + t2*dir2
    // Retourne t1

    private static double resoudreMoindresCarres(double[] dir1, double[] dir2, double[] b) {
        // A = [dir1 | -dir2], systeme A*[t1,t2]^T = b
        // Normal equations : A^T A x = A^T b
        double a11 = 0, a12 = 0, a22 = 0, b1 = 0, b2 = 0;
        for (int i = 0; i < 3; i++) {
            a11 += dir1[i] * dir1[i];
            a12 += dir1[i] * (-dir2[i]);
            a22 += (-dir2[i]) * (-dir2[i]);
            b1 += dir1[i] * b[i];
            b2 += (-dir2[i]) * b[i];
        }
        double det = a11 * a22 - a12 * a12;
        if (Math.abs(det) < 1e-15) return 0;
        return (a22 * b1 - a12 * b2) / det;
    }
}
