package kinematics;

import static kinematics.RobotConfig.*;

// Analyse de l'espace de travail : echantillonnage, volume Monte Carlo, dexterite, charge utile

public class EspaceTravail {

    // Point atteignable avec son indice de manipulabilite

    public static class PointEspace {
        public double x, y, z;
        public double manipulabilite;

        public PointEspace(double x, double y, double z, double w) {
            this.x = x; this.y = y; this.z = z; this.manipulabilite = w;
        }
    }

    // Statistiques de l'espace de travail

    public static class Statistiques {
        public int nPoints;
        public double xMin, xMax;
        public double yMin, yMax;
        public double zMin, zMax;
        public double porteeMax, porteeMin, porteeMoyenne, porteeEcartType;
        public double manipMax, manipMin, manipMoyenne;
    }

    // Echantillonnage uniforme de l'espace articulaire

    public static PointEspace[] calculerEspaceTravail(int nParAxe, double[] q4Values) {
        double[] q1Range = linspace(Q_MIN[0], Q_MAX[0], nParAxe);
        double[] q2Range = linspace(Q_MIN[1], Q_MAX[1], nParAxe);
        double[] q3Range = linspace(Q_MIN[2], Q_MAX[2], nParAxe);

        java.util.List<PointEspace> points = new java.util.ArrayList<>();

        for (double q4 : q4Values) {
            for (double q2 : q2Range) {
                for (double q3 : q3Range) {
                    for (double q1 : q1Range) {
                        double[] q = { q1, q2, q3, q4 };
                        double[] p = Cinematique.cinematiqueDirecte(q);
                        if (p[2] >= 0) {
                            double w = Cinematique.manipulabilite(q);
                            points.add(new PointEspace(p[0], p[1], p[2], w));
                        }
                    }
                }
            }
        }

        return points.toArray(new PointEspace[0]);
    }

    public static PointEspace[] calculerEspaceTravail(int nParAxe) {
        return calculerEspaceTravail(nParAxe, new double[] { 0.0 });
    }

    // Coupe horizontale a une hauteur donnee

    public static PointEspace[] coupeHauteur(PointEspace[] points, double zRef,
                                              double epaisseur) {
        java.util.List<PointEspace> coupe = new java.util.ArrayList<>();
        for (PointEspace p : points) {
            if (Math.abs(p.z - zRef) < epaisseur) {
                coupe.add(p);
            }
        }
        return coupe.toArray(new PointEspace[0]);
    }

    // Estimation du volume par methode Monte Carlo

    public static double estimerVolumeMonteCarlo(int nEchantillons) {
        double[] qMin = Q_MIN;
        double[] qMax = Q_MAX;

        java.util.Random rand = new java.util.Random(42);
        double[][] pointsValides = new double[nEchantillons][3];
        int valides = 0;

        for (int i = 0; i < nEchantillons; i++) {
            double[] q = new double[4];
            for (int j = 0; j < 4; j++) {
                q[j] = qMin[j] + rand.nextDouble() * (qMax[j] - qMin[j]);
            }
            double[] p = Cinematique.cinematiqueDirecte(q);
            if (p[2] >= 0) {
                pointsValides[valides] = p;
                valides++;
            }
        }

        if (valides < 10) return 0;

        // discretiser en voxels
        double resolution = 20;
        double xMin = Double.MAX_VALUE, xMax = -Double.MAX_VALUE;
        double yMin = Double.MAX_VALUE, yMax = -Double.MAX_VALUE;
        double zMin = Double.MAX_VALUE, zMax = -Double.MAX_VALUE;

        for (int i = 0; i < valides; i++) {
            if (pointsValides[i][0] < xMin) xMin = pointsValides[i][0];
            if (pointsValides[i][0] > xMax) xMax = pointsValides[i][0];
            if (pointsValides[i][1] < yMin) yMin = pointsValides[i][1];
            if (pointsValides[i][1] > yMax) yMax = pointsValides[i][1];
            if (pointsValides[i][2] < zMin) zMin = pointsValides[i][2];
            if (pointsValides[i][2] > zMax) zMax = pointsValides[i][2];
        }

        int nxBins = (int) ((xMax - xMin) / resolution) + 1;
        int nyBins = (int) ((yMax - yMin) / resolution) + 1;
        int nzBins = (int) ((zMax - zMin) / resolution) + 1;

        boolean[][][] occupe = new boolean[nxBins][nyBins][nzBins];
        int nVoxels = 0;

        for (int i = 0; i < valides; i++) {
            int ix = (int) ((pointsValides[i][0] - xMin) / resolution);
            int iy = (int) ((pointsValides[i][1] - yMin) / resolution);
            int iz = (int) ((pointsValides[i][2] - zMin) / resolution);
            ix = Math.min(ix, nxBins - 1);
            iy = Math.min(iy, nyBins - 1);
            iz = Math.min(iz, nzBins - 1);
            if (!occupe[ix][iy][iz]) {
                occupe[ix][iy][iz] = true;
                nVoxels++;
            }
        }

        return nVoxels * resolution * resolution * resolution;
    }

    // Carte de dexterite dans le plan XZ (q1 = 0)

    public static class PointDexterite {
        public double x, z, manipulabilite;

        public PointDexterite(double x, double z, double w) {
            this.x = x; this.z = z; this.manipulabilite = w;
        }
    }

    public static PointDexterite[] carteDexterite(int n) {
        double[] q2Range = linspace(Q_MIN[1], Q_MAX[1], n);
        double[] q3Range = linspace(Q_MIN[2], Q_MAX[2], n);
        int nq4 = Math.max(n / 3, 5);
        double[] q4Range = linspace(Q_MIN[3], Q_MAX[3], nq4);

        java.util.List<PointDexterite> points = new java.util.ArrayList<>();

        for (double q2 : q2Range) {
            for (double q3 : q3Range) {
                for (double q4 : q4Range) {
                    double[] q = { 0.0, q2, q3, q4 };
                    double[] p = Cinematique.cinematiqueDirecte(q);
                    if (p[2] >= 0) {
                        double w = Cinematique.manipulabilite(q);
                        points.add(new PointDexterite(p[0], p[2], w));
                    }
                }
            }
        }

        return points.toArray(new PointDexterite[0]);
    }

    // Ellipsoide de manipulabilite (valeurs singulieres de la jacobienne)

    public static class Ellipsoide {
        public double[] valeursSingulieres;  // demi-axes
        public double ratioMaxMin;
    }

    public static Ellipsoide ellipsoideManipulabilite(double[] q) {
        double[][] J = Cinematique.jacobienneAnalytique(q);

        // SVD par decomposition de J*J^T (3x3)
        double[][] jjt = Mat.multABt(J, J, 3, 4);

        // valeurs propres de la matrice 3x3 symetrique
        // on utilise une decomposition simplifiee
        double[] eigenvalues = valeursPropres3x3(jjt);

        Ellipsoide ell = new Ellipsoide();
        ell.valeursSingulieres = new double[3];
        for (int i = 0; i < 3; i++) {
            ell.valeursSingulieres[i] = Math.sqrt(Math.max(0, eigenvalues[i]));
        }

        // trier par ordre decroissant
        java.util.Arrays.sort(ell.valeursSingulieres);
        double tmp = ell.valeursSingulieres[0];
        ell.valeursSingulieres[0] = ell.valeursSingulieres[2];
        ell.valeursSingulieres[2] = tmp;

        ell.ratioMaxMin = ell.valeursSingulieres[2] > 1e-10
            ? ell.valeursSingulieres[0] / ell.valeursSingulieres[2]
            : Double.POSITIVE_INFINITY;

        return ell;
    }

    // Estimation de la charge utile maximale

    public static double chargeMaxEstimee(double[] q, double[] couplesMax) {
        double[][] J = Cinematique.jacobienneAnalytique(q);
        double[][] pinv = Mat.pseudoInverse3x4(J);

        // force a l'effecteur : F = (J^T)^{-1} * tau
        // approximation par la pseudo-inverse de J^T
        // Jt_pinv = (J * J^T)^{-1} * J
        double[][] jjt = Mat.multABt(J, J, 3, 4);
        double[][] jjtInv;
        try {
            jjtInv = Mat.inv3x3(jjt);
        } catch (ArithmeticException e) {
            return 0;
        }

        // Jt_pinv (3x4) = jjtInv (3x3) * J (3x4)
        double[] forcesMax = new double[3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                double jijv = 0;
                for (int k = 0; k < 3; k++) {
                    jijv += jjtInv[i][k] * J[k][j];
                }
                forcesMax[i] += jijv * couplesMax[j];
            }
        }

        double forceVerticale = Math.abs(forcesMax[2]);
        // convertir en masse (distances en mm, force en N)
        return forceVerticale / 9.81 * 1000;
    }

    public static double chargeMaxEstimee(double[] q) {
        return chargeMaxEstimee(q, new double[] { 1.5, 1.5, 1.5, 1.5 });
    }

    // Carte de charge utile dans le plan XZ

    public static PointDexterite[] carteChargeUtile(int n) {
        double[] q2Range = linspace(Q_MIN[1], Q_MAX[1], n);
        double[] q3Range = linspace(Q_MIN[2], Q_MAX[2], n);

        java.util.List<PointDexterite> points = new java.util.ArrayList<>();

        for (double q2 : q2Range) {
            for (double q3 : q3Range) {
                double[] q = { 0.0, q2, q3, 0.0 };
                double[] p = Cinematique.cinematiqueDirecte(q);
                if (p[2] >= 0 && p[0] > 10) {
                    double ch = chargeMaxEstimee(q);
                    if (ch > 0 && ch < 100) {
                        points.add(new PointDexterite(p[0], p[2], ch));
                    }
                }
            }
        }

        return points.toArray(new PointDexterite[0]);
    }

    // Statistiques de l'espace de travail

    public static Statistiques statistiques(PointEspace[] points) {
        Statistiques s = new Statistiques();
        s.nPoints = points.length;

        if (points.length == 0) return s;

        s.xMin = Double.MAX_VALUE; s.xMax = -Double.MAX_VALUE;
        s.yMin = Double.MAX_VALUE; s.yMax = -Double.MAX_VALUE;
        s.zMin = Double.MAX_VALUE; s.zMax = -Double.MAX_VALUE;
        s.manipMax = -Double.MAX_VALUE; s.manipMin = Double.MAX_VALUE;

        double sommePortee = 0, sommePortee2 = 0;
        double sommeManip = 0;

        for (PointEspace p : points) {
            if (p.x < s.xMin) s.xMin = p.x;
            if (p.x > s.xMax) s.xMax = p.x;
            if (p.y < s.yMin) s.yMin = p.y;
            if (p.y > s.yMax) s.yMax = p.y;
            if (p.z < s.zMin) s.zMin = p.z;
            if (p.z > s.zMax) s.zMax = p.z;

            double portee = Math.sqrt(p.x * p.x + p.y * p.y);
            sommePortee += portee;
            sommePortee2 += portee * portee;
            if (portee > s.porteeMax) s.porteeMax = portee;
            if (portee > 10 && portee < s.porteeMin) s.porteeMin = portee;

            if (p.manipulabilite > s.manipMax) s.manipMax = p.manipulabilite;
            if (p.manipulabilite > 0 && p.manipulabilite < s.manipMin) {
                s.manipMin = p.manipulabilite;
            }
            sommeManip += p.manipulabilite;
        }

        int n = points.length;
        s.porteeMoyenne = sommePortee / n;
        double variance = sommePortee2 / n - s.porteeMoyenne * s.porteeMoyenne;
        s.porteeEcartType = Math.sqrt(Math.max(0, variance));
        s.manipMoyenne = sommeManip / n;

        return s;
    }

    // Export du nuage de points en CSV

    public static String exporterCSV(PointEspace[] points, String fichier)
            throws java.io.IOException {
        try (java.io.PrintWriter pw = new java.io.PrintWriter(fichier)) {
            pw.println("x,y,z");
            for (PointEspace p : points) {
                pw.printf("%.2f,%.2f,%.2f%n", p.x, p.y, p.z);
            }
        }
        return fichier;
    }

    // Utilitaires internes

    private static double[] linspace(double start, double end, int n) {
        double[] result = new double[n];
        if (n == 1) {
            result[0] = (start + end) / 2;
            return result;
        }
        for (int i = 0; i < n; i++) {
            result[i] = start + i * (end - start) / (n - 1);
        }
        return result;
    }

    // Valeurs propres d'une matrice 3x3 symetrique (Cardano)

    private static double[] valeursPropres3x3(double[][] m) {
        double a = m[0][0], b = m[1][1], c = m[2][2];
        double d = m[0][1], e = m[0][2], f = m[1][2];

        double p1 = d * d + e * e + f * f;
        if (p1 < 1e-20) {
            // matrice diagonale
            double[] vals = { a, b, c };
            java.util.Arrays.sort(vals);
            return new double[] { vals[2], vals[1], vals[0] };
        }

        double q = (a + b + c) / 3.0;
        double p2 = (a - q) * (a - q) + (b - q) * (b - q) + (c - q) * (c - q) + 2 * p1;
        double p = Math.sqrt(p2 / 6.0);

        double[][] B = {
            { (a - q) / p, d / p, e / p },
            { d / p, (b - q) / p, f / p },
            { e / p, f / p, (c - q) / p }
        };
        double detB = Mat.det3x3(B);
        double r = detB / 2.0;
        r = Math.max(-1, Math.min(1, r));

        double phi = Math.acos(r) / 3.0;

        double eig1 = q + 2 * p * Math.cos(phi);
        double eig3 = q + 2 * p * Math.cos(phi + 2 * Math.PI / 3);
        double eig2 = 3 * q - eig1 - eig3;

        return new double[] { eig1, eig2, eig3 };
    }
}
