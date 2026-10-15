package kinematics;

// Utilitaires matriciels sans dependances externes

public final class Mat {

    private Mat() {}

    // Multiplication 4x4

    public static double[][] mult4x4(double[][] a, double[][] b) {
        double[][] r = new double[4][4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    r[i][j] += a[i][k] * b[k][j];
        return r;
    }

    // Produit A * B^T pour A(m x n) et B(m x n) -> resultat m x m

    public static double[][] multABt(double[][] a, double[][] b, int m, int n) {
        double[][] r = new double[m][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                for (int k = 0; k < n; k++)
                    r[i][j] += a[i][k] * b[j][k];
        return r;
    }

    // Determinant 3x3

    public static double det3x3(double[][] m) {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }

    // Inverse 3x3

    public static double[][] inv3x3(double[][] m) {
        double det = det3x3(m);
        if (Math.abs(det) < 1e-15) {
            throw new ArithmeticException("Matrice singuliere");
        }
        double invDet = 1.0 / det;
        double[][] r = new double[3][3];
        r[0][0] =  (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * invDet;
        r[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * invDet;
        r[0][2] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invDet;
        r[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * invDet;
        r[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invDet;
        r[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * invDet;
        r[2][0] =  (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * invDet;
        r[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * invDet;
        r[2][2] =  (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * invDet;
        return r;
    }

    // Norme euclidienne d'un vecteur

    public static double norme(double[] v) {
        double s = 0;
        for (double x : v) s += x * x;
        return Math.sqrt(s);
    }

    // Distance entre deux vecteurs

    public static double distance(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < a.length; i++) {
            double d = a[i] - b[i];
            s += d * d;
        }
        return Math.sqrt(s);
    }

    // Pseudo-inverse par la methode de Moore-Penrose (via J^T * (J*J^T)^-1)
    // Pour une matrice 3x4

    public static double[][] pseudoInverse3x4(double[][] j) {
        double[][] jjt = multABt(j, j, 3, 4);
        double[][] jjtInv = inv3x3(jjt);
        double[][] pinv = new double[4][3];
        for (int i = 0; i < 4; i++)
            for (int col = 0; col < 3; col++)
                for (int k = 0; k < 3; k++)
                    pinv[i][col] += j[k][i] * jjtInv[k][col];
        return pinv;
    }

    // Resolution d'un systeme tridiagonal (Thomas algorithm)

    public static double[] resoudreTridiagonal(double[] lower, double[] diag,
                                               double[] upper, double[] rhs) {
        int n = diag.length;
        double[] c = new double[n];
        double[] d = new double[n];
        double[] x = new double[n];

        c[0] = upper[0] / diag[0];
        d[0] = rhs[0] / diag[0];
        for (int i = 1; i < n; i++) {
            double m = diag[i] - lower[i] * c[i - 1];
            c[i] = i < n - 1 ? upper[i] / m : 0;
            d[i] = (rhs[i] - lower[i] * d[i - 1]) / m;
        }
        x[n - 1] = d[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            x[i] = d[i] - c[i] * x[i + 1];
        }
        return x;
    }
}
