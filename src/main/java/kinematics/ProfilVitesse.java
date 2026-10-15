package kinematics;

// Profils de vitesse trapezoidal et en S pour la planification de trajectoire

public class ProfilVitesse {

    // Resultat d'un profil mono-axe
    public double[] temps;
    public double[] position;
    public double[] vitesse;
    public double[] acceleration;
    public int longueur;

    // Profil trapezoidal classique

    public static ProfilVitesse trapeze(double qStart, double qEnd,
                                        double vmax, double amax, double dt) {
        ProfilVitesse p = new ProfilVitesse();
        double delta = qEnd - qStart;
        double dir = Math.signum(delta);
        double dist = Math.abs(delta);

        if (dist < 1e-6) {
            p.temps = new double[] { 0 };
            p.position = new double[] { qStart };
            p.vitesse = new double[] { 0 };
            p.acceleration = new double[] { 0 };
            p.longueur = 1;
            return p;
        }

        double tAcc = vmax / amax;
        double dAcc = 0.5 * amax * tAcc * tAcc;
        double tCruise, vPeak;

        if (2 * dAcc > dist) {
            tAcc = Math.sqrt(dist / amax);
            tCruise = 0;
            vPeak = amax * tAcc;
        } else {
            tCruise = (dist - 2 * dAcc) / vmax;
            vPeak = vmax;
        }

        double tTotal = 2 * tAcc + tCruise;
        int n = (int) (tTotal / dt) + 1;
        p.temps = new double[n];
        p.position = new double[n];
        p.vitesse = new double[n];
        p.acceleration = new double[n];
        p.longueur = n;

        double dAccReal = 0.5 * amax * tAcc * tAcc;

        for (int i = 0; i < n; i++) {
            double ti = i * dt;
            p.temps[i] = ti;

            if (ti <= tAcc) {
                p.position[i] = qStart + dir * 0.5 * amax * ti * ti;
                p.vitesse[i] = dir * amax * ti;
                p.acceleration[i] = dir * amax;
            } else if (ti <= tAcc + tCruise) {
                double tc = ti - tAcc;
                p.position[i] = qStart + dir * (dAccReal + vPeak * tc);
                p.vitesse[i] = dir * vPeak;
                p.acceleration[i] = 0;
            } else {
                double td = ti - tAcc - tCruise;
                p.position[i] = qStart + dir * (dAccReal + vPeak * tCruise
                                + vPeak * td - 0.5 * amax * td * td);
                p.vitesse[i] = dir * (vPeak - amax * td);
                p.acceleration[i] = -dir * amax;
            }
        }

        p.position[n - 1] = qEnd;
        p.vitesse[n - 1] = 0;
        p.acceleration[n - 1] = 0;
        return p;
    }

    // Profil en S (jerk limite) par integration numerique

    public static ProfilVitesse sCurve(double qStart, double qEnd,
                                       double vmax, double amax,
                                       double jmax, double dt) {
        ProfilVitesse p = new ProfilVitesse();
        double delta = qEnd - qStart;
        double dir = Math.signum(delta);
        double dist = Math.abs(delta);

        if (dist < 1e-6) {
            p.temps = new double[] { 0 };
            p.position = new double[] { qStart };
            p.vitesse = new double[] { 0 };
            p.acceleration = new double[] { 0 };
            p.longueur = 1;
            return p;
        }

        double tEstimate = 2 * dist / (0.5 * Math.min(vmax, 0.1)) + 4 * (amax / jmax);
        tEstimate = Math.max(tEstimate, 4 * (amax / jmax) + 0.1);
        int nMax = (int) (tEstimate / dt) + 1;

        double[] posArr = new double[nMax];
        double[] velArr = new double[nMax];
        double[] accArr = new double[nMax];
        double[] tArr   = new double[nMax];

        for (int i = 0; i < nMax; i++) tArr[i] = i * dt;

        String phase = "accel_up";
        int idxFin = nMax - 1;

        for (int i = 1; i < nMax; i++) {
            double v = velArr[i - 1];
            double a = accArr[i - 1];
            double distRestante = dist - posArr[i - 1];

            double dDecel = 0;
            if (v > 0) {
                double tDecel = v / amax + amax / jmax;
                dDecel = v * tDecel * 0.5;
            }

            double aNew, vNew;

            switch (phase) {
                case "accel_up":
                    aNew = a + jmax * dt;
                    if (aNew >= amax) { aNew = amax; phase = "accel_const"; }
                    vNew = v + aNew * dt;
                    if (vNew >= vmax) { vNew = vmax; phase = "cruise"; aNew = 0; }
                    if (distRestante <= dDecel + vNew * dt * 2) { phase = "decel_up"; aNew = 0; }
                    break;
                case "accel_const":
                    aNew = amax;
                    vNew = v + aNew * dt;
                    if (vNew >= vmax) { vNew = vmax; phase = "cruise"; aNew = 0; }
                    if (distRestante <= dDecel + vNew * dt * 2) { phase = "decel_up"; aNew = 0; }
                    break;
                case "cruise":
                    aNew = 0;
                    vNew = vmax;
                    if (distRestante <= dDecel + vNew * dt) { phase = "decel_up"; aNew = 0; }
                    break;
                case "decel_up":
                    aNew = a - jmax * dt;
                    if (aNew <= -amax) { aNew = -amax; phase = "decel_const"; }
                    vNew = v + aNew * dt;
                    break;
                case "decel_const":
                    aNew = -amax;
                    vNew = v + aNew * dt;
                    if (vNew <= 0.5 * amax * (amax / jmax)) { phase = "decel_down"; }
                    break;
                case "decel_down":
                    aNew = a + jmax * dt;
                    if (aNew >= 0) aNew = 0;
                    vNew = v + aNew * dt;
                    break;
                default:
                    aNew = 0;
                    vNew = 0;
            }

            vNew = Math.max(vNew, 0);
            double pNew = posArr[i - 1] + vNew * dt;

            if (pNew >= dist) {
                posArr[i] = dist;
                velArr[i] = 0;
                accArr[i] = 0;
                idxFin = i;
                break;
            }

            velArr[i] = vNew;
            accArr[i] = aNew;
            posArr[i] = pNew;
        }

        int n = idxFin + 1;
        p.temps = new double[n];
        p.position = new double[n];
        p.vitesse = new double[n];
        p.acceleration = new double[n];
        p.longueur = n;

        for (int i = 0; i < n; i++) {
            p.temps[i] = tArr[i];
            p.position[i] = qStart + dir * posArr[i];
            p.vitesse[i] = dir * velArr[i];
            p.acceleration[i] = dir * accArr[i];
        }
        return p;
    }

    // Duree d'un mouvement mono-axe (profil trapezoidal)

    public static double calculerDureeAxe(double delta, double vmax, double amax) {
        double d = Math.abs(delta);
        if (d < 1e-6) return 0;
        double tAcc = vmax / amax;
        double dAcc = 0.5 * amax * tAcc * tAcc;
        if (2 * dAcc > d) return 2 * Math.sqrt(d / amax);
        return 2 * tAcc + (d - 2 * dAcc) / vmax;
    }
}
