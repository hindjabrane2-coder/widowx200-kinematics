package kinematics;

// Parametres geometriques et limites articulaires du WidowX-200

public final class RobotConfig {

    private RobotConfig() {}

    // Dimensions DH (mm)
    public static final double L1 = 110.0;   // hauteur base -> epaule
    public static final double L2 = 200.0;   // bras (epaule -> coude)
    public static final double L3 = 200.0;   // avant-bras (coude -> poignet)
    public static final double L4 = 130.0;   // poignet -> extremite pince

    public static final double PORTEE_MAX = L2 + L3 + L4;
    public static final double PORTEE_MIN = 50.0;

    // Limites articulaires (radians)
    public static final double[] Q_MIN = { -Math.PI, -1.85, -2.15, -1.75 };
    public static final double[] Q_MAX = {  Math.PI,  1.60,  1.60,  2.10 };

    public static final String[] NOM_AXE = { "Base", "Epaule", "Coude", "Poignet" };

    // Vitesses max des Dynamixel XL430 (rad/s)
    public static final double VITESSE_MAX_SERVO = 5.76;

    // Parametres dynamiques par defaut
    public static final double VMAX = 1.5;       // rad/s
    public static final double AMAX = 3.0;        // rad/s^2
    public static final double JMAX = 30.0;       // rad/s^3
    public static final double DT   = 0.01;       // s

    // Inertie simplifiee par axe (kg.m^2)
    public static final double[] INERTIE_AXES = { 0.015, 0.045, 0.025, 0.008 };

    // Parametres cartesiens
    public static final double VITESSE_LINEAIRE = 100.0;  // mm/s
    public static final double ACCEL_LINEAIRE   = 500.0;  // mm/s^2
    public static final double SEUIL_SINGULARITE = 5e4;
    public static final double FACTEUR_RALENTISSEMENT = 0.2;
}
