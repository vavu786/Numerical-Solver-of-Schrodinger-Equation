// An example of solving a boundary-value problem via
// the shooting method. The Runge-Kutta and secant
// methods are used for integration and root search.

import java.lang.*;
public class Shooting {
    static final int n = 100, ni = 10, m = 5;
    static final double h = 1.0 / n;

    public static void main(String argv[]) {
        double del = 1e-6, alpha0 = 1, dalpha = 0.01;
        double y1[] = new double[n + 1];
        double y2[] = new double[n + 1];
        double y[] = new double[2];
        // Search for the proper solution of the equation
        y1[0] = y[0] = 0;
        y2[0] = y[1] = secant(ni, del, alpha0, dalpha);
        for (int i = 0; i < n; ++i) {
            double x = h * i;
            y = rungeKutta(y, x, h);
            y1[i + 1] = y[0];
            y2[i + 1] = y[1];
        }
        // Output the result in every m points
        for (int i = 0; i <= n; i += m)
            System.out.println(y1[i]);
    }

    public static double secant(int n, double del, double x, double dx) {
        int k = 0;
        double x1 = x + dx;
        while ((Math.abs(dx) > del) && (k < n)) {
            double d = f(x1) - f(x);
            double x2 = x1 - f(x1) * (x1 - x) / d;
            x = x1;
            x1 = x2;
            dx = x1 - x;
            k++;
        }
        if (k == n) System.out.println("Convergence not" +
                " found after" + n + " iterations");
        return x1;
    }

    public static double[] rungeKutta(double y[], double t, double dt) {
        int l = y.length;
        double c1[] = new double[l];
        double c2[] = new double[l];
        double c3[] = new double[l];
        double c4[] = new double[l];

        c1 = g(y, t);
        for (int i = 0; i < l; ++i)
            c2[i] = y[i] + dt * c1[i] / 2;

        c2 = g(c2, t + dt / 2);
        for (int i = 0; i < l; ++i)
            c3[i] = y[i] + dt * c2[i] / 2;

        c3 = g(c3, t + dt / 2);
        for (int i = 0; i < l; ++i)
            c4[i] = y[i] + dt * c3[i];

        c4 = g(c4, t + dt);

        for (int i = 0; i < l; ++i)
            c1[i] = y[i] + dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6;

        return c1;
    }

    // Method to provide the function for the root search.
    public static double f(double x) {
        double y[] = new double[2];
        y[0] = 0;
        y[1] = x;
        for (int i = 0; i < n - 1; ++i) {
            double xi = h * i;
            y = rungeKutta(y, xi, h);
        }
        return y[0] - 1;
    }

    // Method to provide the generalized velocity vector.
    public static double[] g(double y[], double t) {
        int k = y.length;
        double v[] = new double[k];
        v[0] = y[1];
        v[1] = -Math.PI*Math.PI*(y[0]+1)/4;
        return v;
    }
}