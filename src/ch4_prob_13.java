// Problem 4.13 in Introduction to Computational Physics
// Solving Schroedinger equation using bisection and 4th order Runge-Kutta method

import java.io.FileNotFoundException;
import java.lang.*;
import java.io.FileOutputStream;
import java.io.PrintWriter;

public class ch4_prob_13 {
    static final int n = 300; // Number of iterations for Runge-Kutta
    static final int ni = 10; // Number of iterations for secant
    static final int m = 1; // Every m indices, output of wavefunction is printed to the file
    static final double x1 = -10.0, x2 = 10.0, h = (x2-x1)/n; // To have n steps, step size in x must be h
    static double[] ul0 = new double[n + 1]; // phi integrated from the left to the turning point
    static double[] ul1 = new double[n + 1]; // first derivative of phi integrated from the left to the turning point
    static double[] ur0 = new double[n + 1]; // phi integrated from the right to the turning point
    static double[] ur1 = new double[n + 1]; // first derivative of phi integrated from the right to the turning point
    static double[] ul = new double[2]; // values of ul0[] and ul1[] at a particular point
    static double[] ur = new double[2]; // values of ul0[] and ul1[] at a particular point
    static double[] u = new double[n + 1]; // The actual wavefunction phi, using ul for the left and ur for the right

    public static void main(String[] args) throws FileNotFoundException {
        double del = 1e-6; // Tolerance for root search (below this value = 0)
        double e = 13.0; // Initial guess for eigenenergy
        double de = 0.1; // Initial change in e for root search

        /* The root search. Returns the eigenvalue. The function of the root search depends on phi left anf phi right.
         *  phi left and phi right depend on the Runge-Kutta algorithm to create them.
         * The Runge-Kutta depends on the generalized velocity vector g().
         * g() depends on the second derivative of phi, which depends on the eigenenergy.
         * Everything is called through secant(). */

        e = secant(ni, del, e, de);

        // Prints wavefunction data to new file
        PrintWriter w = new PrintWriter(new FileOutputStream("ch4prob13.txt"), true);
        double x = x1;
        double mh = m*h;
        for (int i = 0; i < n+1; i += m) {
            w.println(x+" "+ u[i]);
            x += mh;
        }

        // Output the eigenvalue obtained
        // 8.317 (6.2 - 10.6) , 14.621 (12.7 - 15.9), 19.726 (17.4 - 21.3), 24.292 (22.0 - 25.3)
        System.out.println("The eigenvalue: " + e);
        System.out.println();

//        double testEigen = 0.0;
//        double guessEigen = 0.0;
//        for (int i = 0; i < 50; i++) {
//            try {
//                testEigen = secant(ni, del, guessEigen, de);
//                System.out.println(testEigen);
//                guessEigen = guessEigen + 1.0;
//            }
//            catch (ArrayIndexOutOfBoundsException ex) {
//                System.out.println(i + ": Array index out of bounds");
//            }
//        }
    }

    // Method to complete one Runge-Kutta step.
    public static double[] rungeKutta(double[] u, double x, double e, double dx) {
        int l = u.length;
        double[] c1 = new double[l];
        double[] c2 = new double[l];
        double[] c3 = new double[l];
        double[] c4 = new double[l];

        c1 = g(u, x, e);
        for (int i=0; i<l; ++i)
            c2[i] = u[i] + dx *c1[i]/2;

        c2 = g(c2, x + dx /2, e);
        for (int i=0; i<l; ++i)
            c3[i] = u[i] + dx *c2[i]/2;

        c3 = g(c3, x + dx /2, e);
        for (int i=0; i<l; ++i)
            c4[i] = u[i] + dx *c3[i];

        c4 = g(c4, x + dx, e);

        for (int i=0; i<l; ++i)
            c1[i] = u[i] + dx *(c1[i]+2*(c2[i]+c3[i])+c4[i])/6;

        return c1;
    }

    // Function used for numerical integration (integral of y with a step of h)
    public static double simpson(double[] y, double h) {
        int n = y.length-1;
        double s0 = 0, s1 = 0, s2 = 0;
        for (int i=1; i<n; i+=2) {
            s0 += y[i];
            s1 += y[i-1];
            s2 += y[i+1];
        }
        double s = (s1+4*s0+s2)/3;
        // Add the last slice separately for an even n+1
        if ((n+1)%2 == 0)
            return h*(s+(5*y[n]+8*y[n-1]-y[n-2])/12);
        else
            return h*s;
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
                " found after " + n + " iterations");
        return x1;
    }

    public static double f(double e0) {
        // All of the code until the line is for using e0 to construct a normalized wavefunction phi
        // Initial values (u represents phi, so "ul" = "phi subscript l", so on)
        ul0[0] = ul[0] = 0; // Left phi starts at 0
        ul1[0] = ul[1] = 0.01; // Left phi derivative starts at 0.01

        ur0[0] = ur[0] = 0; // Right phi starts at 0
        ur1[0] = ur[1] = -0.01; // Right phi derivative starts at -0.01

        double[] u_squared = new double[n + 1]; // Square of the wavefunction phi. Used for normalization

        // Convert turning point to an index i as the upper bound for rungeKutta loop
        int turnIndex = turning_index(e0); // ONLY COMMENT FOR FINITE OR INFINITE WELL
        //int turnIndex = 190;
        System.out.println(x1 + (h*turnIndex));

        // Runge Kutta for ul, from -5 -> turning point + 1
        for (int i = 0; i < turnIndex+2; i++) {
            double x = x1 + i * h;
            ul = rungeKutta(ul, x, e0, h);
            ul0[i+1] = ul[0];
            ul1[i+1] = ul[1];
        }

        // Runge Kutta for ur, from 10 -> turning point - 1
        for (int i = n; i > turnIndex-2; i--) {
            double x = x2 - ((n-i) * h);
            ur = rungeKutta(ur, x, e0, -1 * h);
            ur0[i] = ur[0];
            ur1[i] = ur[1];
        }

        // Rescale
        for (int i = 0; i < turnIndex+1; i++) {
            ul0[i] *= ur0[turnIndex] / ul0[turnIndex];
        }

        // Normalize
        for (int i = 0; i < turnIndex; i++) {
            u_squared[i] = ul0[i] * ul0[i];
            u[i] = ul0[i];
        }
        for (int i = turnIndex; i < n+1; i++) {
            u_squared[i] = ur0[i] * ur0[i];
            u[i] = ur0[i];
        }
        double sum = simpson(u_squared, h);
        sum = Math.sqrt(sum);
        for (int i = 0; i < n+1; i++) u[i] /= sum; // Edits the global value u (which represents the wavefunction phi)

        return (ur0[turnIndex-1]-ur0[turnIndex+1]-ul0[turnIndex-1]+ul0[turnIndex+1]) / (2*h*ur0[turnIndex]);
    }
    // Linear Potential
    public static double V(double x) {
        if (x < 0.0 || x > 5.0) {
            return 50.0;
        }
        else {
            return 10.0 * x;
        }
    }

    // Book Potential
//    public static double V(double x) {
//        double alpha = 1, lambda = 4;
//        return alpha*alpha*lambda*(lambda-1)
//                *(0.5-1/Math.pow(cosh(alpha*x),2))/2;
//    }

    // Infinite well potential
//    public static double V(double x) {
//        if (x < 0 || x > 5.0) {
//            return 100.0;
//        }
//        else {
//            return 0.0;
//        }
//    }
    // Method to provide the hyperbolic cosine needed.
    public static double cosh(double x) {
        return (Math.exp(x)+Math.exp(-x))/2;
    }

    // Method to provide the generalized velocity vector. This is what actually uses the differential equation
    public static double[] g(double[] y, double x, double energy) {
        int l = y.length;
        double[] v = new double[l];
        v[0] = y[1];
        v[1] = (-2 * y[0]) * (energy - V(x));
        return v;
    }

    public static int turning_index(double energy) {
        double xr;
        int return_i = 0;
        for (int i = 0; i < n; i++) {
            xr = x1 + (i * h);
            if (xr > 0.0 && energy - V(xr) <= 0.0) {
                return_i = i;
                break;
            }
        }
        return return_i;
    }
}
