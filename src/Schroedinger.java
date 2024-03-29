// An example of solving the eigenvalue problem of the
// one-dimensional Schroedinger equation via the secant
// and Numerov methods.
import java.lang.*;
import java.io.*;

public class Schroedinger {
    static final int nx = 500, m = 1, ni = 10;
    static final double x1 = -10, x2 = 10, h = (x2-x1)/nx;
    static int nr, nl;
    static double[] ul = new double[nx+1];
    static double[] ur = new double[nx+1];
    static double ql[] = new double[nx+1];
    static double qr[] = new double[nx+1];
    static double s[] = new double[nx+1];
    static double u[] = new double[nx+1];

    public static void main(String argv[]) throws FileNotFoundException {
        double del = 1e-6, e = 2.4, de = 0.1;
        // Find the eigenvalue via the secant search
        e = secant(ni, del, e, de);
        // Output the wavefunction to a file
        PrintWriter w = new PrintWriter(new FileOutputStream("wave.data"), true);
        double x = x1;
        double mh = m*h;
        for (int i=0; i<=nx; i+=m) {
            w.println(x+" "+ u[i]);
            x += mh;
        }
        // Output the eigenvalue obtained
        // 6.094, 12.840
        System.out.println("The eigenvalue: " + e);

    }
    public static double secant(int n, double del, double x, double dx) {
        int k = 0;
        double x1 = x+dx;
        while ((Math.abs(dx)>del) && (k<n)) {
            double d = f(x1)-f(x);
            double x2 = x1-f(x1)*(x1-x)/d;
            x = x1;
            x1 = x2;
            dx = x1-x;
            k++;
        }
        if (k==n) System.out.println("Convergence not" +
                " found after"+n+" iterations");
        return x1;
    }
    // Method to provide the function for the root search.
    public static double f(double x) {
        wave(x);
        double f0 = ur[nr-1]+ul[nl-1]-ur[nr-3]-ul[nl-3];
        return f0/(2*h*ur[nr-2]);
    }
    // Method to calculate the wavefunction.
    public static void wave(double energy) {
        double y[] = new double [nx+1];
        double u0 = 0, u1 = 0.01;
        // Set up function q(x) in the equation
        for (int i=0; i<=nx; ++i) {
            double x = x1+i*h;
            ql[i] = 2*(energy-v(x));
            qr[nx-i] = ql[i];
        }
        // Find the matching point at the right turning point
        int im = 0;
        for (int i=0; i<nx; ++i)
            if (((ql[i]*ql[i+1])<0) && (ql[i]>0)) im = i;

        //System.out.println("Turning point: " + (x1 + im*h) + " for eigenvalue guess " + energy);

        // Carry out the Numerov integrations
        nl = im+2;
        nr = nx-im+2;
        ul = numerov(nl, h, u0, u1, ql, s);
        ur = numerov(nr, h, u0, u1, qr, s);
        // Find the wavefunction on the left
        double ratio = ur[nr-2]/ul[im];
        for (int i=0; i<=im; ++i) {
            u[i] = ratio*ul[i];
            y[i] = u[i]*u[i];
        }
        // Find the wavefunction on the right
        for (int i=0; i<nr-1; ++i) {
            u[i+im] = ur[nr-i-2];
            y[i+im] = u[i+im]*u[i+im];
        }
        // Normalize the wavefunction
        double sum = simpson(y, h);
        sum = Math.sqrt(sum);
        for (int i=0; i<=nx; ++i) u[i] /= sum;
    }
    // Method to perform the Numerov integration.
    public static double[] numerov(int m, double h,
                                   double u0, double u1, double q[], double s[]) {
        double u[] = new double[m];
        u[0] = u0;
        u[1] = u1;
        double g = h*h/12;
        for (int i=1; i<m-1; ++i) {
            double c0 = 1+g*q[i-1];
            double c1 = 2-10*g*q[i];
            double c2 = 1+g*q[i+1];
            double d = g*(s[i+1]+s[i-1]+10*s[i]);
            u[i+1] = (c1*u[i]-c0*u[i-1]+d)/c2;
        }
        return u;
    }
    public static double simpson(double y[], double h) {
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
    // Method to provide the given potential in the problem.
//    public static double v(double x) {
//        double alpha = 1, lambda = 4;
//        return alpha*alpha*lambda*(lambda-1)
//                *(0.5-1/Math.pow(cosh(alpha*x),2))/2;
//    }


//    public static double v(double x) {
//        if (x < 0.0 || x > 5.0) {
//            return 50.0;
//        }
//        else {
//            return 10.0 * x;
//        }
//    }

    // Harmonic Oscillator
    public static double v(double x) {
        return 0.5 * x*x;
    }

//    public static double v(double x) {
//        if (x < 0.0 || x > 5.0) {
//            return 50.0;
//        }
//        else {
//            return 0;
//        }
//    }
    // Method to provide the hyperbolic cosine needed.
    public static double cosh(double x) {
        return (Math.exp(x)+Math.exp(-x))/2;
    }
}
