// A program to study the motion of a particle under an
// elastic force in one dimension through the simplest
// predictor-corrector scheme.
import java.lang.*;
public class Predictor_Corrector {
    static final int n = 100, j = 5;

    public static void main(String argv[]) {
        double x[] = new double[n + 1];
        double v[] = new double[n + 1];
        // Assign time step and initial position and velocity
        double dt = 2 * Math.PI / n;
        x[0] = 0;
        v[0] = 1;
        // Calculate other position and velocity recursively
        for (int i = 0; i < n; ++i) {
        // Predict the next position and velocity
            x[i + 1] = x[i] + v[i] * dt;
            v[i + 1] = v[i] - x[i] * dt;
        // Correct the new position and velocity
            x[i + 1] = x[i] + (v[i] + v[i + 1]) * dt / 2;
            v[i + 1] = v[i] - (x[i] + x[i + 1]) * dt / 2;
        }
        // Output the result in every j time steps
        double t = 0;
        double jdt = j * dt;
        for (int i = 0; i <= n; i += j) {
            System.out.println(t + " " + x[i] + " " + v[i]);
            t += jdt;
        }
    }
}