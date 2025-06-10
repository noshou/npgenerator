public class ScatteringVector {
}

//import java.util.Scanner;
//
///**
// * Represents the magnitude of the scattering vector q used in
// * small-angle X-ray scattering (SAXS) calculations.
// *
// * <p>The scattering vector is computed using the formula:
// * <pre>
// *     q = (4π / λ) * sin(θ / 2)
// * </pre>
// * where λ is the X-ray wavelength and θ is the full scattering angle.
// *
// * <p>This implementation assumes:
// * <ul>
// *   <li>Wavelength λ is provided in nanometers (nm)</li>
// *   <li>The resulting q is converted from nm⁻¹ to Å⁻¹ by multiplying by 0.1</li>
// * </ul>
// */
//public class ScatteringVector {
//    /** The computed scattering vector magnitude in Å⁻¹ */
//    private double value;
//
//    /**
//     * Constructs a scattering vector from wavelength and angle.
//     *
//     * @param lambda_nm    the wavelength of incident X-rays in nanometers (nm)
//     * @param theta_rads   the full scattering angle in radians
//     */
//    public Q(double lambda_nm, double theta_rads) {
//        this.value = ((4 * Math.PI) / lambda_nm) * Math.sin(theta_rads / 2) * 0.1;
//    }
//
//    /**
//     * Returns the magnitude of the scattering vector q in Å⁻¹.
//     *
//     * @return the computed q value
//     */
//    public double value() {
//        return this.value;
//    }
//
//    /**
//     * Returns string representation of q
//     *
//     * @return the computed q value
//     */
//    public String toString() {
//        StringBuilder Q = new StringBuilder();
//        Q.append(String.format("%f Å⁻¹"), this.value());
//        return Q.toString();
//    }
//    /**
//     * Main method to prompt the user for wavelength and angle, and compute q.
//     *
//     * @param args command-line arguments (not used)
//     */
//    public static void main(String[] args) {
//        Scanner scanner = new Scanner(System.in);
//
//        System.out.print("Enter X-ray wavelength (nm): ");
//        double lambda_nm = scanner.nextDouble();
//
//        System.out.print("Enter scattering angle θ (radians): ");
//        double theta_degrees = scanner.nextDouble();
//
//        Q scatteringVector = new Q(lambda_nm, theta_degrees);
//
//        System.out.printf("Scattering vector q = %.6f Å⁻¹\n", scatteringVector.value());
//
//        scanner.close();
//    }
//}
