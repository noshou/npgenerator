import com.oson.tuple.*;
import org.jetbrains.annotations.*;

/**
 * {@code TestDriverAu} is a standalone test driver for benchmarking gold (Au)
 * nanoparticle structures using an FCC lattice.
 * <p>
 * It constructs a 4-atom FCC basis for gold, sets lattice parameters,
 * and runs both {@link TestSphere} and {@link TestCube} simulations.
 * Execution times are formatted using {@code FormatExecTime}.
 */
public class TestDriverAu {

    /**
     * Creates a default FCC atomic basis for gold (Au) using four atoms
     * positioned at the standard fractional coordinates for an FCC unit cell.
     *
     * @return a {@link Polyad} containing the four FCC-positioned gold atoms
     */
    private static @NotNull Polyad<Atom> getBasis() {
        Atom _1 = new Atom(
                "Au",
                "1.44",
                new Triad<>("0", "0", "0"),
                0,
                100
        );
        Atom _2 = new Atom(
                "Au",
                "1.44",
                new Triad<>("0.5", "0.5", "0"),
                0,
                100
        );
        Atom _3 = new Atom(
                "Au",
                "1.44",
                new Triad<>("0.5", "0", "0.5"),
                0,
                100
        );
        Atom _4 = new Atom(
                "Au",
                "1.44",
                new Triad<>("0", "0.5", "0.5"),
                0,
                100
        );
        Atom[] atoms = new Atom[]{_1, _2, _3, _4};
        return new Polyad<>(atoms);
    }

    /**
     * Main method for running FCC-based gold nanoparticle tests.
     * It builds both a sphere and cube using the Au FCC basis and
     * reports their construction times.
     *
     * @param args should contain radius of nanoparticle in nanometers
     */
    public static void main(String @NotNull [] args) {

        // Create FCC basis set of gold atoms
        Polyad<Atom> basis = getBasis();

        // Lattice parameters
        String lattice_constant = "4.33";    // FCC lattice constant (in angstroms)
        String lattice_type = "FCC";         // Lattice type
        String test_id = "Au";               // Unique test name

        // Argument array in expected format for test methods
        System.out.println("\nNanoparticle Radius:  \t"+args[0]+" "+args[1]);
        String[] _args = new String[]{
                args[0],
                lattice_constant,
                lattice_type,
                test_id,
                args[1]
        };

        int precision = 100; // Apfloat numerical precision

        // TODO: Add unit test for FccUnitCell here

        // --- Sphere test ---
        TestSphere ts = new TestSphere();
        long ts_time = ts.runTest(_args, basis, precision);
        System.out.println(
                "Sphere execution time:\t" + FormatExecTime.formatDuration(ts_time)
        );

        // --- Cube test ---
        TestCube tc = new TestCube();
        long tc_time = tc.runTest(_args, basis, precision);
        System.out.println(
                "Cube execution time:  \t" + FormatExecTime.formatDuration(tc_time) + "\n"
        );
    }
}
