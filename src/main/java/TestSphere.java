import com.oson.tuple.*;
import org.jetbrains.annotations.*;

/**
 * {@code TestSphere} is a utility class used to test and benchmark the construction
 * of a {@link Sphere} shape from given lattice parameters and atomic basis.
 * <p>
 * This class is intended for programmatic use by test drivers or benchmarking tools.
 */
public class TestSphere {

    /**
     * Builds a {@link Sphere} structure using the provided parameters and returns
     * the elapsed build time in minutes.
     *
     * @param args          an array of exactly 4 string arguments:
     *                      <ul>
     *                          <li>{@code args[0]} — radius in specified units (e.g., {@code "10"})</li>
     *                          <li>{@code args[1]} — lattice constant (e.g., {@code "0.408"})</li>
     *                          <li>{@code args[2]} — lattice type (currently only {@code "FCC"})</li>
     *                          <li>{@code args[3]} — test name (used for file and structure metadata)</li>
     *                          <li>{@code args[4]} - units of radius of nanoparticle (e.g., {@code "nm"})</li>
     *                      </ul>
     * @param atoms         the atomic basis to use, must not be {@code null}
     * @param precision     numerical precision to use for {@code Apfloat} calculations (e.g., 100)
     * @return the elapsed build time in nanoseconds
     * @throws IllegalArgumentException if {@code args[2]} is not a recognized lattice type
     */
    @Contract(pure = true)
    public long runTest(
            @NotNull String @NotNull [] args,
            @NotNull Polyad<Atom> atoms,
            int precision
    ) {
        // Validate argument count
        if (args.length != 5) {
            System.err.println(
                    """
                    Usage:   TestSphere.runTest({<radius‑nm>, <lattice‑constant>, <lattice‑type>, <name>, <units>})
                    Example: TestSphere.runTest({10, 0.408, FCC, sphere_au, nanometer})
                    """
            );
            return -1;
        }

        // Determine lattice type from args
        LatticeType lattice_type;
        if (args[2].equalsIgnoreCase("FCC")) {
            lattice_type = LatticeType.FCC;
        } else {
            throw new IllegalArgumentException("Unsupported lattice type: " + args[2]);
        }

        // Time the build process
        long start_time = System.nanoTime();

        Shape test = new Sphere(
                args[0],                  // radius (string)
                args[4],                  // units
                lattice_type,             // lattice type
                precision,                // Apfloat precision
                atoms,                    // atom basis
                args[1],                  // lattice constant
                "sphere_test_" + args[3], // file name
                "sphere_test_" + args[3], // structure name
                "sphere_test_" + args[3]  // structure index
        );
        test.build(true);
        long end_time = System.nanoTime();
        return end_time - start_time;
    }
}
