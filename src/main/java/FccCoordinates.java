import java.util.concurrent.ConcurrentLinkedQueue;
import com.oson.tuple.Triad;
import org.apfloat.*;

/**
 * {@code FccCoordinates} generates a grid of fractional coordinates based on a face-centered cubic (FCC) lattice.
 * <p>
 * It constructs a cubic region of lattice points (with half-step spacing) centered at the origin,
 * based on a specified lattice constant and nanoparticle radius.
 * Each grid point represents a potential unit cell location, but no atoms are stored or assigned yet.
 * </p>
 *
 * <p>
 * This implementation is thread-safe due to its use of {@link ConcurrentLinkedQueue}.
 * </p>
 */
public class FccCoordinates extends AtomicCoordinates {
    /**
     * helper function for constructor
     *
     * @param lattice_constant      the FCC lattice constant, in the same units as the radius (e.g., Ångstroms)
     * @param nanoparticle_radius   the target nanoparticle radius, in same units as the lattice constant
     * @param precision             the number of digits for Apfloat precision
     */
    private static ConcurrentLinkedQueue<Triad<String>> build(
            String lattice_constant,
            String nanoparticle_radius,
            int precision
    ) {
        ConcurrentLinkedQueue<Triad<String>> grid_positions = new ConcurrentLinkedQueue<>();

        // Set up arbitrary-precision floating point constants
        Apfloat LC = new Apfloat(lattice_constant, precision);
        Apfloat NR = new Apfloat(nanoparticle_radius, precision);
        Apfloat TWO = new Apfloat("2", precision);

        // Compute the discrete search radius: how far we go from the center in each direction
        int DR = ApfloatMath.ceil(NR.divide(LC)).multiply(TWO).intValue();

        // Loop through the full range in 3D half-step increments (0.5 unit)
        for (int i = -DR; i <= DR; i++) {
            for (int j = -DR; j <= DR; j++) {
                for (int k = -DR; k <= DR; k++) {
                    Apfloat x = new Apfloat(i, precision).divide(TWO);
                    Apfloat y = new Apfloat(j, precision).divide(TWO);
                    Apfloat z = new Apfloat(k, precision).divide(TWO);
                    grid_positions.add(new Triad<>(x.toString(), y.toString(), z.toString()));
                }
            }
        }

        return grid_positions;
    }

    /**
     * Constructs an FCC grid of 3D positions based on a lattice constant and particle radius.
     *
     * @param lattice_constant      the FCC lattice constant, in the same units as the radius (e.g., Ångstroms)
     * @param nanoparticle_radius   the target nanoparticle radius, in same units as the lattice constant
     * @param precision             the number of digits for Apfloat precision
     */
    public FccCoordinates(String lattice_constant, String nanoparticle_radius, int precision) {
        super(build(lattice_constant, nanoparticle_radius, precision));
    }
}
