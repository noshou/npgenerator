import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;

/**
 * Represents a spherical shape constructed on a specified lattice.
 * <p>
 * Extends the abstract {@link Shape} class to implement spherical boundary checks.
 * * </p>
 */
public class Sphere extends Shape {

    /**
     * Constructs a new {@code Sphere} instance with the given parameters.
     *
     * @param radius          the radius of the sphere as a string representation of a number
     * @param radius_type     the type of radius (e.g., "angstrom", "nm")
     * @param lattice_type    the lattice type enumeration specifying the crystal lattice structure
     * @param precision       the precision level for Apfloat calculations
     * @param basis           the {@link Polyad} of {@link Atom} objects representing the atomic basis of the lattice
     * @param lattice_constant the lattice constant as a string representation of a number
     * @param file_name       the output file name associated with this shape
     * @param structure_name  the name of the structure represented by this shape
     * @param structure_index the structure index identifier
     */
    public Sphere(
            @NotNull String radius,
            @NotNull String radius_type,
            @NotNull LatticeType lattice_type,
            int precision,
            @NotNull Polyad<Atom> basis,
            @NotNull String lattice_constant,
            @NotNull String file_name,
            @NotNull String structure_name,
            @NotNull String structure_index
    ) {
        super(
                radius,
                radius_type,
                lattice_type,
                precision,
                basis,
                lattice_constant,
                file_name,
                structure_name,
                structure_index
        );
    }

    /**
     * Determines whether the given Cartesian coordinates are inside the spherical boundary.
     *
     * @param x_cart the x-coordinate in Cartesian space, must not be null
     * @param y_cart the y-coordinate in Cartesian space, must not be null
     * @param z_cart the z-coordinate in Cartesian space, must not be null
     * @return {@code true} if the point is within or on the sphere boundary, {@code false} otherwise
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(
            @NotNull Apfloat x_cart,
            @NotNull Apfloat y_cart,
            @NotNull Apfloat z_cart
    ) {
        // Calculate squared distance from origin: x² + y² + z²
        // noinspection SuspiciousNameCombination
        Apfloat distance_squared =
                ApfloatMath.pow(x_cart, 2)
                .add(ApfloatMath.pow(y_cart, 2)) // IDEA throwing warning for y_cart hence the noinspection
                .add(ApfloatMath.pow(z_cart, 2));

        // if distance_squared <= radius², point within sphere
        return distance_squared.compareTo(ApfloatMath.pow(super.getRadius(), 2)) <= 0;
    }
}