import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;

/**
 * Represents a cubic shape constructed on a specified lattice.
 * <p>
 * Extends the abstract {@link Shape} class to implement cubic boundary checks.
 * </p>
 */
public class Cube extends Shape {

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
    public Cube(
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
     * Determines whether the given Cartesian coordinates are inside the cubic boundary.
     *
     * @param x_cart the x-coordinate in Cartesian space, must not be null
     * @param y_cart the y-coordinate in Cartesian space, must not be null
     * @param z_cart the z-coordinate in Cartesian space, must not be null
     * @return {@code true} if the point is within or on the cube boundary, {@code false} otherwise
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(
            @NotNull Apfloat x_cart,
            @NotNull Apfloat y_cart,
            @NotNull Apfloat z_cart
    ) {
        // noinspection SuspiciousNameCombination
        return      ApfloatMath.abs(x_cart).compareTo(super.getRadius()) <= 0
                &&  ApfloatMath.abs(y_cart).compareTo(super.getRadius()) <= 0
                &&  ApfloatMath.abs(z_cart).compareTo(super.getRadius()) <= 0;
    }
}