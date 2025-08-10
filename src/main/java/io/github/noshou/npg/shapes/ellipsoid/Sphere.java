package io.github.noshou.npg.shapes.ellipsoid;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;

/**
 * Represents a <b>Sphere</b>.
 */
public class Sphere extends Shape {

    /**
     * Constructs a new {@code Shapes.Sphere} instance with the given parameters.
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
     * @param point_cart the coordinates of a point in Cartesian space, must not be null
     * @return {@code true} if the point is within or on the sphere boundary, {@code false} otherwise
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        // get x,y,z coords
        Apfloat x_cart = point_cart.fetch(0);
        Apfloat y_cart = point_cart.fetch(1);
        Apfloat z_cart = point_cart.fetch(2);


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