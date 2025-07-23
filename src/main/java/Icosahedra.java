import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;
import java.util.*;

public class Icosahedra extends Shape {

    // 12 vertices of the icosahedra
    private final Dodecad<Tuple<Apfloat>> vertices;

    /**
     * Constructs a new {@code Icosahedra} instance with the given parameters.
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
    public Icosahedra (
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
        // constants
        Apfloat ZERO = new Apfloat("0", precision);
        Apfloat ONE = new Apfloat("1", precision);
        Apfloat TWO = new Apfloat ("2", precision);
        Apfloat FIVE = new Apfloat("5", precision);
        Apfloat gold_ratio = ONE.add(ApfloatMath.sqrt(FIVE)).divide(TWO);

        // basis points of vertices
        Stack<Triad<Apfloat>> unscaled_vertices = new Stack<>();
        unscaled_vertices.push(new Triad<>(ONE, gold_ratio, ZERO));                     // +1, +φ, 0
        unscaled_vertices.push(new Triad<>(ONE.negate(), gold_ratio, ZERO));            // -1, +φ, 0
        unscaled_vertices.push(new Triad<>(ONE.negate(), gold_ratio.negate(), ZERO));   // -1, -φ, 0
        unscaled_vertices.push(new Triad<>(ONE, gold_ratio.negate(), ZERO));            // +1, -φ, 0
        unscaled_vertices.push(new Triad<>(gold_ratio, ZERO, ONE));                     // +φ, 0, +1
        unscaled_vertices.push(new Triad<>(gold_ratio, ZERO, ONE.negate()));            // +φ, 0, -1
        unscaled_vertices.push(new Triad<>(gold_ratio.negate(), ZERO, ONE.negate()));   // -φ, 0, -1
        unscaled_vertices.push(new Triad<>(gold_ratio.negate(), ZERO, ONE));            // -φ, 0, +1
        unscaled_vertices.push(new Triad<>(ZERO, ONE, gold_ratio));                     // 0, +1, +φ
        unscaled_vertices.push(new Triad<>(ZERO, ONE.negate(), gold_ratio));            // 0, -1, +φ
        unscaled_vertices.push(new Triad<>(ZERO, ONE.negate(), gold_ratio.negate()));   // 0, -1, -φ
        unscaled_vertices.push(new Triad<>(ZERO, ONE, gold_ratio.negate()));            // 0, +1, -φ

        // dist from origin = sqrt(x²+y²+z²) = sqrt(φ²+1²) = sqrt(φ²+1)
        Apfloat dist = ApfloatMath.sqrt(ONE.add(ApfloatMath.pow(gold_ratio, 2)));
        Apfloat r = super.getRadius();
        Iterator<Triad<Apfloat>> basis_point = unscaled_vertices.iterator();
        ArrayList<Triad<Apfloat>> _vertices = new ArrayList<>();

        // point (p_x, p_y, p_z) = (r * (basis_x)/dist, r * (basis_y)/length, r * (basis_z)/dist)
        while (basis_point.hasNext()) {

            // get basis points
            Triad<Apfloat> _basis = basis_point.next();
            Apfloat basis_x = _basis.fetch(0);
            Apfloat basis_y = _basis.fetch(1);
            Apfloat basis_z = _basis.fetch(2);

            // calculate cartesian point
            Apfloat p_x = r.multiply(basis_x.divide(dist));
            Apfloat p_y = r.multiply(basis_y.divide(dist));
            Apfloat p_z = r.multiply(basis_z.divide(dist));

            // add to list of vertices
            Triad<Apfloat> vertex = new Triad<>(p_x, p_y, p_z);
            _vertices.add(vertex);
        }

        // post-cond: 12 vertices
        assert(_vertices.size() == 12);
        vertices = new Dodecad<>(
                _vertices.get(0),
                _vertices.get(1),
                _vertices.get(2),
                _vertices.get(3),
                _vertices.get(4),
                _vertices.get(5),
                _vertices.get(6),
                _vertices.get(7),
                _vertices.get(8),
                _vertices.get(9),
                _vertices.get(10),
                _vertices.get(11)
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
        return false;
    }
}
