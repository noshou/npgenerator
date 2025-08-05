package Utilities;

import org.jetbrains.annotations.*;
import com.oson.tuple.*;
import org.apfloat.*;

/**
 * Utility class for 3D vector operations using {@link Apfloat} components.
 *
 * <p>This class provides static methods for common vector operations in 3D space,
 * including basic arithmetic operations, geometric calculations, and normal vector
 * computations. All operations use high-precision {@link Apfloat} arithmetic to
 * maintain numerical accuracy.</p>
 *
 * <p>Vectors are represented as {@link Triad} objects containing three {@link Apfloat}
 * components corresponding to x, y, and z coordinates.</p>
 *
 * @since 1.0
 * @author [Author Name]
 */
public final class VectorMath {

    /**
     * Private constructor to prevent instantiation of this utility class.
     */
    private VectorMath() {
        throw new AssertionError("Utility class should not be instantiated");
    }

    /**
     * Computes the Euclidean norm (magnitude) of a 3D vector.
     *
     * <p>The norm is calculated as √(x² + y² + z²) where x, y, z are the
     * vector components.</p>
     *
     * @param u the 3D vector with {@link Apfloat} components, must not be null
     * @return the norm as an {@link Apfloat}, always non-negative
     * @throws NullPointerException if {@code u} or any component is {@code null}
     * @see #normalize(Triad)
     */
    @Contract(pure = true)
    public static @NotNull Apfloat norm(@NotNull Triad<@NotNull Apfloat> u) {
        return ApfloatMath.sqrt(
                ApfloatMath.pow(u.fetch(0), 2)
                        .add(ApfloatMath.pow(u.fetch(1), 2))
                        .add(ApfloatMath.pow(u.fetch(2), 2))
        );
    }

    /**
     * Normalizes a 3D vector to a unit vector.
     *
     * <p>A unit vector has the same direction as the original vector but has
     * magnitude 1. The normalized vector is computed by dividing each component
     * by the vector's norm.</p>
     *
     * @param u the vector to normalize, must not be null and must have non-zero magnitude
     * @return a new unit vector in the same direction as {@code u}
     * @throws ArithmeticException if the norm is zero (cannot normalize zero vector)
     * @throws NullPointerException if {@code u} or any component is {@code null}
     * @see #norm(Triad)
     */
    @Contract(pure = true)
    public static @NotNull Triad<@NotNull Apfloat> normalize(
            @NotNull Triad<@NotNull Apfloat> u
    ) {
        Apfloat n = norm(u);
        return new Triad<>(
                u.fetch(0).divide(n),
                u.fetch(1).divide(n),
                u.fetch(2).divide(n)
        );
    }

    /**
     * Computes the dot product of two 3D vectors.
     *
     * <p>The dot product is calculated as u·v = u_x*v_x + u_y*v_y + u_z*v_z.
     * The result is a scalar that represents the projection of one vector onto
     * another, scaled by their magnitudes.</p>
     *
     * <p>Geometric interpretation:
     * <ul>
     *   <li>Positive result: vectors point in generally the same direction</li>
     *   <li>Zero result: vectors are perpendicular</li>
     *   <li>Negative result: vectors point in generally opposite directions</li>
     * </ul></p>
     *
     * @param u the first vector, must not be null
     * @param v the second vector, must not be null
     * @return the dot product as an {@link Apfloat}
     * @throws NullPointerException if any argument or component is {@code null}
     * @see #cross_prod(Triad, Triad)
     */
    @Contract(pure = true)
    public static @NotNull Apfloat dot_prod(@NotNull Triad<@NotNull Apfloat> u,
                                            @NotNull Triad<@NotNull Apfloat> v) {
        return u.fetch(0).multiply(v.fetch(0))
                .add(u.fetch(1).multiply(v.fetch(1)))
                .add(u.fetch(2).multiply(v.fetch(2)));
    }

    /**
     * Multiplies each component of a vector by a scalar.
     *
     * <p>This operation scales the vector by the given factor. The direction
     * remains unchanged unless the scalar is negative (which reverses direction)
     * or zero (which produces the zero vector).</p>
     *
     * @param u the vector to scale, must not be null
     * @param cons the scalar factor as a string representation, must not be null
     * @return a new vector scaled by {@code cons}
     * @throws NumberFormatException if {@code cons} is not a valid number format
     * @throws NullPointerException if any argument or vector component is {@code null}
     * @see #div(Triad, String)
     */
    @Contract(pure = true)
    public static @NotNull Triad<@NotNull Apfloat> mult(
            @NotNull Triad<@NotNull Apfloat> u,
            @NotNull String cons) {
        long p = u.fetch(0).precision();
        Apfloat c = new Apfloat(cons, p);
        return new Triad<>(
                u.fetch(0).multiply(c),
                u.fetch(1).multiply(c),
                u.fetch(2).multiply(c)
        );
    }

    /**
     * Divides each component of a vector by a scalar.
     *
     * <p>This operation scales the vector by the reciprocal of the given factor.
     * Equivalent to multiplying by 1/cons.</p>
     *
     * @param u the vector to scale, must not be null
     * @param cons the scalar divisor as a string representation, must not be null and not zero
     * @return a new vector scaled by 1/{@code cons}
     * @throws ArithmeticException if {@code cons} represents zero
     * @throws NumberFormatException if {@code cons} is not a valid number format
     * @throws NullPointerException if any argument or vector component is {@code null}
     * @see #mult(Triad, String)
     */
    @Contract(pure = true)
    public static @NotNull Triad<@NotNull Apfloat> div(
            @NotNull Triad<@NotNull Apfloat> u,
            @NotNull String cons) {
        long p = u.fetch(0).precision();
        Apfloat c = new Apfloat(cons, p);
        return new Triad<>(
                u.fetch(0).divide(c),
                u.fetch(1).divide(c),
                u.fetch(2).divide(c)
        );
    }

    /**
     * Computes the cross product of two 3D vectors.
     *
     * <p>The cross product u × v produces a vector that is orthogonal to both
     * input vectors. The magnitude of the result equals the area of the parallelogram
     * formed by the two vectors. The direction follows the right-hand rule.</p>
     *
     * <p>Mathematical formula: u × v = (u_y*v_z - u_z*v_y, u_z*v_x - u_x*v_z, u_x*v_y - u_y*v_x)</p>
     *
     * <p>Properties:
     * <ul>
     *   <li>Anti-commutative: u × v = -(v × u)</li>
     *   <li>Result is zero if vectors are parallel</li>
     *   <li>Magnitude equals |u||v|sin(θ) where θ is the angle between vectors</li>
     * </ul></p>
     *
     * @param u the first vector, must not be null
     * @param v the second vector, must not be null
     * @return a new vector orthogonal to both {@code u} and {@code v}
     * @throws NullPointerException if any argument or component is {@code null}
     * @see #dot_prod(Triad, Triad)
     */
    @Contract(pure = true)
    public static @NotNull Triad<@NotNull Apfloat> cross_prod(
            @NotNull Triad<@NotNull Apfloat> u,
            @NotNull Triad<@NotNull Apfloat> v
    ) {
        Apfloat u_x = u.fetch(0), u_y = u.fetch(1), u_z = u.fetch(2);
        Apfloat v_x = v.fetch(0), v_y = v.fetch(1), v_z = v.fetch(2);

        Apfloat n_x = u_y.multiply(v_z).subtract(u_z.multiply(v_y));
        Apfloat n_y = u_z.multiply(v_x).subtract(u_x.multiply(v_z));
        Apfloat n_z = u_x.multiply(v_y).subtract(u_y.multiply(v_x));

        return new Triad<>(n_x, n_y, n_z);
    }

    /**
     * Computes the difference of two 3D vectors (vector subtraction).
     *
     * <p>Calculates u - v by subtracting corresponding components:
     * (u_x - v_x, u_y - v_y, u_z - v_z). The result represents the vector
     * from point v to point u when vectors are interpreted as position vectors.</p>
     *
     * @param u the minuend vector, must not be null
     * @param v the subtrahend vector, must not be null
     * @return a new vector representing u - v
     * @throws NullPointerException if any argument or component is {@code null}
     * @see #add(Triad, Triad)
     */
    @Contract(pure = true)
    public static @NotNull Triad<@NotNull Apfloat> subs(
            @NotNull Triad<@NotNull Apfloat> u,
            @NotNull Triad<@NotNull Apfloat> v
    ) {
        Apfloat u_x = u.fetch(0), u_y = u.fetch(1), u_z = u.fetch(2);
        Apfloat v_x = v.fetch(0), v_y = v.fetch(1), v_z = v.fetch(2);

        Apfloat n_x = u_x.subtract(v_x);
        Apfloat n_y = u_y.subtract(v_y);
        Apfloat n_z = u_z.subtract(v_z);
        return new Triad<>(n_x, n_y, n_z);
    }

    /**
     * Computes the sum of two 3D vectors (vector addition).
     *
     * <p>Calculates u + v by adding corresponding components:
     * (u_x + v_x, u_y + v_y, u_z + v_z). Vector addition follows the
     * parallelogram rule and is commutative and associative.</p>
     *
     * @param u the first vector, must not be null
     * @param v the second vector, must not be null
     * @return a new vector representing u + v
     * @throws NullPointerException if any argument or component is {@code null}
     * @see #subs(Triad, Triad)
     */
    @Contract(pure = true)
    public static @NotNull Triad<@NotNull Apfloat> add(
            @NotNull Triad<@NotNull Apfloat> u,
            @NotNull Triad<@NotNull Apfloat> v
    ) {
        Apfloat u_x = u.fetch(0), u_y = u.fetch(1), u_z = u.fetch(2);
        Apfloat v_x = v.fetch(0), v_y = v.fetch(1), v_z = v.fetch(2);

        Apfloat n_x = u_x.add(v_x);
        Apfloat n_y = u_y.add(v_y);
        Apfloat n_z = u_z.add(v_z);
        return new Triad<>(n_x, n_y, n_z);
    }

    /**
     * Computes a normal vector for a quadrilateral face.
     *
     * <p>This method calculates the normal vector for a quadrilateral face defined by
     * four vertices.  If out set to true, the method ensures
     * that the normal points outwards by checking against
     * the face centroid. </p>

     *
     * @param v0 the first vertex of the quadrilateral, must not be null
     * @param v1 the second vertex of the quadrilateral, must not be null
     * @param v2 the third vertex of the quadrilateral, must not be null
     * @param v3 the fourth vertex of the quadrilateral, must not be null
     * @param out sets whether the vector should be forced to point outwards or not
     * @return a unit normal vector pointing outward from the face
     * @throws NullPointerException if any vertex or component is {@code null}
     * @throws ArithmeticException if the vertices are collinear (cannot compute normal)
     * @see #cross_prod(Triad, Triad)
     * @see #normalize(Triad)
     */
    @Contract(pure = true)
    public static @NotNull Triad<@NotNull Apfloat> normalQuad(
            @NotNull Triad<@NotNull Apfloat> v0,
            @NotNull Triad<@NotNull Apfloat> v1,
            @NotNull Triad<@NotNull Apfloat> v2,
            @NotNull Triad<@NotNull Apfloat> v3,
            boolean out
    ) {

        // Constants for computation
        long precision = v0.fetch(0).precision();
        Apfloat FOUR = new Apfloat("4", precision);

        // Use first three vertices to compute normal
        Triad<Apfloat> u = subs(v1, v0);
        Triad<Apfloat> v = subs(v2, v0);

        Triad<Apfloat> n = cross_prod(u, v);
        Triad<Apfloat> n_norm = normalize(n);


        if (out) {
            // Compute face centroid
            Apfloat c_x = (v0.fetch(0).add(v1.fetch(0)).add(v2.fetch(0)).add(v3.fetch(0))).divide(FOUR);
            Apfloat c_y = (v0.fetch(1).add(v1.fetch(1)).add(v2.fetch(1)).add(v3.fetch(1))).divide(FOUR);
            Apfloat c_z = (v0.fetch(2).add(v1.fetch(2)).add(v2.fetch(2)).add(v3.fetch(2))).divide(FOUR);
            Triad<Apfloat> centroid = new Triad<>(c_x, c_y, c_z);
            Apfloat ZERO = new Apfloat("0", precision);
            Apfloat dot = dot_prod(n_norm, centroid);
            if (dot.compareTo(ZERO) < 0) {
                return mult(n_norm, "-1");
            }
        }
        return n_norm;
    }

    /**
     * Computes an outward-pointing normal vector for a triangular face.
     *
    <p>This method calculates the normal vector for a triangular face defined by
     * three vertices.  If out set to true, the method ensures
     * that the normal points outwards by checking against
     * the face centroid.</p>
     *
     * @param v0 the first vertex of the triangle, must not be null
     * @param v1 the second vertex of the triangle, must not be null
     * @param v2 the third vertex of the triangle, must not be null
     * @param out sets whether the vector should be forced to point outwards or not
     * @return a unit normal vector pointing outward from the face
     * @throws NullPointerException if any vertex or component is {@code null}
     * @throws ArithmeticException if the vertices are collinear (cannot compute normal)
     * @see #normalQuad(Triad, Triad, Triad, Triad)
     * @see #cross_prod(Triad, Triad)
     * @see #normalize(Triad)
     */
    @Contract(pure = true)
    public static @NotNull Triad<@NotNull Apfloat> normalTriple(
            @NotNull Triad<@NotNull Apfloat> v0,
            @NotNull Triad<@NotNull Apfloat> v1,
            @NotNull Triad<@NotNull Apfloat> v2,
            boolean out
    ) {

        // Edge vectors (ordered for CCW when viewed from outside)
        Triad<Apfloat> u = subs(v1, v0);
        Triad<Apfloat> v = subs(v2, v0);

        // Cross product u × v gives normal
        Triad<Apfloat> n = cross_prod(u, v);
        Triad<Apfloat> n_norm = normalize(n);


        if (out) {
            // Compute face centroid
            long precision = v0.fetch(0).precision();
            Apfloat ZERO = new Apfloat("0", precision);
            Apfloat THREE = new Apfloat("3", precision);
            Triad<Apfloat> sum = add(add(v0, v1), v2);
            Triad<Apfloat> centroid = new Triad<>(
                    sum.fetch(0).divide(THREE),
                    sum.fetch(1).divide(THREE),
                    sum.fetch(2).divide(THREE)
            );
            Apfloat dot = dot_prod(n_norm, centroid);
            if (dot.compareTo(ZERO) < 0) {
                return mult(n_norm, "-1");
            }
        }
        return n_norm;
    }
}