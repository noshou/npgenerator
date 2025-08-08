package Shapes;

import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;
import java.util.*;
import Lattice.*;
import Atom.*;
import Utilities.VectorMath;

@SuppressWarnings("FieldCanBeLocal")
public class Icosahedron extends Shape {

    // 20 triangle faces of the icosahedra
    private final Icosad<Tuple<Tuple<Apfloat>>> faces;

    // 20 normalized face vectors of the icosahedra
    private final Icosad<Tuple<Apfloat>> face_norms;

    // constants
    private final Apfloat NEG_N1 = new Apfloat("-1", precision);
    private final Apfloat N0 = new Apfloat("0", precision);
    private final Apfloat N1 = new Apfloat("1", precision);
    private final Apfloat N2 = new Apfloat ("2", precision);
    private final Apfloat N3 = new Apfloat ("3", precision);
    private final Apfloat N5 = new Apfloat("5", precision);


    /**
     * Constructs a new {@code Shapes.Icosahedron} instance with the given parameters.
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
    public Icosahedron(
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
        Apfloat gold_ratio = N1.add(ApfloatMath.sqrt(N5)).divide(N2);

        // indexed vertex basis points
        Triad<Apfloat> vB0 = new Triad<>(N1, gold_ratio, N0);                                 // +1, +φ, 0
        Triad<Apfloat> vB1 = new Triad<>(N1.multiply(NEG_N1), gold_ratio, N0);                   // -1, +φ, 0
        Triad<Apfloat> vB2 = new Triad<>(N1.multiply(NEG_N1), gold_ratio.multiply(NEG_N1), N0);     // -1, -φ, 0
        Triad<Apfloat> vB3 = new Triad<>(N1, gold_ratio.multiply(NEG_N1), N0);                   // +1, -φ, 0
        Triad<Apfloat> vB4 = new Triad<>(gold_ratio, N0, N1);                                 // +φ, 0, +1
        Triad<Apfloat> vB5 = new Triad<>(gold_ratio, N0, N1.multiply(NEG_N1));                   // +φ, 0, -1
        Triad<Apfloat> vB6 = new Triad<>(gold_ratio.multiply(NEG_N1), N0, N1.multiply(NEG_N1));     // -φ, 0, -1
        Triad<Apfloat> vB7 = new Triad<>(gold_ratio.multiply(NEG_N1), N0, N1);                   // -φ, 0, +1
        Triad<Apfloat> vB8 = new Triad<>(N0, N1, gold_ratio);                                 // 0, +1, +φ
        Triad<Apfloat> vB9 = new Triad<>(N0, N1.multiply(NEG_N1), gold_ratio);                   // 0, -1, +φ
        Triad<Apfloat> vB10 = new Triad<>(N0, N1.multiply(NEG_N1), gold_ratio.multiply(NEG_N1));    // 0, -1, -φ
        Triad<Apfloat> vB11 = new Triad<>(N0, N1, gold_ratio.multiply(NEG_N1));                  // 0, +1, -φ


        // dist from origin = sqrt(x²+y²+z²) = sqrt(φ²+1²) = sqrt(φ²+1)
        Apfloat dist = ApfloatMath.sqrt(N1.add(ApfloatMath.pow(gold_ratio, 2)));
        Apfloat r = super.getRadius();
        String r_dist = r.divide(dist).toString();

        // cartesian vertex = (r_dist * vB_x, r_dist * vB_y, r_dist * vB_z)
        Triad<Apfloat> vC0 = VectorMath.mult(vB0, r_dist);
        Triad<Apfloat> vC1 = VectorMath.mult(vB1, r_dist);
        Triad<Apfloat> vC2 = VectorMath.mult(vB2, r_dist);
        Triad<Apfloat> vC3 = VectorMath.mult(vB3, r_dist);
        Triad<Apfloat> vC4 = VectorMath.mult(vB4, r_dist);
        Triad<Apfloat> vC5 = VectorMath.mult(vB5, r_dist);
        Triad<Apfloat> vC6 = VectorMath.mult(vB6, r_dist);
        Triad<Apfloat> vC7 = VectorMath.mult(vB7, r_dist);
        Triad<Apfloat> vC8 = VectorMath.mult(vB8, r_dist);
        Triad<Apfloat> vC9 = VectorMath.mult(vB9, r_dist);
        Triad<Apfloat> vC10 = VectorMath.mult(vB10, r_dist);
        Triad<Apfloat> vC11 = VectorMath.mult(vB11, r_dist);

        faces = new Icosad<> (
                // Upper cap faces (around vertex vC8)
                new Triad<>(vC8, vC0, vC4),
                new Triad<>(vC8, vC4, vC9),
                new Triad<>(vC8, vC9, vC7),
                new Triad<>(vC8, vC7, vC1),
                new Triad<>(vC8, vC1, vC0),

                // Upper belt faces
                new Triad<>(vC0, vC1, vC11),
                new Triad<>(vC1, vC7, vC6),
                new Triad<>(vC7, vC9, vC2),
                new Triad<>(vC9, vC4, vC3),
                new Triad<>(vC4, vC0, vC5),

                // Lower belt faces
                new Triad<>(vC11, vC1, vC6),
                new Triad<>(vC6, vC7, vC2),
                new Triad<>(vC2, vC9, vC3),
                new Triad<>(vC3, vC4, vC5),
                new Triad<>(vC5, vC0, vC11),

                // Lower cap faces (around vertex vC10)
                new Triad<>(vC10, vC11, vC6),
                new Triad<>(vC10, vC6, vC2),
                new Triad<>(vC10, vC2, vC3),
                new Triad<>(vC10, vC3, vC5),
                new Triad<>(vC10, vC5, vC11)
        );


        // compute face normals
        ArrayList<Triad<Apfloat>> norms = new ArrayList<>();
        for (int i = 0; i < faces.fetchSize(); i++) {

            // get points of face
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> vertB = (Triad<Apfloat>) face.fetch(1);
            Triad<Apfloat> vertC = (Triad<Apfloat>) face.fetch(2);

            // compute u edge vector
            // u = B - A = (vertB_x - vertA_x, vertB_y - vertA_y, vertB_z - vertA_z)
            Triad<Apfloat> u = VectorMath.subs(vertB, vertA);

            // compute v edge vector
            // v = C - A = (vertC_x - vertA_x, vertC_y - vertA_y, vertC_z - vertA_z)
            Triad<Apfloat> v = VectorMath.subs(vertC, vertA);

            // calculate face vector
            // n = u×v = (u_y*v_z - u_z*v_y, u_z*v_x - u_x * v_z, u_x*v_y - u_y*v_x)
            Triad<Apfloat> n = VectorMath.cross_prod(u, v);

            // normalize face vector
            // n_norm = (n_x / ||n||, n_y / ||n||, n_z / ||n||)
            Triad<Apfloat> n_norm = VectorMath.normalize(n);

            // compute centroid vector c (midpoint of triangle)
            Apfloat c_x =   vertA.fetch(0)
                    .add(vertB.fetch(0))
                    .add(vertC.fetch(0))
                    .divide(N3);
            Apfloat c_y =   vertA.fetch(1)
                    .add(vertB.fetch(1))
                    .add(vertC.fetch(1))
                    .divide(N3);
            Apfloat c_z =   vertA.fetch(2)
                    .add(vertB.fetch(2))
                    .add(vertC.fetch(2))
                    .divide(N3);
            Triad<Apfloat> c = new Triad<>(c_x, c_y, c_z);

            // compute dot product of c and n_norm
            // w = n_norm ⋅ centroid vector
            Apfloat w = VectorMath.dot_prod(n_norm, c);

            // w < 0 =>  n_norm points inwards;
            // if above is true, flip such that  n_norm points outwards (negate  n_norm)
            if (w.compareTo(N0) < 0) {
                norms.add(VectorMath.mult(n_norm, "-1"));
            }
            else {
                norms.add(n_norm);
            }
        }

        // post-cond: 20 face norms
        assert(norms.size() == 20);

        // assign to face_norms
        face_norms = new Icosad<>(
                norms.get(0),
                norms.get(1),
                norms.get(2),
                norms.get(3),
                norms.get(4),
                norms.get(5),
                norms.get(6),
                norms.get(7),
                norms.get(8),
                norms.get(9),
                norms.get(10),
                norms.get(11),
                norms.get(12),
                norms.get(13),
                norms.get(14),
                norms.get(15),
                norms.get(16),
                norms.get(17),
                norms.get(18),
                norms.get(19)
        );
    }


    /**
     * Determines whether the given Cartesian coordinates are inside the icosahedral boundary.
     *
     * @param point_cart the coordinates of a point in Cartesian space, must not be null
     * @return {@code true} if the point is within bounds, {@code false} otherwise
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        // check if each point lies within bounds of each face
        for (int i = 0; i < faces.fetchSize(); i++) {

            // get vertices of vertA
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x - vertA_z)
            Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> face_norm = (Triad<Apfloat>) face_norms.fetch(i);
            Apfloat d = VectorMath.dot_prod(face_norm, m);

            // if d < 0,  point lies behind face (within bounds)
            // if d == 0, point lies on face (within bounds)
            // if d > 0,  point lies in front of face (out of bounds)
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
