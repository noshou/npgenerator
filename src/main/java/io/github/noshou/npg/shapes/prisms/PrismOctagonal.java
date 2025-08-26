package io.github.noshou.npg.shapes.prisms;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.*;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents an <b>Octagonal Prism</b>.
 * <p> This solid is a prism with two parallel, congruent octagonal bases connected
 * by rectangular faces. It has 10 faces in total — 2 octagonal bases and 8
 * rectangular lateral faces. The shape is a uniform polyhedron with 16 vertices
 * and 24 edges, exhibiting prismatic symmetry along its axis.
 * @see <a href="https://dmccooey.com/polyhedra/OctagonalPrism.html">
 *      Octagonal Prism (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class PrismOctagonal extends Shape {

    // 10 faces: 2 octagons + 8 rectangles
    private final Decad<Tuple<Tuple<Apfloat>>> faces;
    private final Decad<Tuple<Apfloat>> face_norms;

    // constants
    private final Apfloat NEG_N1 = new Apfloat("-1", precision);
    private final Apfloat N1 = new Apfloat("1", precision);
    private final Apfloat N2 = new Apfloat("2", precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N8 = new Apfloat("8", super.precision);


    /**
     * Constructs a new {@code Shapes.PrismOctagonal} instance with the given parameters.
     */
    public PrismOctagonal(
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

        // Constants
        Apfloat SQRT2 = ApfloatMath.sqrt(N2);
        Apfloat r = super.getRadius();

        // For rectangular side faces with all vertices on sphere of radius r:
        // If octagon has radius r_oct and height is h, then r² = r_oct² + h²
        // For rectangular faces: edge_length = height = 2h
        // Octagon edge = 2*r_oct*sin(π/8) = r_oct*√(2-√2)
        // Setting equal: 2h = r_oct*√(2-√2)
        // So: h = r_oct*√(2-√2)/2
        // And: r² = r_oct² + (r_oct*√(2-√2)/2)² = r_oct²(1 + (2-√2)/4)

        Apfloat PI = ApfloatMath.pi(precision);
        Apfloat N4 = new Apfloat("4", precision);

        // Calculate r_oct from the constraint that all vertices lie on sphere
        Apfloat sqrt2_minus_sqrt2 = ApfloatMath.sqrt(N2.subtract(SQRT2));
        Apfloat factor = N1.add(N2.subtract(SQRT2).divide(N4));
        Apfloat r_oct = r.divide(ApfloatMath.sqrt(factor));
        Apfloat h = r_oct.multiply(sqrt2_minus_sqrt2).divide(N2);

        Apfloat[] cos_vals = new Apfloat[8];
        Apfloat[] sin_vals = new Apfloat[8];

        // Calculate octagon vertex positions (0°, 45°, 90°, 135°, 180°, 225°, 270°, 315°)
        for (int i = 0; i < 8; i++) {
            Apfloat angle = new Apfloat(i * 45, precision).multiply(PI).divide(new Apfloat("180", precision));
            cos_vals[i] = ApfloatMath.cos(angle);
            sin_vals[i] = ApfloatMath.sin(angle);
        }

        // Top octagon vertices (z = +h)
        Triad<Apfloat> vTop0 = new Triad<>(r_oct.multiply(cos_vals[0]), r_oct.multiply(sin_vals[0]), h);
        Triad<Apfloat> vTop1 = new Triad<>(r_oct.multiply(cos_vals[1]), r_oct.multiply(sin_vals[1]), h);
        Triad<Apfloat> vTop2 = new Triad<>(r_oct.multiply(cos_vals[2]), r_oct.multiply(sin_vals[2]), h);
        Triad<Apfloat> vTop3 = new Triad<>(r_oct.multiply(cos_vals[3]), r_oct.multiply(sin_vals[3]), h);
        Triad<Apfloat> vTop4 = new Triad<>(r_oct.multiply(cos_vals[4]), r_oct.multiply(sin_vals[4]), h);
        Triad<Apfloat> vTop5 = new Triad<>(r_oct.multiply(cos_vals[5]), r_oct.multiply(sin_vals[5]), h);
        Triad<Apfloat> vTop6 = new Triad<>(r_oct.multiply(cos_vals[6]), r_oct.multiply(sin_vals[6]), h);
        Triad<Apfloat> vTop7 = new Triad<>(r_oct.multiply(cos_vals[7]), r_oct.multiply(sin_vals[7]), h);

        // Bottom octagon vertices (z = -h)
        Triad<Apfloat> vBot0 = new Triad<>(r_oct.multiply(cos_vals[0]), r_oct.multiply(sin_vals[0]), h.multiply(NEG_N1));
        Triad<Apfloat> vBot1 = new Triad<>(r_oct.multiply(cos_vals[1]), r_oct.multiply(sin_vals[1]), h.multiply(NEG_N1));
        Triad<Apfloat> vBot2 = new Triad<>(r_oct.multiply(cos_vals[2]), r_oct.multiply(sin_vals[2]), h.multiply(NEG_N1));
        Triad<Apfloat> vBot3 = new Triad<>(r_oct.multiply(cos_vals[3]), r_oct.multiply(sin_vals[3]), h.multiply(NEG_N1));
        Triad<Apfloat> vBot4 = new Triad<>(r_oct.multiply(cos_vals[4]), r_oct.multiply(sin_vals[4]), h.multiply(NEG_N1));
        Triad<Apfloat> vBot5 = new Triad<>(r_oct.multiply(cos_vals[5]), r_oct.multiply(sin_vals[5]), h.multiply(NEG_N1));
        Triad<Apfloat> vBot6 = new Triad<>(r_oct.multiply(cos_vals[6]), r_oct.multiply(sin_vals[6]), h.multiply(NEG_N1));
        Triad<Apfloat> vBot7 = new Triad<>(r_oct.multiply(cos_vals[7]), r_oct.multiply(sin_vals[7]), h.multiply(NEG_N1));

        // Top octagon face (vertices in counterclockwise order when viewed from above)
        Octad<Tuple<Apfloat>> top_face = new Octad<>(vTop0, vTop1, vTop2, vTop3, vTop4, vTop5, vTop6, vTop7);
        Triad<Apfloat> top_norm = new Triad<>(new Apfloat("0", precision), new Apfloat("0", precision), N1);

        // Bottom octagon face (vertices in clockwise order when viewed from above, so CCW from below)
        Octad<Tuple<Apfloat>> bottom_face = new Octad<>(vBot7, vBot6, vBot5, vBot4, vBot3, vBot2, vBot1, vBot0);
        Triad<Apfloat> bottom_norm = new Triad<>(new Apfloat("0", precision), new Apfloat("0", precision), NEG_N1);

        // 8 rectangular side faces (each connects corresponding top/bottom edges)
        Tetrad<Tuple<Apfloat>> side0 = new Tetrad<>(vTop0, vBot0, vBot1, vTop1);
        Triad<Apfloat> side0_norm = normalQuad(vTop0, vBot0, vBot1, vTop1, true);

        Tetrad<Tuple<Apfloat>> side1 = new Tetrad<>(vTop1, vBot1, vBot2, vTop2);
        Triad<Apfloat> side1_norm = normalQuad(vTop1, vBot1, vBot2, vTop2, true);

        Tetrad<Tuple<Apfloat>> side2 = new Tetrad<>(vTop2, vBot2, vBot3, vTop3);
        Triad<Apfloat> side2_norm = normalQuad(vTop2, vBot2, vBot3, vTop3, true);

        Tetrad<Tuple<Apfloat>> side3 = new Tetrad<>(vTop3, vBot3, vBot4, vTop4);
        Triad<Apfloat> side3_norm = normalQuad(vTop3, vBot3, vBot4, vTop4, true);

        Tetrad<Tuple<Apfloat>> side4 = new Tetrad<>(vTop4, vBot4, vBot5, vTop5);
        Triad<Apfloat> side4_norm = normalQuad(vTop4, vBot4, vBot5, vTop5, true);

        Tetrad<Tuple<Apfloat>> side5 = new Tetrad<>(vTop5, vBot5, vBot6, vTop6);
        Triad<Apfloat> side5_norm = normalQuad(vTop5, vBot5, vBot6, vTop6, true);

        Tetrad<Tuple<Apfloat>> side6 = new Tetrad<>(vTop6, vBot6, vBot7, vTop7);
        Triad<Apfloat> side6_norm = normalQuad(vTop6, vBot6, vBot7, vTop7, true);

        Tetrad<Tuple<Apfloat>> side7 = new Tetrad<>(vTop7, vBot7, vBot0, vTop0);
        Triad<Apfloat> side7_norm = normalQuad(vTop7, vBot7, vBot0, vTop0, true);

        faces = new Decad<>(
                top_face,
                bottom_face,
                side0,
                side1,
                side2,
                side3,
                side4,
                side5,
                side6,
                side7
        );

        face_norms = new Decad<>(
                top_norm,
                bottom_norm,
                side0_norm,
                side1_norm,
                side2_norm,
                side3_norm,
                side4_norm,
                side5_norm,
                side6_norm,
                side7_norm
        );
    }

    /**
     * Determines whether the given Cartesian coordinates are inside the boundary.
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < faces.fetchSize(); i++) {
            Tuple<Tuple<Apfloat>> face = faces.fetch(i);
            Triad<Apfloat> face_norm = (Triad<Apfloat>) face_norms.fetch(i);

            // Calculate face centroid based on number of vertices
            Triad<Apfloat> centroid;
            if (i < 2) {  // Octagonal faces (top and bottom)
                Octad<Tuple<Apfloat>> octFace = (Octad<Tuple<Apfloat>>) face;
                Apfloat sum_x = new Apfloat("0", super.precision);
                Apfloat sum_y = new Apfloat("0", super.precision);
                Apfloat sum_z = new Apfloat("0", super.precision);

                for (int j = 0; j < 8; j++) {
                    Triad<Apfloat> vertex = (Triad<Apfloat>) octFace.fetch(j);
                    sum_x = sum_x.add(vertex.fetch(0));
                    sum_y = sum_y.add(vertex.fetch(1));
                    sum_z = sum_z.add(vertex.fetch(2));
                }
                centroid = new Triad<>(sum_x.divide(N8), sum_y.divide(N8), sum_z.divide(N8));
            } else {  // Rectangular faces
                Tetrad<Tuple<Apfloat>> rectFace = (Tetrad<Tuple<Apfloat>>) face;
                Triad<Apfloat> v0 = (Triad<Apfloat>) rectFace.fetch(0);
                Triad<Apfloat> v1 = (Triad<Apfloat>) rectFace.fetch(1);
                Triad<Apfloat> v2 = (Triad<Apfloat>) rectFace.fetch(2);
                Triad<Apfloat> v3 = (Triad<Apfloat>) rectFace.fetch(3);

                centroid = new Triad<>(
                        v0.fetch(0).add(v1.fetch(0)).add(v2.fetch(0)).add(v3.fetch(0)).divide(N4),
                        v0.fetch(1).add(v1.fetch(1)).add(v2.fetch(1)).add(v3.fetch(1)).divide(N4),
                        v0.fetch(2).add(v1.fetch(2)).add(v2.fetch(2)).add(v3.fetch(2)).divide(N4)
                );
            }

            // Vector from centroid to point
            Triad<Apfloat> m = subs(point_cart, centroid);

            // Dot product with face normal
            Apfloat d = dot_prod(face_norm, m);

            // If point is in front of face (outside), reject
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}