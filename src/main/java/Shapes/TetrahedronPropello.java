package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/PropelloTetrahedron.txt...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class TetrahedronPropello extends Shape {

    // 12 kite faces
    Dodecad<Tuple<Tuple<Apfloat>>> faces_kte;
    Dodecad<Tuple<Apfloat>> face_norms_kte;

    // 4 triangular faces
    Tetrad<Tuple<Tuple<Apfloat>>> faces_tri;
    Tetrad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N11 = new Apfloat("11", super.precision);
    private final Apfloat N25 = new Apfloat("25", super.precision);
    private final Apfloat N33 = new Apfloat("33", super.precision);
    private final Apfloat N69 = new Apfloat("69", super.precision);
    private final Apfloat SQRT69 = ApfloatMath.sqrt(N69);
    private final Apfloat N371 = new Apfloat("371", super.precision); // 371

    // C0 = (cbrt(4 * (11 + 3 * sqrt(69))) - cbrt(4 * (3 * sqrt(69) - 11)) - 1) / 3
    private final Apfloat C0 = (
            ApfloatMath.cbrt(N4.multiply(N11.add(N3.multiply(SQRT69))))
            .subtract(
            ApfloatMath.cbrt(N4.multiply(N3.multiply(SQRT69).subtract(N11)))
            ).subtract(N1))
            .divide(N3);

    // C1 = (cbrt(4 * (25 + 3 * sqrt(69))) + cbrt(4 * (25 - 3 * sqrt(69))) - 5) / 3
    private final Apfloat C1 = (
            ApfloatMath.cbrt(N4.multiply(N25.add(N3.multiply(SQRT69))))
            .add(
            ApfloatMath.cbrt(N4.multiply(N25.subtract(N3.multiply(SQRT69))))
            ).subtract(N5))
            .divide(N3);

    // C2 = (cbrt(4 * (371 + 33*sqrt(69))) + cbrt(4 * (371 - 33*sqrt(69))) - 1) / 33
    private final Apfloat C2 = (
            ApfloatMath.cbrt(N4.multiply(N371.add(N33.multiply(SQRT69))))
            .add(
            ApfloatMath.cbrt(N4.multiply(N371.subtract(N33.multiply(SQRT69))))
            ).subtract(N1))
            .divide(N33);

    public TetrahedronPropello(
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

        // ==== CANONICAL BASIS VERTICES ====
        Triad<Apfloat> vB0 = new Triad<>(C1, C0, N1);
        Triad<Apfloat> vB1 = new Triad<>(C1, C0.multiply(NEG_N1), NEG_N1);
        Triad<Apfloat> vB2 = new Triad<>(C1.multiply(NEG_N1), C0.multiply(NEG_N1), N1);
        Triad<Apfloat> vB3 = new Triad<>(C1.multiply(NEG_N1), C0, NEG_N1);
        Triad<Apfloat> vB4 = new Triad<>(N1, C1, C0);
        Triad<Apfloat> vB5 = new Triad<>(N1, C1.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB6 = new Triad<>(NEG_N1, C1.multiply(NEG_N1), C0);
        Triad<Apfloat> vB7 = new Triad<>(NEG_N1, C1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB8 = new Triad<>(C0, N1, C1);
        Triad<Apfloat> vB9 = new Triad<>(C0, NEG_N1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(C0.multiply(NEG_N1), NEG_N1, C1);
        Triad<Apfloat> vB11 = new Triad<>(C0.multiply(NEG_N1), N1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(C2, C2.multiply(NEG_N1), C2);
        Triad<Apfloat> vB13 = new Triad<>(C2, C2, C2.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(C2.multiply(NEG_N1), C2, C2);
        Triad<Apfloat> vB15 = new Triad<>(C2.multiply(NEG_N1), C2.multiply(NEG_N1), C2.multiply(NEG_N1));

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = VectorMath.mult(VectorMath.normalize(vB0),  super.getRadius().toString());
        Triad<Apfloat> vC1  = VectorMath.mult(VectorMath.normalize(vB1),  super.getRadius().toString());
        Triad<Apfloat> vC2  = VectorMath.mult(VectorMath.normalize(vB2),  super.getRadius().toString());
        Triad<Apfloat> vC3  = VectorMath.mult(VectorMath.normalize(vB3),  super.getRadius().toString());
        Triad<Apfloat> vC4  = VectorMath.mult(VectorMath.normalize(vB4),  super.getRadius().toString());
        Triad<Apfloat> vC5  = VectorMath.mult(VectorMath.normalize(vB5),  super.getRadius().toString());
        Triad<Apfloat> vC6  = VectorMath.mult(VectorMath.normalize(vB6),  super.getRadius().toString());
        Triad<Apfloat> vC7  = VectorMath.mult(VectorMath.normalize(vB7),  super.getRadius().toString());
        Triad<Apfloat> vC8  = VectorMath.mult(VectorMath.normalize(vB8),  super.getRadius().toString());
        Triad<Apfloat> vC9  = VectorMath.mult(VectorMath.normalize(vB9),  super.getRadius().toString());
        Triad<Apfloat> vC10 = VectorMath.mult(VectorMath.normalize(vB10), super.getRadius().toString());
        Triad<Apfloat> vC11 = VectorMath.mult(VectorMath.normalize(vB11), super.getRadius().toString());
        Triad<Apfloat> vC12 = VectorMath.mult(VectorMath.normalize(vB12), super.getRadius().toString());
        Triad<Apfloat> vC13 = VectorMath.mult(VectorMath.normalize(vB13), super.getRadius().toString());
        Triad<Apfloat> vC14 = VectorMath.mult(VectorMath.normalize(vB14), super.getRadius().toString());
        Triad<Apfloat> vC15 = VectorMath.mult(VectorMath.normalize(vB15), super.getRadius().toString());

        // ==== KITE FACES ====
        Tetrad<Tuple<Apfloat>> kte0 = new Tetrad<>(vC12, vC0, vC2, vC10);
        Triad<Apfloat> kte0_norm = VectorMath.normalQuad(vC12, vC0, vC2, vC10, true);
        Tetrad<Tuple<Apfloat>> kte1 = new Tetrad<>(vC12, vC10, vC9, vC5);
        Triad<Apfloat> kte1_norm = VectorMath.normalQuad(vC12, vC10, vC9, vC5, true);
        Tetrad<Tuple<Apfloat>> kte2 = new Tetrad<>(vC12, vC5, vC4, vC0);
        Triad<Apfloat> kte2_norm = VectorMath.normalQuad(vC12, vC5, vC4, vC0, true);
        Tetrad<Tuple<Apfloat>> kte3 = new Tetrad<>(vC13, vC1, vC3, vC11);
        Triad<Apfloat> kte3_norm = VectorMath.normalQuad(vC13, vC1, vC3, vC11, true);
        Tetrad<Tuple<Apfloat>> kte4 = new Tetrad<>(vC13, vC11, vC8, vC4);
        Triad<Apfloat> kte4_norm = VectorMath.normalQuad(vC13, vC11, vC8, vC4, true);
        Tetrad<Tuple<Apfloat>> kte5 = new Tetrad<>(vC13, vC4, vC5, vC1);
        Triad<Apfloat> kte5_norm = VectorMath.normalQuad(vC13, vC4, vC5, vC1, true);
        Tetrad<Tuple<Apfloat>> kte6 = new Tetrad<>(vC14, vC2, vC0, vC8);
        Triad<Apfloat> kte6_norm = VectorMath.normalQuad(vC14, vC2, vC0, vC8, true);
        Tetrad<Tuple<Apfloat>> kte7 = new Tetrad<>(vC14, vC8, vC11, vC7);
        Triad<Apfloat> kte7_norm = VectorMath.normalQuad(vC14, vC8, vC11, vC7, true);
        Tetrad<Tuple<Apfloat>> kte8 = new Tetrad<>(vC14, vC7, vC6, vC2);
        Triad<Apfloat> kte8_norm = VectorMath.normalQuad(vC14, vC7, vC6, vC2, true);
        Tetrad<Tuple<Apfloat>> kte9 = new Tetrad<>(vC15, vC3, vC1, vC9);
        Triad<Apfloat> kte9_norm = VectorMath.normalQuad(vC15, vC3, vC1, vC9, true);
        Tetrad<Tuple<Apfloat>> kte10 = new Tetrad<>(vC15, vC9, vC10, vC6);
        Triad<Apfloat> kte10_norm = VectorMath.normalQuad(vC15, vC9, vC10, vC6, true);
        Tetrad<Tuple<Apfloat>> kte11 = new Tetrad<>(vC15, vC6, vC7, vC3);
        Triad<Apfloat> kte11_norm = VectorMath.normalQuad(vC15, vC6, vC7, vC3, true);
        faces_kte = new Dodecad<>(
                kte0, kte1, kte2, kte3,
                kte4, kte5, kte6, kte7,
                kte8, kte9, kte10, kte11
        );
        face_norms_kte = new Dodecad<>(
                kte0_norm, kte1_norm, kte2_norm, kte3_norm,
                kte4_norm, kte5_norm, kte6_norm, kte7_norm,
                kte8_norm, kte9_norm, kte10_norm, kte11_norm
        );

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC0, vC4, vC8);
        Triad<Apfloat> tri0_norm = VectorMath.normalTriple(vC0, vC4, vC8, true);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC1, vC5, vC9);
        Triad<Apfloat> tri1_norm = VectorMath.normalTriple(vC1, vC5, vC9, true);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC2, vC6, vC10);
        Triad<Apfloat> tri2_norm = VectorMath.normalTriple(vC2, vC6, vC10, true);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC3, vC7, vC11);
        Triad<Apfloat> tri3_norm = VectorMath.normalTriple(vC3, vC7, vC11, true);
        faces_tri = new Tetrad<>(tri0, tri1, tri2, tri3);
        face_norms_tri = new Tetrad<>(tri0_norm, tri1_norm, tri2_norm, tri3_norm);
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < 12; i++) {

            // triangular faces
            if (i < 4) {
                Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.fetch(i);

                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Apfloat d = VectorMath.dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // square faces
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_kte.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_kte.fetch(i);

            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
            Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Apfloat d = VectorMath.dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
