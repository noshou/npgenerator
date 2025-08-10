package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b>6-Truncated Triakis Tetrahedron</b>.
 * <p>This polyhedron is a modified form of the TetrahedronTriakis}
 * in which each face is truncated in a way that produces additional hexagonal and
 * triangular faces.
 * <p>The name "6-truncated" refers to the degree of truncation,
 * which alters the edge lengths while preserving overall symmetry.
 * It retains tetrahedral symmetry and can be classified as a convex,
 * vertex-transitive polyhedron with a combination of regular and irregular faces.
 * @see <a href="https://dmccooey.com/polyhedra/6TruncatedTriakisTetrahedron2.html">
 *      6-Truncated Triakis Tetrahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class TetrahedronTriakisTruncated6 extends Shape {

    // 12 pentagonal faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_pnt;
    private final Dodecad<Tuple<Apfloat>> face_norms_pnt;

    // 4 hexagonal faces
    private final Tetrad<Tuple<Tuple<Apfloat>>> faces_hex;
    private final Tetrad<Tuple<Apfloat>> face_norms_hex;

    // Apfloat Constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N12 = new Apfloat("12", super.precision);
    private final Apfloat N20 = new Apfloat("20", super.precision);
    private final Apfloat SQRT2 = ApfloatMath.sqrt(N2);
    private final Apfloat SQRT6 = ApfloatMath.sqrt(new Apfloat("6", super.precision));

    // C0 = (2 * sqrt(6) - 3 * sqrt(2)) / 12
    private final Apfloat C0 = (N2.multiply(SQRT6).subtract(N3.multiply(SQRT2))).divide(N12);

    // C1 = sqrt(2) / 4
    private final Apfloat C1 = SQRT2.divide(N4);

    // C2 = 3 * (sqrt(2) + 2 * sqrt(6)) / 20
    private final Apfloat C2 = N3.multiply(SQRT2.add(N2.multiply(SQRT6))).divide(N20);

    // C3 = (3 * sqrt(2) + 4 * sqrt(6)) / 12
    private final Apfloat C3 = (N3.multiply(SQRT2).add(N4.multiply(SQRT6))).divide(N12);

    // C4 = (sqrt(2) + 2 * sqrt(6)) / 4
    private final Apfloat C4 = (SQRT2.add(N2.multiply(SQRT6))).divide(N4);

    public TetrahedronTriakisTruncated6(
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

        // ==== BASIS VERTICES ====
        Triad<Apfloat> vB0 = new Triad<>( C1,  C1,  C4);
        Triad<Apfloat> vB1 = new Triad<>( C1, NEG_N1.multiply(C1), NEG_N1.multiply(C4));
        Triad<Apfloat> vB2 = new Triad<>(NEG_N1.multiply(C1), NEG_N1.multiply(C1),  C4);
        Triad<Apfloat> vB3 = new Triad<>(NEG_N1.multiply(C1),  C1, NEG_N1.multiply(C4));
        Triad<Apfloat> vB4 = new Triad<>( C4,  C1,  C1);
        Triad<Apfloat> vB5 = new Triad<>( C4, NEG_N1.multiply(C1), NEG_N1.multiply(C1));
        Triad<Apfloat> vB6 = new Triad<>(NEG_N1.multiply(C4), NEG_N1.multiply(C1),  C1);
        Triad<Apfloat> vB7 = new Triad<>(NEG_N1.multiply(C4),  C1, NEG_N1.multiply(C1));
        Triad<Apfloat> vB8 = new Triad<>( C1,  C4,  C1);
        Triad<Apfloat> vB9 = new Triad<>( C1, NEG_N1.multiply(C4), NEG_N1.multiply(C1));
        Triad<Apfloat> vB10 = new Triad<>(NEG_N1.multiply(C1), NEG_N1.multiply(C4),  C1);
        Triad<Apfloat> vB11 = new Triad<>(NEG_N1.multiply(C1),  C4, NEG_N1.multiply(C1));
        Triad<Apfloat> vB12 = new Triad<>( C3, NEG_N1.multiply(C0),  C3);
        Triad<Apfloat> vB13 = new Triad<>( C3,  C0, NEG_N1.multiply(C3));
        Triad<Apfloat> vB14 = new Triad<>(NEG_N1.multiply(C3),  C0,  C3);
        Triad<Apfloat> vB15 = new Triad<>(NEG_N1.multiply(C3), NEG_N1.multiply(C0), NEG_N1.multiply(C3));
        Triad<Apfloat> vB16 = new Triad<>( C3, NEG_N1.multiply(C3),  C0);
        Triad<Apfloat> vB17 = new Triad<>( C3,  C3, NEG_N1.multiply(C0));
        Triad<Apfloat> vB18 = new Triad<>(NEG_N1.multiply(C3),  C3,  C0);
        Triad<Apfloat> vB19 = new Triad<>(NEG_N1.multiply(C3), NEG_N1.multiply(C3), NEG_N1.multiply(C0));
        Triad<Apfloat> vB20 = new Triad<>( C0, NEG_N1.multiply(C3),  C3);
        Triad<Apfloat> vB21 = new Triad<>( C0,  C3, NEG_N1.multiply(C3));
        Triad<Apfloat> vB22 = new Triad<>(NEG_N1.multiply(C0),  C3,  C3);
        Triad<Apfloat> vB23 = new Triad<>(NEG_N1.multiply(C0), NEG_N1.multiply(C3), NEG_N1.multiply(C3));
        Triad<Apfloat> vB24 = new Triad<>( C2, NEG_N1.multiply(C2),  C2);
        Triad<Apfloat> vB25 = new Triad<>( C2,  C2, NEG_N1.multiply(C2));
        Triad<Apfloat> vB26 = new Triad<>(NEG_N1.multiply(C2),  C2,  C2);
        Triad<Apfloat> vB27 = new Triad<>(NEG_N1.multiply(C2), NEG_N1.multiply(C2), NEG_N1.multiply(C2));

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = mult(normalize(vB0),  super.getRadius().toString());
        Triad<Apfloat> vC1  = mult(normalize(vB1),  super.getRadius().toString());
        Triad<Apfloat> vC2  = mult(normalize(vB2),  super.getRadius().toString());
        Triad<Apfloat> vC3  = mult(normalize(vB3),  super.getRadius().toString());
        Triad<Apfloat> vC4  = mult(normalize(vB4),  super.getRadius().toString());
        Triad<Apfloat> vC5  = mult(normalize(vB5),  super.getRadius().toString());
        Triad<Apfloat> vC6  = mult(normalize(vB6),  super.getRadius().toString());
        Triad<Apfloat> vC7  = mult(normalize(vB7),  super.getRadius().toString());
        Triad<Apfloat> vC8  = mult(normalize(vB8),  super.getRadius().toString());
        Triad<Apfloat> vC9  = mult(normalize(vB9),  super.getRadius().toString());
        Triad<Apfloat> vC10 = mult(normalize(vB10), super.getRadius().toString());
        Triad<Apfloat> vC11 = mult(normalize(vB11), super.getRadius().toString());
        Triad<Apfloat> vC12 = mult(normalize(vB12), super.getRadius().toString());
        Triad<Apfloat> vC13 = mult(normalize(vB13), super.getRadius().toString());
        Triad<Apfloat> vC14 = mult(normalize(vB14), super.getRadius().toString());
        Triad<Apfloat> vC15 = mult(normalize(vB15), super.getRadius().toString());
        Triad<Apfloat> vC16 = mult(normalize(vB16), super.getRadius().toString());
        Triad<Apfloat> vC17 = mult(normalize(vB17), super.getRadius().toString());
        Triad<Apfloat> vC18 = mult(normalize(vB18), super.getRadius().toString());
        Triad<Apfloat> vC19 = mult(normalize(vB19), super.getRadius().toString());
        Triad<Apfloat> vC20 = mult(normalize(vB20), super.getRadius().toString());
        Triad<Apfloat> vC21 = mult(normalize(vB21), super.getRadius().toString());
        Triad<Apfloat> vC22 = mult(normalize(vB22), super.getRadius().toString());
        Triad<Apfloat> vC23 = mult(normalize(vB23), super.getRadius().toString());
        Triad<Apfloat> vC24 = mult(normalize(vB24), super.getRadius().toString());
        Triad<Apfloat> vC25 = mult(normalize(vB25), super.getRadius().toString());
        Triad<Apfloat> vC26 = mult(normalize(vB26), super.getRadius().toString());
        Triad<Apfloat> vC27 = mult(normalize(vB27), super.getRadius().toString());

        // ==== PENTAGONAL FACES ====
        Pentad<Tuple<Apfloat>> pnt0 = new Pentad<>(vC24, vC12, vC0, vC2, vC20);
        Triad<Apfloat> pnt0_norm = normalPent(vC24, vC12, vC0, vC2, vC20, true);
        Pentad<Tuple<Apfloat>> pnt1 = new Pentad<>(vC24, vC20, vC10, vC9, vC16);
        Triad<Apfloat> pnt1_norm = normalPent(vC24, vC20, vC10, vC9, vC16, true);
        Pentad<Tuple<Apfloat>> pnt2 = new Pentad<>(vC24, vC16, vC5, vC4, vC12);
        Triad<Apfloat> pnt2_norm = normalPent(vC24, vC16, vC5, vC4, vC12, true);
        Pentad<Tuple<Apfloat>> pnt3 = new Pentad<>(vC25, vC13, vC1, vC3, vC21);
        Triad<Apfloat> pnt3_norm = normalPent(vC25, vC13, vC1, vC3, vC21, true);
        Pentad<Tuple<Apfloat>> pnt4 = new Pentad<>(vC25, vC21, vC11, vC8, vC17);
        Triad<Apfloat> pnt4_norm = normalPent(vC25, vC21, vC11, vC8, vC17, true);
        Pentad<Tuple<Apfloat>> pnt5 = new Pentad<>(vC25, vC17, vC4, vC5, vC13);
        Triad<Apfloat> pnt5_norm = normalPent(vC25, vC17, vC4, vC5, vC13, true);
        Pentad<Tuple<Apfloat>> pnt6 = new Pentad<>(vC26, vC14, vC2, vC0, vC22);
        Triad<Apfloat> pnt6_norm = normalPent(vC26, vC14, vC2, vC0, vC22, true);
        Pentad<Tuple<Apfloat>> pnt7 = new Pentad<>(vC26, vC22, vC8, vC11, vC18);
        Triad<Apfloat> pnt7_norm = normalPent(vC26, vC22, vC8, vC11, vC18, true);
        Pentad<Tuple<Apfloat>> pnt8 = new Pentad<>(vC26, vC18, vC7, vC6, vC14);
        Triad<Apfloat> pnt8_norm = normalPent(vC26, vC18, vC7, vC6, vC14, true);
        Pentad<Tuple<Apfloat>> pnt9 = new Pentad<>(vC27, vC15, vC3, vC1, vC23);
        Triad<Apfloat> pnt9_norm = normalPent(vC27, vC15, vC3, vC1, vC23, true);
        Pentad<Tuple<Apfloat>> pnt10 = new Pentad<>(vC27, vC23, vC9, vC10, vC19);
        Triad<Apfloat> pnt10_norm = normalPent(vC27, vC23, vC9, vC10, vC19, true);
        Pentad<Tuple<Apfloat>> pnt11 = new Pentad<>(vC27, vC19, vC6, vC7, vC15);
        Triad<Apfloat> pnt11_norm = normalPent(vC27, vC19, vC6, vC7, vC15, true);
        faces_pnt = new Dodecad<>(
                pnt0, pnt1, pnt2, pnt3,
                pnt4, pnt5, pnt6, pnt7,
                pnt8, pnt9, pnt10, pnt11
        );
        face_norms_pnt = new Dodecad<>(
                pnt0_norm, pnt1_norm, pnt2_norm, pnt3_norm,
                pnt4_norm, pnt5_norm, pnt6_norm, pnt7_norm,
                pnt8_norm, pnt9_norm, pnt10_norm, pnt11_norm
        );

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC0, vC12, vC4, vC17, vC8, vC22);
        Triad<Apfloat> hex0_norm = normalHex(vC0, vC12, vC4, vC17, vC8, vC22, true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC1, vC13, vC5, vC16, vC9, vC23);
        Triad<Apfloat> hex1_norm = normalHex(vC1, vC13, vC5, vC16, vC9, vC23, true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC2, vC14, vC6, vC19, vC10, vC20);
        Triad<Apfloat> hex2_norm = normalHex(vC2, vC14, vC6, vC19, vC10, vC20, true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC3, vC15, vC7, vC18, vC11, vC21);
        Triad<Apfloat> hex3_norm = normalHex(vC3, vC15, vC7, vC18, vC11, vC21, true);
        faces_hex = new Tetrad<>(hex0, hex1, hex2, hex3);
        face_norms_hex = new Tetrad<>(hex0_norm, hex1_norm, hex2_norm, hex3_norm);

    }
    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 12; i++) {
            // hexagonal faces
            if (i < 4) {
                Hexad<Tuple<Apfloat>> face = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_hex.fetch(i);

                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Apfloat d = dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // pentagonal faces  faces
            Pentad<Tuple<Apfloat>> face = (Pentad<Tuple<Apfloat>>) faces_pnt.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_pnt.fetch(i);

            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
            Triad<Apfloat> m = subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
