package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.archimedean.Icosidodecahedron;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;
import java.util.ArrayList;

/**
 * Represents a <b>Rhombic Triacontahedron</b>.
 * <p> A convex polyhedron composed of 30 congruent rhombic faces,
 * 60 edges, and 32 vertices. All vertices lie on a sphere, making it a <i>zonohedron</i>
 * and the dual of the {@link Icosidodecahedron}.
 * <p> The rhombus face angles are 2 × arctan(φ²) and 2 × arctan(1/φ²),
 * where φ is the golden ratio.
 * @see <a href="https://dmccooey.com/polyhedra/RhombicTriacontahedron.html">
 *      Rhombic Triacontahedron (David McCooey)</a>
 */
public class TriacontahedronRhombic extends Shape {

    // 30 rhomboidal faces
    private final ArrayList<Tetrad<Tuple<Apfloat>>> faces_rho = new ArrayList<>();
    private final ArrayList<Triad<Apfloat>> face_norms_rho = new ArrayList<>();

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N8 = new Apfloat("8", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(N5);

    // C0 = sqrt(5) / 4
    final Apfloat C0 = SQRT5.divide(N4);

    // NEG_C0 = -(sqrt(5) / 4)
    final Apfloat NEG_C0 = C0.multiply(NEG_N1);

    // C1 = (5 + sqrt(5)) / 8
    final Apfloat C1 = (N5.add(SQRT5)).divide(N8);

    // NEG_C1 = -((5 + sqrt(5)) / 8)
    final Apfloat NEG_C1 = C1.multiply(NEG_N1);

    // C2 = (5 + 3 * sqrt(5)) / 8
    final Apfloat C2 = (N5.add(N3.multiply(SQRT5))).divide(N8);

    // NEG_C2 = -((5 + 3 * sqrt(5)) / 8)
    final Apfloat NEG_C2 = C2.multiply(NEG_N1);

    public TriacontahedronRhombic(
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
        Triad<Apfloat> vB0 = new Triad<>(C1, N0, C2);
        Triad<Apfloat> vB1 = new Triad<>(C1, N0, NEG_C2);
        Triad<Apfloat> vB2 = new Triad<>(NEG_C1, N0, C2);
        Triad<Apfloat> vB3 = new Triad<>(NEG_C1, N0, NEG_C2);
        Triad<Apfloat> vB4 = new Triad<>(C2, C1, N0);
        Triad<Apfloat> vB5 = new Triad<>(C2, NEG_C1, N0);
        Triad<Apfloat> vB6 = new Triad<>(NEG_C2, C1, N0);
        Triad<Apfloat> vB7 = new Triad<>(NEG_C2, NEG_C1, N0);
        Triad<Apfloat> vB8 = new Triad<>(N0, C2, C1);
        Triad<Apfloat> vB9 = new Triad<>(N0, C2, NEG_C1);
        Triad<Apfloat> vB10 = new Triad<>(N0, NEG_C2, C1);
        Triad<Apfloat> vB11 = new Triad<>(N0, NEG_C2, NEG_C1);
        Triad<Apfloat> vB12 = new Triad<>(N0, C0, C2);
        Triad<Apfloat> vB13 = new Triad<>(N0, C0, NEG_C2);
        Triad<Apfloat> vB14 = new Triad<>(N0, NEG_C0, C2);
        Triad<Apfloat> vB15 = new Triad<>(N0, NEG_C0, NEG_C2);
        Triad<Apfloat> vB16 = new Triad<>(C2, N0, C0);
        Triad<Apfloat> vB17 = new Triad<>(C2, N0, NEG_C0);
        Triad<Apfloat> vB18 = new Triad<>(NEG_C2, N0, C0);
        Triad<Apfloat> vB19 = new Triad<>(NEG_C2, N0, NEG_C0);
        Triad<Apfloat> vB20 = new Triad<>(C0, C2, N0);
        Triad<Apfloat> vB21 = new Triad<>(C0, NEG_C2, N0);
        Triad<Apfloat> vB22 = new Triad<>(NEG_C0, C2, N0);
        Triad<Apfloat> vB23 = new Triad<>(NEG_C0, NEG_C2, N0);
        Triad<Apfloat> vB24 = new Triad<>(C1, C1, C1);
        Triad<Apfloat> vB25 = new Triad<>(C1, C1, NEG_C1);
        Triad<Apfloat> vB26 = new Triad<>(C1, NEG_C1, C1);
        Triad<Apfloat> vB27 = new Triad<>(C1, NEG_C1, NEG_C1);
        Triad<Apfloat> vB28 = new Triad<>(NEG_C1, C1, C1);
        Triad<Apfloat> vB29 = new Triad<>(NEG_C1, C1, NEG_C1);
        Triad<Apfloat> vB30 = new Triad<>(NEG_C1, NEG_C1, C1);
        Triad<Apfloat> vB31 = new Triad<>(NEG_C1, NEG_C1, NEG_C1);

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = mult(normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1  = mult(normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2  = mult(normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3  = mult(normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4  = mult(normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5  = mult(normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6  = mult(normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7  = mult(normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8  = mult(normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9  = mult(normalize(vB9), super.getRadius().toString());
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
        Triad<Apfloat> vC28 = mult(normalize(vB28), super.getRadius().toString());
        Triad<Apfloat> vC29 = mult(normalize(vB29), super.getRadius().toString());
        Triad<Apfloat> vC30 = mult(normalize(vB30), super.getRadius().toString());
        Triad<Apfloat> vC31 = mult(normalize(vB31), super.getRadius().toString());

        // ==== RHOMBOIDAL FACES ====
        Tetrad<Tuple<Apfloat>> rho0 = new Tetrad<>(vC0, vC12, vC2, vC14);
        Triad<Apfloat> rho0_norm = normalQuad(vC0, vC12,vC2,vC14,true);
        faces_rho.add(rho0);
        face_norms_rho.add(rho0_norm);
        Tetrad<Tuple<Apfloat>> rho1 = new Tetrad<>(vC0, vC14, vC10, vC26);
        Triad<Apfloat> rho1_norm = normalQuad(vC0, vC14, vC10, vC26, true);
        faces_rho.add(rho1);
        face_norms_rho.add(rho1_norm);

        Tetrad<Tuple<Apfloat>> rho2 = new Tetrad<>(vC0, vC26, vC5, vC16);
        Triad<Apfloat> rho2_norm = normalQuad(vC0, vC26, vC5, vC16, true);
        faces_rho.add(rho2);
        face_norms_rho.add(rho2_norm);

        Tetrad<Tuple<Apfloat>> rho3 = new Tetrad<>(vC1, vC13, vC9, vC25);
        Triad<Apfloat> rho3_norm = normalQuad(vC1, vC13, vC9, vC25, true);
        faces_rho.add(rho3);
        face_norms_rho.add(rho3_norm);

        Tetrad<Tuple<Apfloat>> rho4 = new Tetrad<>(vC1, vC25, vC4, vC17);
        Triad<Apfloat> rho4_norm = normalQuad(vC1, vC25, vC4, vC17, true);
        faces_rho.add(rho4);
        face_norms_rho.add(rho4_norm);

        Tetrad<Tuple<Apfloat>> rho5 = new Tetrad<>(vC1, vC17, vC5, vC27);
        Triad<Apfloat> rho5_norm = normalQuad(vC1, vC17, vC5, vC27, true);
        faces_rho.add(rho5);
        face_norms_rho.add(rho5_norm);

        Tetrad<Tuple<Apfloat>> rho6 = new Tetrad<>(vC2, vC28, vC6, vC18);
        Triad<Apfloat> rho6_norm = normalQuad(vC2, vC28, vC6, vC18, true);
        faces_rho.add(rho6);
        face_norms_rho.add(rho6_norm);

        Tetrad<Tuple<Apfloat>> rho7 = new Tetrad<>(vC2, vC18, vC7, vC30);
        Triad<Apfloat> rho7_norm = normalQuad(vC2, vC18, vC7, vC30, true);
        faces_rho.add(rho7);
        face_norms_rho.add(rho7_norm);

        Tetrad<Tuple<Apfloat>> rho8 = new Tetrad<>(vC2, vC30, vC10, vC14);
        Triad<Apfloat> rho8_norm = normalQuad(vC2, vC30, vC10, vC14, true);
        faces_rho.add(rho8);
        face_norms_rho.add(rho8_norm);

        Tetrad<Tuple<Apfloat>> rho9 = new Tetrad<>(vC3, vC19, vC6, vC29);
        Triad<Apfloat> rho9_norm = normalQuad(vC3, vC19, vC6, vC29, true);
        faces_rho.add(rho9);
        face_norms_rho.add(rho9_norm);

        Tetrad<Tuple<Apfloat>> rho10 = new Tetrad<>(vC3, vC29, vC9, vC13);
        Triad<Apfloat> rho10_norm = normalQuad(vC3, vC29, vC9, vC13, true);
        faces_rho.add(rho10);
        face_norms_rho.add(rho10_norm);

        Tetrad<Tuple<Apfloat>> rho11 = new Tetrad<>(vC3, vC13, vC1, vC15);
        Triad<Apfloat> rho11_norm = normalQuad(vC3, vC13, vC1, vC15, true);
        faces_rho.add(rho11);
        face_norms_rho.add(rho11_norm);

        Tetrad<Tuple<Apfloat>> rho12 = new Tetrad<>(vC4, vC20, vC8, vC24);
        Triad<Apfloat> rho12_norm = normalQuad(vC4, vC20, vC8, vC24, true);
        faces_rho.add(rho12);
        face_norms_rho.add(rho12_norm);

        Tetrad<Tuple<Apfloat>> rho13 = new Tetrad<>(vC4, vC24, vC0, vC16);
        Triad<Apfloat> rho13_norm = normalQuad(vC4, vC24, vC0, vC16, true);
        faces_rho.add(rho13);
        face_norms_rho.add(rho13_norm);

        Tetrad<Tuple<Apfloat>> rho14 = new Tetrad<>(vC4, vC16, vC5, vC17);
        Triad<Apfloat> rho14_norm = normalQuad(vC4, vC16, vC5, vC17, true);
        faces_rho.add(rho14);
        face_norms_rho.add(rho14_norm);

        Tetrad<Tuple<Apfloat>> rho15 = new Tetrad<>(vC7, vC18, vC6, vC19);
        Triad<Apfloat> rho15_norm = normalQuad(vC7, vC18, vC6, vC19, true);
        faces_rho.add(rho15);
        face_norms_rho.add(rho15_norm);

        Tetrad<Tuple<Apfloat>> rho16 = new Tetrad<>(vC7, vC19, vC3, vC31);
        Triad<Apfloat> rho16_norm = normalQuad(vC7, vC19, vC3, vC31, true);
        faces_rho.add(rho16);
        face_norms_rho.add(rho16_norm);

        Tetrad<Tuple<Apfloat>> rho17 = new Tetrad<>(vC7, vC31, vC11, vC23);
        Triad<Apfloat> rho17_norm = normalQuad(vC7, vC31, vC11, vC23, true);
        faces_rho.add(rho17);
        face_norms_rho.add(rho17_norm);

        Tetrad<Tuple<Apfloat>> rho18 = new Tetrad<>(vC8, vC22, vC6, vC28);
        Triad<Apfloat> rho18_norm = normalQuad(vC8, vC22, vC6, vC28, true);
        faces_rho.add(rho18);
        face_norms_rho.add(rho18_norm);

        Tetrad<Tuple<Apfloat>> rho19 = new Tetrad<>(vC8, vC28, vC2, vC12);
        Triad<Apfloat> rho19_norm = normalQuad(vC8, vC28, vC2, vC12, true);
        faces_rho.add(rho19);
        face_norms_rho.add(rho19_norm);

        Tetrad<Tuple<Apfloat>> rho20 = new Tetrad<>(vC8, vC12, vC0, vC24);
        Triad<Apfloat> rho20_norm = normalQuad(vC8, vC12, vC0, vC24, true);
        faces_rho.add(rho20);
        face_norms_rho.add(rho20_norm);

        Tetrad<Tuple<Apfloat>> rho21 = new Tetrad<>(vC9, vC29, vC6, vC22);
        Triad<Apfloat> rho21_norm = normalQuad(vC9, vC29, vC6, vC22, true);
        faces_rho.add(rho21);
        face_norms_rho.add(rho21_norm);

        Tetrad<Tuple<Apfloat>> rho22 = new Tetrad<>(vC9, vC22, vC8, vC20);
        Triad<Apfloat> rho22_norm = normalQuad(vC9, vC22, vC8, vC20, true);
        faces_rho.add(rho22);
        face_norms_rho.add(rho22_norm);
        Tetrad<Tuple<Apfloat>> rho23 = new Tetrad<>(vC9, vC20, vC4, vC25);
        Triad<Apfloat> rho23_norm = normalQuad(vC9, vC20, vC4, vC25, true);
        faces_rho.add(rho23);
        face_norms_rho.add(rho23_norm);
        Tetrad<Tuple<Apfloat>> rho24 = new Tetrad<>(vC10, vC30, vC7, vC23);
        Triad<Apfloat> rho24_norm = normalQuad(vC10, vC30, vC7, vC23, true);
        faces_rho.add(rho24);
        face_norms_rho.add(rho24_norm);
        Tetrad<Tuple<Apfloat>> rho25 = new Tetrad<>(vC10, vC23, vC11, vC21);
        Triad<Apfloat> rho25_norm = normalQuad(vC10, vC23, vC11, vC21, true);
        faces_rho.add(rho25);
        face_norms_rho.add(rho25_norm);
        Tetrad<Tuple<Apfloat>> rho26 = new Tetrad<>(vC10, vC21, vC5, vC26);
        Triad<Apfloat> rho26_norm = normalQuad(vC10, vC21, vC5, vC26, true);
        faces_rho.add(rho26);
        face_norms_rho.add(rho26_norm);
        Tetrad<Tuple<Apfloat>> rho27 = new Tetrad<>(vC11, vC31, vC3, vC15);
        Triad<Apfloat> rho27_norm = normalQuad(vC11, vC31, vC3, vC15, true);
        faces_rho.add(rho27);
        face_norms_rho.add(rho27_norm);
        Tetrad<Tuple<Apfloat>> rho28 = new Tetrad<>(vC11, vC15, vC1, vC27);
        Triad<Apfloat> rho28_norm = normalQuad(vC11, vC15, vC1, vC27, true);
        faces_rho.add(rho28);
        face_norms_rho.add(rho28_norm);
        Tetrad<Tuple<Apfloat>> rho29 = new Tetrad<>(vC11, vC27, vC5, vC21);
        Triad<Apfloat> rho29_norm = normalQuad(vC11, vC27, vC5, vC21, true);
        faces_rho.add(rho29);
        face_norms_rho.add(rho29_norm);
    }

    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_rho.size(); i++){
            // === Check square faces ===
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_rho.get(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
            Triad<Apfloat> m = subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_rho.get(i);
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }

        return true;
    }
}
