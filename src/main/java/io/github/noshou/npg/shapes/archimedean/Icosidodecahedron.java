package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.platonic.Dodecahedron;
import io.github.noshou.npg.shapes.platonic.Icosahedron;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents an <b>Icosidodecahedron</b>.
 * <p> This Archimedean solid features 32 faces:
 * 20 equilateral triangles and 12 regular pentagons.
 * <p> It is formed as the convex hull of the vertices shared by an
 * {@link Icosahedron} and a {@link Dodecahedron}.
 * @see <a href="https://dmccooey.com/polyhedra/Icosidodecahedron.html">
 *      Icosidodecahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class Icosidodecahedron extends Shape {

    // 20 triangular faces
    private final Icosad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Icosad<Tuple<Apfloat>> face_norms_tri;

    // 12 pentagonal faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_pnt;
    private final Dodecad<Tuple<Apfloat>> face_norms_pnt;


    // Apfloat Constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(new Apfloat("5", super.precision));
    private final Apfloat HALF = new Apfloat("0.5", super.precision);

    // C0 = (1 + sqrt(5)) / 4
    private final Apfloat C0 = (N1.add(SQRT5)).divide(N4);

    // C1 = (3 + sqrt(5)) / 4
    private final Apfloat C1 = (N3.add(SQRT5)).divide(N4);

      // C2 = (1 + sqrt(5)) / 2
    private final Apfloat C2 = (N1.add(SQRT5)).divide(N2);

    public Icosidodecahedron(
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
        Triad<Apfloat> vB0 = new Triad<>(N0, N0, C2);
        Triad<Apfloat> vB1 = new Triad<>(N0, N0, C2.multiply(NEG_N1));
        Triad<Apfloat> vB2 = new Triad<>(C2, N0, N0);
        Triad<Apfloat> vB3 = new Triad<>(C2.multiply(NEG_N1), N0, N0);
        Triad<Apfloat> vB4 = new Triad<>(N0, C2, N0);
        Triad<Apfloat> vB5 = new Triad<>(N0, C2.multiply(NEG_N1), N0);
        Triad<Apfloat> vB6 = new Triad<>(HALF, C0, C1);
        Triad<Apfloat> vB7 = new Triad<>(HALF, C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB8 = new Triad<>(HALF, C0.multiply(NEG_N1), C1);
        Triad<Apfloat> vB9 = new Triad<>(HALF, C0.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(HALF.multiply(NEG_N1), C0, C1);
        Triad<Apfloat> vB11 = new Triad<>(HALF.multiply(NEG_N1), C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(HALF.multiply(NEG_N1), C0.multiply(NEG_N1), C1);
        Triad<Apfloat> vB13 = new Triad<>(HALF.multiply(NEG_N1), C0.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>( C1, HALF, C0);
        Triad<Apfloat> vB15 = new Triad<>( C1, HALF, C0.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>( C1, HALF.multiply(NEG_N1), C0);
        Triad<Apfloat> vB17 = new Triad<>( C1, HALF.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>( C1.multiply(NEG_N1), HALF, C0);
        Triad<Apfloat> vB19 = new Triad<>( C1.multiply(NEG_N1), HALF, C0.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>( C1.multiply(NEG_N1), HALF.multiply(NEG_N1), C0);
        Triad<Apfloat> vB21 = new Triad<>( C1.multiply(NEG_N1), HALF.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB22 = new Triad<>( C0, C1, HALF);
        Triad<Apfloat> vB23 = new Triad<>( C0, C1, HALF.multiply(NEG_N1));
        Triad<Apfloat> vB24 = new Triad<>( C0, C1.multiply(NEG_N1), HALF);
        Triad<Apfloat> vB25 = new Triad<>( C0, C1.multiply(NEG_N1), HALF.multiply(NEG_N1));
        Triad<Apfloat> vB26 = new Triad<>( C0.multiply(NEG_N1), C1, HALF);
        Triad<Apfloat> vB27 = new Triad<>( C0.multiply(NEG_N1), C1, HALF.multiply(NEG_N1));
        Triad<Apfloat> vB28 = new Triad<>( C0.multiply(NEG_N1), C1.multiply(NEG_N1), HALF);
        Triad<Apfloat> vB29 = new Triad<>( C0.multiply(NEG_N1), C1.multiply(NEG_N1), HALF.multiply(NEG_N1));


        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0 = mult(normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1 = mult(normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2 = mult(normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3 = mult(normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4 = mult(normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5 = mult(normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6 = mult(normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7 = mult(normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8 = mult(normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9 = mult(normalize(vB9), super.getRadius().toString());
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

        // ==== Triangular FACES ====
        Triad<Tuple<Apfloat>> tri0  = new Triad<>(vC0,  vC12, vC8);
        Triad<Apfloat>        tri0_norm  = normalTriple(vC0,  vC12, vC8,  true);
        Triad<Tuple<Apfloat>> tri1  = new Triad<>(vC1,  vC9,  vC13);
        Triad<Apfloat>        tri1_norm  = normalTriple(vC1,  vC9,  vC13, true);
        Triad<Tuple<Apfloat>> tri2  = new Triad<>(vC1,  vC11, vC7);
        Triad<Apfloat>        tri2_norm  = normalTriple(vC1,  vC11, vC7,  true);
        Triad<Tuple<Apfloat>> tri3  = new Triad<>(vC2,  vC14, vC16);
        Triad<Apfloat>        tri3_norm  = normalTriple(vC2,  vC14, vC16, true);
        Triad<Tuple<Apfloat>> tri4  = new Triad<>(vC2,  vC17, vC15);
        Triad<Apfloat>        tri4_norm  = normalTriple(vC2,  vC17, vC15, true);
        Triad<Tuple<Apfloat>> tri5  = new Triad<>(vC3,  vC19, vC21);
        Triad<Apfloat>        tri5_norm  = normalTriple(vC3,  vC19, vC21, true);
        Triad<Tuple<Apfloat>> tri6  = new Triad<>(vC3,  vC20, vC18);
        Triad<Apfloat>        tri6_norm  = normalTriple(vC3,  vC20, vC18, true);
        Triad<Tuple<Apfloat>> tri7  = new Triad<>(vC4,  vC22, vC23);
        Triad<Apfloat>        tri7_norm  = normalTriple(vC4,  vC22, vC23, true);
        Triad<Tuple<Apfloat>> tri8  = new Triad<>(vC4,  vC27, vC26);
        Triad<Apfloat>        tri8_norm  = normalTriple(vC4,  vC27, vC26, true);
        Triad<Tuple<Apfloat>> tri9  = new Triad<>(vC5,  vC25, vC24);
        Triad<Apfloat>        tri9_norm  = normalTriple(vC5,  vC25, vC24, true);
        Triad<Tuple<Apfloat>> tri10 = new Triad<>(vC5,  vC28, vC29);
        Triad<Apfloat>        tri10_norm = normalTriple(vC5,  vC28, vC29, true);
        Triad<Tuple<Apfloat>> tri11 = new Triad<>(vC6,  vC14, vC22);
        Triad<Apfloat>        tri11_norm = normalTriple(vC6,  vC14, vC22, true);
        Triad<Tuple<Apfloat>> tri12 = new Triad<>(vC7,  vC23, vC15);
        Triad<Apfloat>        tri12_norm = normalTriple(vC7,  vC23, vC15, true);
        Triad<Tuple<Apfloat>> tri13 = new Triad<>(vC8,  vC24, vC16);
        Triad<Apfloat>        tri13_norm = normalTriple(vC8,  vC24, vC16, true);
        Triad<Tuple<Apfloat>> tri14 = new Triad<>(vC9,  vC17, vC25);
        Triad<Apfloat>        tri14_norm = normalTriple(vC9,  vC17, vC25, true);
        Triad<Tuple<Apfloat>> tri15 = new Triad<>(vC10, vC26, vC18);
        Triad<Apfloat>        tri15_norm = normalTriple(vC10, vC26, vC18, true);
        Triad<Tuple<Apfloat>> tri16 = new Triad<>(vC11, vC19, vC27);
        Triad<Apfloat>        tri16_norm = normalTriple(vC11, vC19, vC27, true);
        Triad<Tuple<Apfloat>> tri17 = new Triad<>(vC12, vC20, vC28);
        Triad<Apfloat>        tri17_norm = normalTriple(vC12, vC20, vC28, true);
        Triad<Tuple<Apfloat>> tri18 = new Triad<>(vC13, vC29, vC21);
        Triad<Apfloat>        tri18_norm = normalTriple(vC13, vC29, vC21, true);
        Triad<Tuple<Apfloat>> tri19 = new Triad<>(vC0, vC6, vC10);
        Triad<Apfloat> tri19_norm = normalTriple(vC0, vC6, vC10, true);
        faces_tri = new Icosad<>(
            tri0, tri1, tri2, tri3,
            tri4, tri5, tri6, tri7,
            tri8, tri9, tri10, tri11,
            tri12, tri13, tri14, tri15,
            tri16, tri17, tri18, tri19
        );
        face_norms_tri = new Icosad<>(
            tri0_norm, tri1_norm, tri2_norm, tri3_norm,
            tri4_norm, tri5_norm, tri6_norm, tri7_norm,
            tri8_norm, tri9_norm, tri10_norm, tri11_norm,
            tri12_norm, tri13_norm, tri14_norm, tri15_norm,
            tri16_norm, tri17_norm, tri18_norm, tri19_norm
        );

        // ==== PENTAGONAL FACES ====
        Pentad<Tuple<Apfloat>> pnt0 = new Pentad<>(vC0, vC8, vC16, vC14, vC6);
        Triad<Apfloat> pnt0_norm = normalPent(vC0, vC8, vC16, vC14, vC6, true);
        Pentad<Tuple<Apfloat>> pnt1  = new Pentad<>(vC1,  vC7,  vC15, vC17, vC9);
        Triad<Apfloat>         pnt1_norm  = normalPent(vC1,  vC7,  vC15, vC17, vC9, true);
        Pentad<Tuple<Apfloat>> pnt2  = new Pentad<>(vC1,  vC13, vC21, vC19, vC11);
        Triad<Apfloat>         pnt2_norm  = normalPent(vC1,  vC13, vC21, vC19, vC11, true);
        Pentad<Tuple<Apfloat>> pnt3  = new Pentad<>(vC2,  vC15, vC23, vC22, vC14);
        Triad<Apfloat>         pnt3_norm  = normalPent(vC2,  vC15, vC23, vC22, vC14, true);
        Pentad<Tuple<Apfloat>> pnt4  = new Pentad<>(vC2,  vC16, vC24, vC25, vC17);
        Triad<Apfloat>         pnt4_norm  = normalPent(vC2,  vC16, vC24, vC25, vC17, true);
        Pentad<Tuple<Apfloat>> pnt5  = new Pentad<>(vC3,  vC18, vC26, vC27, vC19);
        Triad<Apfloat>         pnt5_norm  = normalPent(vC3,  vC18, vC26, vC27, vC19, true);
        Pentad<Tuple<Apfloat>> pnt6  = new Pentad<>(vC3,  vC21, vC29, vC28, vC20);
        Triad<Apfloat>         pnt6_norm  = normalPent(vC3,  vC21, vC29, vC28, vC20, true);
        Pentad<Tuple<Apfloat>> pnt7  = new Pentad<>(vC4,  vC23, vC7,  vC11, vC27);
        Triad<Apfloat>         pnt7_norm  = normalPent(vC4,  vC23, vC7,  vC11, vC27, true);
        Pentad<Tuple<Apfloat>> pnt8  = new Pentad<>(vC4,  vC26, vC10, vC6,  vC22);
        Triad<Apfloat>         pnt8_norm  = normalPent(vC4,  vC26, vC10, vC6,  vC22, true);
        Pentad<Tuple<Apfloat>> pnt9  = new Pentad<>(vC5,  vC24, vC8,  vC12, vC28);
        Triad<Apfloat>         pnt9_norm  = normalPent(vC5,  vC24, vC8,  vC12, vC28, true);
        Pentad<Tuple<Apfloat>> pnt10 = new Pentad<>(vC5,  vC29, vC13, vC9,  vC25);
        Triad<Apfloat>         pnt10_norm = normalPent(vC5,  vC29, vC13, vC9,  vC25, true);
        Pentad<Tuple<Apfloat>> pnt11 = new Pentad<>(vC0,  vC10, vC18, vC20,  vC12);
        Triad<Apfloat>         pnt11_norm = normalPent(vC0,  vC10, vC18, vC20,  vC12, true);

        faces_pnt = new Dodecad<>(
            pnt0, pnt1, pnt2,
            pnt3, pnt4, pnt5,
            pnt6, pnt7, pnt8,
            pnt9, pnt10, pnt11
        );
        face_norms_pnt = new Dodecad<>(
            pnt0_norm, pnt1_norm, pnt2_norm,
            pnt3_norm, pnt4_norm, pnt5_norm,
            pnt6_norm, pnt7_norm, pnt8_norm,
            pnt9_norm, pnt10_norm, pnt11_norm
        );
    }

    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_tri.fetchSize(); i++){
            // === Check pentagonal faces ===
            if (i < faces_pnt.fetchSize()) {
                Pentad<Tuple<Apfloat>> face = (Pentad<Tuple<Apfloat>>) faces_pnt.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //  = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //  = face_norm_x*(p_x-vertA_x)
                //  + face_norm_y*(p_y-vertA_y)
                //  + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_pnt.fetch(i);
                Apfloat d = dot_prod(norm, m);

                  // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                    }
            }

            // === Check triangular faces ===
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
            Triad<Apfloat> m = subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //  = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //  = face_norm_x*(p_x-vertA_x)
            //  + face_norm_y*(p_y-vertA_y)
            //  + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.fetch(i);
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
