package Shapes;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;
import Lattice.*;
import Atom.*;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/Icosidodecahedron.txt>...</a>}
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
        Triad<Apfloat> vC0 = VectorMath.mult(VectorMath.normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1 = VectorMath.mult(VectorMath.normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2 = VectorMath.mult(VectorMath.normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3 = VectorMath.mult(VectorMath.normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4 = VectorMath.mult(VectorMath.normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5 = VectorMath.mult(VectorMath.normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6 = VectorMath.mult(VectorMath.normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7 = VectorMath.mult(VectorMath.normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8 = VectorMath.mult(VectorMath.normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9 = VectorMath.mult(VectorMath.normalize(vB9), super.getRadius().toString());
        Triad<Apfloat> vC10 = VectorMath.mult(VectorMath.normalize(vB10), super.getRadius().toString());
        Triad<Apfloat> vC11 = VectorMath.mult(VectorMath.normalize(vB11), super.getRadius().toString());
        Triad<Apfloat> vC12 = VectorMath.mult(VectorMath.normalize(vB12), super.getRadius().toString());
        Triad<Apfloat> vC13 = VectorMath.mult(VectorMath.normalize(vB13), super.getRadius().toString());
        Triad<Apfloat> vC14 = VectorMath.mult(VectorMath.normalize(vB14), super.getRadius().toString());
        Triad<Apfloat> vC15 = VectorMath.mult(VectorMath.normalize(vB15), super.getRadius().toString());
        Triad<Apfloat> vC16 = VectorMath.mult(VectorMath.normalize(vB16), super.getRadius().toString());
        Triad<Apfloat> vC17 = VectorMath.mult(VectorMath.normalize(vB17), super.getRadius().toString());
        Triad<Apfloat> vC18 = VectorMath.mult(VectorMath.normalize(vB18), super.getRadius().toString());
        Triad<Apfloat> vC19 = VectorMath.mult(VectorMath.normalize(vB19), super.getRadius().toString());
        Triad<Apfloat> vC20 = VectorMath.mult(VectorMath.normalize(vB20), super.getRadius().toString());
        Triad<Apfloat> vC21 = VectorMath.mult(VectorMath.normalize(vB21), super.getRadius().toString());
        Triad<Apfloat> vC22 = VectorMath.mult(VectorMath.normalize(vB22), super.getRadius().toString());
        Triad<Apfloat> vC23 = VectorMath.mult(VectorMath.normalize(vB23), super.getRadius().toString());
        Triad<Apfloat> vC24 = VectorMath.mult(VectorMath.normalize(vB24), super.getRadius().toString());
        Triad<Apfloat> vC25 = VectorMath.mult(VectorMath.normalize(vB25), super.getRadius().toString());
        Triad<Apfloat> vC26 = VectorMath.mult(VectorMath.normalize(vB26), super.getRadius().toString());
        Triad<Apfloat> vC27 = VectorMath.mult(VectorMath.normalize(vB27), super.getRadius().toString());
        Triad<Apfloat> vC28 = VectorMath.mult(VectorMath.normalize(vB28), super.getRadius().toString());
        Triad<Apfloat> vC29 = VectorMath.mult(VectorMath.normalize(vB29), super.getRadius().toString());

        // ==== Triangular FACES ====
        Triad<Tuple<Apfloat>> tri0  = new Triad<>(vC0,  vC12, vC8);
        Triad<Apfloat>        tri0_norm  = VectorMath.normalTriple(vC0,  vC12, vC8,  true);
        Triad<Tuple<Apfloat>> tri1  = new Triad<>(vC1,  vC9,  vC13);
        Triad<Apfloat>        tri1_norm  = VectorMath.normalTriple(vC1,  vC9,  vC13, true);
        Triad<Tuple<Apfloat>> tri2  = new Triad<>(vC1,  vC11, vC7);
        Triad<Apfloat>        tri2_norm  = VectorMath.normalTriple(vC1,  vC11, vC7,  true);
        Triad<Tuple<Apfloat>> tri3  = new Triad<>(vC2,  vC14, vC16);
        Triad<Apfloat>        tri3_norm  = VectorMath.normalTriple(vC2,  vC14, vC16, true);
        Triad<Tuple<Apfloat>> tri4  = new Triad<>(vC2,  vC17, vC15);
        Triad<Apfloat>        tri4_norm  = VectorMath.normalTriple(vC2,  vC17, vC15, true);
        Triad<Tuple<Apfloat>> tri5  = new Triad<>(vC3,  vC19, vC21);
        Triad<Apfloat>        tri5_norm  = VectorMath.normalTriple(vC3,  vC19, vC21, true);
        Triad<Tuple<Apfloat>> tri6  = new Triad<>(vC3,  vC20, vC18);
        Triad<Apfloat>        tri6_norm  = VectorMath.normalTriple(vC3,  vC20, vC18, true);
        Triad<Tuple<Apfloat>> tri7  = new Triad<>(vC4,  vC22, vC23);
        Triad<Apfloat>        tri7_norm  = VectorMath.normalTriple(vC4,  vC22, vC23, true);
        Triad<Tuple<Apfloat>> tri8  = new Triad<>(vC4,  vC27, vC26);
        Triad<Apfloat>        tri8_norm  = VectorMath.normalTriple(vC4,  vC27, vC26, true);
        Triad<Tuple<Apfloat>> tri9  = new Triad<>(vC5,  vC25, vC24);
        Triad<Apfloat>        tri9_norm  = VectorMath.normalTriple(vC5,  vC25, vC24, true);
        Triad<Tuple<Apfloat>> tri10 = new Triad<>(vC5,  vC28, vC29);
        Triad<Apfloat>        tri10_norm = VectorMath.normalTriple(vC5,  vC28, vC29, true);
        Triad<Tuple<Apfloat>> tri11 = new Triad<>(vC6,  vC14, vC22);
        Triad<Apfloat>        tri11_norm = VectorMath.normalTriple(vC6,  vC14, vC22, true);
        Triad<Tuple<Apfloat>> tri12 = new Triad<>(vC7,  vC23, vC15);
        Triad<Apfloat>        tri12_norm = VectorMath.normalTriple(vC7,  vC23, vC15, true);
        Triad<Tuple<Apfloat>> tri13 = new Triad<>(vC8,  vC24, vC16);
        Triad<Apfloat>        tri13_norm = VectorMath.normalTriple(vC8,  vC24, vC16, true);
        Triad<Tuple<Apfloat>> tri14 = new Triad<>(vC9,  vC17, vC25);
        Triad<Apfloat>        tri14_norm = VectorMath.normalTriple(vC9,  vC17, vC25, true);
        Triad<Tuple<Apfloat>> tri15 = new Triad<>(vC10, vC26, vC18);
        Triad<Apfloat>        tri15_norm = VectorMath.normalTriple(vC10, vC26, vC18, true);
        Triad<Tuple<Apfloat>> tri16 = new Triad<>(vC11, vC19, vC27);
        Triad<Apfloat>        tri16_norm = VectorMath.normalTriple(vC11, vC19, vC27, true);
        Triad<Tuple<Apfloat>> tri17 = new Triad<>(vC12, vC20, vC28);
        Triad<Apfloat>        tri17_norm = VectorMath.normalTriple(vC12, vC20, vC28, true);
        Triad<Tuple<Apfloat>> tri18 = new Triad<>(vC13, vC29, vC21);
        Triad<Apfloat>        tri18_norm = VectorMath.normalTriple(vC13, vC29, vC21, true);
        Triad<Tuple<Apfloat>> tri19 = new Triad<>(vC0, vC6, vC10);
        Triad<Apfloat> tri19_norm = VectorMath.normalTriple(vC0, vC6, vC10, true);
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
        Triad<Apfloat> pnt0_norm = VectorMath.normalPent(vC0, vC8, vC16, vC14, vC6, true);
        Pentad<Tuple<Apfloat>> pnt1  = new Pentad<>(vC1,  vC7,  vC15, vC17, vC9);
        Triad<Apfloat>         pnt1_norm  = VectorMath.normalPent(vC1,  vC7,  vC15, vC17, vC9, true);
        Pentad<Tuple<Apfloat>> pnt2  = new Pentad<>(vC1,  vC13, vC21, vC19, vC11);
        Triad<Apfloat>         pnt2_norm  = VectorMath.normalPent(vC1,  vC13, vC21, vC19, vC11, true);
        Pentad<Tuple<Apfloat>> pnt3  = new Pentad<>(vC2,  vC15, vC23, vC22, vC14);
        Triad<Apfloat>         pnt3_norm  = VectorMath.normalPent(vC2,  vC15, vC23, vC22, vC14, true);
        Pentad<Tuple<Apfloat>> pnt4  = new Pentad<>(vC2,  vC16, vC24, vC25, vC17);
        Triad<Apfloat>         pnt4_norm  = VectorMath.normalPent(vC2,  vC16, vC24, vC25, vC17, true);
        Pentad<Tuple<Apfloat>> pnt5  = new Pentad<>(vC3,  vC18, vC26, vC27, vC19);
        Triad<Apfloat>         pnt5_norm  = VectorMath.normalPent(vC3,  vC18, vC26, vC27, vC19, true);
        Pentad<Tuple<Apfloat>> pnt6  = new Pentad<>(vC3,  vC21, vC29, vC28, vC20);
        Triad<Apfloat>         pnt6_norm  = VectorMath.normalPent(vC3,  vC21, vC29, vC28, vC20, true);
        Pentad<Tuple<Apfloat>> pnt7  = new Pentad<>(vC4,  vC23, vC7,  vC11, vC27);
        Triad<Apfloat>         pnt7_norm  = VectorMath.normalPent(vC4,  vC23, vC7,  vC11, vC27, true);
        Pentad<Tuple<Apfloat>> pnt8  = new Pentad<>(vC4,  vC26, vC10, vC6,  vC22);
        Triad<Apfloat>         pnt8_norm  = VectorMath.normalPent(vC4,  vC26, vC10, vC6,  vC22, true);
        Pentad<Tuple<Apfloat>> pnt9  = new Pentad<>(vC5,  vC24, vC8,  vC12, vC28);
        Triad<Apfloat>         pnt9_norm  = VectorMath.normalPent(vC5,  vC24, vC8,  vC12, vC28, true);
        Pentad<Tuple<Apfloat>> pnt10 = new Pentad<>(vC5,  vC29, vC13, vC9,  vC25);
        Triad<Apfloat>         pnt10_norm = VectorMath.normalPent(vC5,  vC29, vC13, vC9,  vC25, true);
        Pentad<Tuple<Apfloat>> pnt11 = new Pentad<>(vC0,  vC10, vC18, vC20,  vC12);
        Triad<Apfloat>         pnt11_norm = VectorMath.normalPent(vC0,  vC10, vC18, vC20,  vC12, true);

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
                Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //  = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //  = face_norm_x*(p_x-vertA_x)
                //  + face_norm_y*(p_y-vertA_y)
                //  + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_pnt.fetch(i);
                Apfloat d = VectorMath.dot_prod(norm, m);

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
            Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //  = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //  = face_norm_x*(p_x-vertA_x)
            //  + face_norm_y*(p_y-vertA_y)
            //  + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.fetch(i);
            Apfloat d = VectorMath.dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
