package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href="https://dmccooey.com/polyhedra/CanonicalTruncatedTruncatedTetrahedron.html">...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class TetrahedronDitruncated extends Shape {

    // 3 dodecahedral faces
    private final Triad<Tuple<Tuple<Apfloat>>> faces_dod;
    private final Triad<Tuple<Apfloat>> face_norms_dod;

    // 4 hexagonal faces
    private final Tetrad<Tuple<Tuple<Apfloat>>> faces_hex;
    private final Tetrad<Tuple<Apfloat>> face_norms_hex;

    // 12 triangular faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Dodecad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N9 = new Apfloat("9", super.precision);
    private final Apfloat N7 = new Apfloat("7", super.precision);
    private final Apfloat FRAC_N5_over_N9 = N5.divide(N9);
    private final Apfloat FRAC_N7_over_N9 = N7.divide(N9);
    private final Apfloat FRAC_N1_over_N3 = N1.divide(N3);
    private final Apfloat FRAC_N1_over_N5 = N1.divide(N5);

    public TetrahedronDitruncated(
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
        Triad<Apfloat> vB0 = new Triad<>(FRAC_N1_over_N5, FRAC_N1_over_N5.multiply(NEG_N1), N1);                                            // (1/5,-1/5,1)
        Triad<Apfloat> vB1 = new Triad<>(FRAC_N1_over_N5, FRAC_N1_over_N5, NEG_N1);                                                          // (1/5,1/5,-1)
        Triad<Apfloat> vB2 = new Triad<>(FRAC_N1_over_N5.multiply(NEG_N1), FRAC_N1_over_N5, N1);                                            // (-1/5,1/5,1)
        Triad<Apfloat> vB3  = new Triad<>(FRAC_N1_over_N5.multiply(NEG_N1), FRAC_N1_over_N5.multiply(NEG_N1), NEG_N1);                             // (-1/5,-1/5,-1)
        Triad<Apfloat> vB4  = new Triad<>(N1, FRAC_N1_over_N5.multiply(NEG_N1), FRAC_N1_over_N5);                                           // (+1,-1/5,+1/5)
        Triad<Apfloat> vB5  = new Triad<>(N1, FRAC_N1_over_N5, FRAC_N1_over_N5.multiply(NEG_N1));                                           // (+1,+1/5,-1/5)
        Triad<Apfloat> vB6  = new Triad<>(N1.multiply(NEG_N1), FRAC_N1_over_N5, FRAC_N1_over_N5);                                           // (-1,+1/5,+1/5)
        Triad<Apfloat> vB7  = new Triad<>(N1.multiply(NEG_N1), FRAC_N1_over_N5.multiply(NEG_N1), FRAC_N1_over_N5.multiply(NEG_N1));               // (-1,-1/5,-1/5)
        Triad<Apfloat> vB8  = new Triad<>(FRAC_N1_over_N5, N1.multiply(NEG_N1), FRAC_N1_over_N5);                                           // (+1/5,-1,+1/5)
        Triad<Apfloat> vB9  = new Triad<>(FRAC_N1_over_N5, N1, FRAC_N1_over_N5.multiply(NEG_N1));                                           // (+1/5,+1,-1/5)
        Triad<Apfloat> vB10 = new Triad<>(FRAC_N1_over_N5.multiply(NEG_N1), N1, FRAC_N1_over_N5);                                           // (-1/5,+1,+1/5)
        Triad<Apfloat> vB11 = new Triad<>(FRAC_N1_over_N5.multiply(NEG_N1), N1.multiply(NEG_N1), FRAC_N1_over_N5.multiply(NEG_N1));               // (-1/5,-1,-1/5)
        Triad<Apfloat> vB12 = new Triad<>(FRAC_N5_over_N9, FRAC_N1_over_N3.multiply(NEG_N1), FRAC_N7_over_N9);                              // (+5/9,-1/3,+7/9)
        Triad<Apfloat> vB13 = new Triad<>(FRAC_N5_over_N9, FRAC_N1_over_N3, FRAC_N7_over_N9.multiply(NEG_N1));                              // (+5/9,+1/3,-7/9)
        Triad<Apfloat> vB14 = new Triad<>(FRAC_N5_over_N9.multiply(NEG_N1), FRAC_N1_over_N3, FRAC_N7_over_N9);                              // (-5/9,+1/3,+7/9)
        Triad<Apfloat> vB15 = new Triad<>(FRAC_N5_over_N9.multiply(NEG_N1), FRAC_N1_over_N3.multiply(NEG_N1), FRAC_N7_over_N9.multiply(NEG_N1));  // (-5/9,-1/3,-7/9)
        Triad<Apfloat> vB16 = new Triad<>(FRAC_N7_over_N9, FRAC_N5_over_N9.multiply(NEG_N1), FRAC_N1_over_N3);                              // (+7/9,-5/9,+1/3)
        Triad<Apfloat> vB17 = new Triad<>(FRAC_N7_over_N9, FRAC_N5_over_N9, FRAC_N1_over_N3.multiply(NEG_N1));                              // (+7/9,+5/9,-1/3)
        Triad<Apfloat> vB18 = new Triad<>(FRAC_N7_over_N9.multiply(NEG_N1), FRAC_N5_over_N9, FRAC_N1_over_N3);                              // (-7/9,+5/9,+1/3)
        Triad<Apfloat> vB19 = new Triad<>(FRAC_N7_over_N9.multiply(NEG_N1), FRAC_N5_over_N9.multiply(NEG_N1), FRAC_N1_over_N3.multiply(NEG_N1));  // (-7/9,-5/9,-1/3)
        Triad<Apfloat> vB20 = new Triad<>(FRAC_N1_over_N3, FRAC_N7_over_N9.multiply(NEG_N1), FRAC_N5_over_N9);                              // (+1/3,-7/9,+5/9)
        Triad<Apfloat> vB21 = new Triad<>(FRAC_N1_over_N3, FRAC_N7_over_N9, FRAC_N5_over_N9.multiply(NEG_N1));                              // (+1/3,+7/9,-5/9)
        Triad<Apfloat> vB22 = new Triad<>(FRAC_N1_over_N3.multiply(NEG_N1), FRAC_N7_over_N9, FRAC_N5_over_N9);                              // (-1/3,+7/9,+5/9)
        Triad<Apfloat> vB23 = new Triad<>(FRAC_N1_over_N3.multiply(NEG_N1), FRAC_N7_over_N9.multiply(NEG_N1), FRAC_N5_over_N9.multiply(NEG_N1));  // (-1/3,-7/9,-5/9)
        Triad<Apfloat> vB24 = new Triad<>(FRAC_N1_over_N3, FRAC_N5_over_N9.multiply(NEG_N1), FRAC_N7_over_N9);                              // (+1/3,-5/9,+7/9)
        Triad<Apfloat> vB25 = new Triad<>(FRAC_N1_over_N3, FRAC_N5_over_N9, FRAC_N7_over_N9.multiply(NEG_N1));                              // (+1/3,+5/9,-7/9)
        Triad<Apfloat> vB26 = new Triad<>(FRAC_N1_over_N3.multiply(NEG_N1), FRAC_N5_over_N9, FRAC_N7_over_N9);                              // (-1/3,+5/9,+7/9)
        Triad<Apfloat> vB27 = new Triad<>(FRAC_N1_over_N3.multiply(NEG_N1), FRAC_N5_over_N9.multiply(NEG_N1), FRAC_N7_over_N9.multiply(NEG_N1));  // (-1/3,-5/9,-7/9)
        Triad<Apfloat> vB28 = new Triad<>(FRAC_N7_over_N9, FRAC_N1_over_N3.multiply(NEG_N1), FRAC_N5_over_N9);                              // (+7/9,-1/3,+5/9)
        Triad<Apfloat> vB29 = new Triad<>(FRAC_N7_over_N9, FRAC_N1_over_N3, FRAC_N5_over_N9.multiply(NEG_N1));                              // (+7/9,+1/3,-5/9)
        Triad<Apfloat> vB30 = new Triad<>(FRAC_N7_over_N9.multiply(NEG_N1), FRAC_N1_over_N3, FRAC_N5_over_N9);                              // (-7/9,+1/3,+5/9)
        Triad<Apfloat> vB31 = new Triad<>(FRAC_N7_over_N9.multiply(NEG_N1), FRAC_N1_over_N3.multiply(NEG_N1), FRAC_N5_over_N9.multiply(NEG_N1));  // (-7/9,-1/3,-5/9)
        Triad<Apfloat> vB32 = new Triad<>(FRAC_N5_over_N9, FRAC_N7_over_N9.multiply(NEG_N1), FRAC_N1_over_N3);                              // (+5/9,-7/9,+1/3)
        Triad<Apfloat> vB33 = new Triad<>(FRAC_N5_over_N9, FRAC_N7_over_N9, FRAC_N1_over_N3.multiply(NEG_N1));                              // (+5/9,+7/9,-1/3)
        Triad<Apfloat> vB34 = new Triad<>(FRAC_N5_over_N9.multiply(NEG_N1), FRAC_N7_over_N9, FRAC_N1_over_N3);                              // (-5/9,+7/9,+1/3)
        Triad<Apfloat> vB35 = new Triad<>(FRAC_N5_over_N9.multiply(NEG_N1), FRAC_N7_over_N9.multiply(NEG_N1), FRAC_N1_over_N3.multiply(NEG_N1));  // (-5/9,-7/9,-1/3)

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0 = VectorMath.mult(VectorMath.normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1 = VectorMath.mult(VectorMath.normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2 = VectorMath.mult(VectorMath.normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3  = VectorMath.mult(VectorMath.normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4  = VectorMath.mult(VectorMath.normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5  = VectorMath.mult(VectorMath.normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6  = VectorMath.mult(VectorMath.normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7  = VectorMath.mult(VectorMath.normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8  = VectorMath.mult(VectorMath.normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9  = VectorMath.mult(VectorMath.normalize(vB9), super.getRadius().toString());
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
        Triad<Apfloat> vC30 = VectorMath.mult(VectorMath.normalize(vB30), super.getRadius().toString());
        Triad<Apfloat> vC31 = VectorMath.mult(VectorMath.normalize(vB31), super.getRadius().toString());
        Triad<Apfloat> vC32 = VectorMath.mult(VectorMath.normalize(vB32), super.getRadius().toString());
        Triad<Apfloat> vC33 = VectorMath.mult(VectorMath.normalize(vB33), super.getRadius().toString());
        Triad<Apfloat> vC34 = VectorMath.mult(VectorMath.normalize(vB34), super.getRadius().toString());
        Triad<Apfloat> vC35 = VectorMath.mult(VectorMath.normalize(vB35), super.getRadius().toString());

        // ==== DODECAGONAL FACES ====
        Dodecad<Tuple<Apfloat>> dod0 = new Dodecad<>(vC12, vC28, vC4, vC5, vC17, vC33, vC9, vC10, vC22, vC26, vC2, vC0);
        Triad<Apfloat> dod0_norm = VectorMath.normalDodeca(vC12, vC28, vC4, vC5, vC17, vC33, vC9, vC10, vC22, vC26, vC2, vC0, true);
        Dodecad<Tuple<Apfloat>> dod1 = new Dodecad<>(vC12, vC29, vC5, vC4, vC16, vC32, vC8, vC11, vC23, vC27, vC3, vC1);
        Triad<Apfloat> dod1_norm = VectorMath.normalDodeca(vC12, vC29, vC5, vC4, vC16, vC32, vC8, vC11, vC23, vC27, vC3, vC1, true);
        Dodecad<Tuple<Apfloat>> dod2 = new Dodecad<>(vC15, vC31, vC7, vC6, vC18, vC34, vC10, vC9, vC21, vC25, vC1, vC3);
        Triad<Apfloat> dod2_norm = VectorMath.normalDodeca(vC15, vC31, vC7, vC6, vC18, vC34, vC10, vC9, vC21, vC25, vC1, vC3, true);
        faces_dod = new Triad<>(dod0, dod1, dod2);
        face_norms_dod = new Triad<>(dod0_norm, dod1_norm, dod2_norm);

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC12, vC24, vC20, vC32, vC16, vC28);
        Triad<Apfloat> hex0_norm = VectorMath.normalHex(vC12, vC24, vC20, vC32, vC16, vC28,true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC13, vC25, vC21, vC33, vC17, vC29);
        Triad<Apfloat> hex1_norm = VectorMath.normalHex(vC13, vC25, vC21, vC33, vC17, vC29,true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC14, vC26, vC22, vC34, vC18, vC30);
        Triad<Apfloat> hex2_norm = VectorMath.normalHex(vC14, vC26, vC22, vC34, vC18, vC30,true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC15, vC27, vC23, vC35, vC19, vC31);
        Triad<Apfloat> hex3_norm = VectorMath.normalHex(vC15, vC27, vC23, vC35, vC19, vC31,true);
        faces_hex = new Tetrad<>(hex0, hex1, hex2, hex3);
        face_norms_hex = new Tetrad<>(hex0_norm, hex1_norm, hex2_norm, hex3_norm);

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0  = new Triad<>(vC12, vC0,  vC24);
        Triad<Apfloat>        tri0_norm  = VectorMath.normalTriple(vC12, vC0,  vC24, true);
        Triad<Tuple<Apfloat>> tri1  = new Triad<>(vC13, vC1,  vC25);
        Triad<Apfloat>        tri1_norm  = VectorMath.normalTriple(vC13, vC1,  vC25, true);
        Triad<Tuple<Apfloat>> tri2  = new Triad<>(vC14, vC2,  vC26);
        Triad<Apfloat>        tri2_norm  = VectorMath.normalTriple(vC14, vC2,  vC26, true);
        Triad<Tuple<Apfloat>> tri3  = new Triad<>(vC15, vC3,  vC27);
        Triad<Apfloat>        tri3_norm  = VectorMath.normalTriple(vC15, vC3,  vC27, true);
        Triad<Tuple<Apfloat>> tri4  = new Triad<>(vC16, vC4,  vC28);
        Triad<Apfloat>        tri4_norm  = VectorMath.normalTriple(vC16, vC4,  vC28, true);
        Triad<Tuple<Apfloat>> tri5  = new Triad<>(vC17, vC5,  vC29);
        Triad<Apfloat>        tri5_norm  = VectorMath.normalTriple(vC17, vC5,  vC29, true);
        Triad<Tuple<Apfloat>> tri6  = new Triad<>(vC18, vC6,  vC30);
        Triad<Apfloat>        tri6_norm  = VectorMath.normalTriple(vC18, vC6,  vC30, true);
        Triad<Tuple<Apfloat>> tri7  = new Triad<>(vC19, vC7,  vC31);
        Triad<Apfloat>        tri7_norm  = VectorMath.normalTriple(vC19, vC7,  vC31, true);
        Triad<Tuple<Apfloat>> tri8  = new Triad<>(vC20, vC8,  vC32);
        Triad<Apfloat>        tri8_norm  = VectorMath.normalTriple(vC20, vC8,  vC32, true);
        Triad<Tuple<Apfloat>> tri9  = new Triad<>(vC21, vC9,  vC33);
        Triad<Apfloat>        tri9_norm  = VectorMath.normalTriple(vC21, vC9,  vC33, true);
        Triad<Tuple<Apfloat>> tri10 = new Triad<>(vC22, vC10, vC34);
        Triad<Apfloat>        tri10_norm = VectorMath.normalTriple(vC22, vC10, vC34, true);
        Triad<Tuple<Apfloat>> tri11 = new Triad<>(vC23, vC11, vC35);
        Triad<Apfloat>        tri11_norm = VectorMath.normalTriple(vC23, vC11, vC35, true);
        faces_tri = new Dodecad<>(
                tri0, tri1, tri2, tri3,
                tri4, tri5, tri6, tri7,
                tri8, tri9, tri10, tri11
        );
        face_norms_tri = new Dodecad<>(
                tri0_norm, tri1_norm, tri2_norm, tri3_norm,
                tri4_norm, tri5_norm, tri6_norm, tri7_norm,
                tri8_norm, tri9_norm, tri10_norm, tri11_norm
        );

    }


    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < faces_tri.fetchSize(); i++){
            // === Check hexagonal faces ===
            if (i < face_norms_hex.fetchSize()) {
                Hexad<Tuple<Apfloat>> face = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_hex.fetch(i);
                Apfloat d = VectorMath.dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check dodecagonal faces ===
            if (i < faces_dod.fetchSize()) {
                Dodecad<Tuple<Apfloat>> face = (Dodecad<Tuple<Apfloat>>) faces_dod.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_dod.fetch(i);
                Apfloat d = VectorMath.dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check triangular faces ===
            // check if each point lies within bounds of each face

            // get vertices of vertA
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.fetch(i);
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
            Triad<Apfloat> face_norm = (Triad<Apfloat>) face_norms_tri.fetch(i);
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
