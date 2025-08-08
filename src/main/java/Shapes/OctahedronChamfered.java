package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/ChamferedOctahedron3.html>...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class OctahedronChamfered extends Shape {

    //  dodecahedral faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_hex;
    private final Dodecad<Tuple<Apfloat>> face_norms_hex;

    // 8 hexagonal faces
    private final Octad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Octad<Tuple<Apfloat>> face_norms_tri;


    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N6 = new Apfloat("6", super.precision);
    private final Apfloat FRAC_N1_over_N2 = N1.divide(N2);
    private final Apfloat SQRT6minFRAC_N1_over_N2 = (ApfloatMath.sqrt(N6).subtract(N1)).divide(N2);
    private final Apfloat FRAC_SQRT6_over_N2 = ApfloatMath.sqrt(N6).divide(N2);

    public OctahedronChamfered(
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
        Triad<Apfloat> vB0 = new Triad<>(N0, N0, SQRT6minFRAC_N1_over_N2);
        Triad<Apfloat> vB1 = new Triad<>(N0, N0, SQRT6minFRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB2  = new Triad<>(SQRT6minFRAC_N1_over_N2, N0, N0);
        Triad<Apfloat> vB3  = new Triad<>(SQRT6minFRAC_N1_over_N2.multiply(NEG_N1), N0, N0);
        Triad<Apfloat> vB4  = new Triad<>(N0, SQRT6minFRAC_N1_over_N2, N0);
        Triad<Apfloat> vB5  = new Triad<>(N0, SQRT6minFRAC_N1_over_N2.multiply(NEG_N1), N0);
        Triad<Apfloat> vB6  = new Triad<>(FRAC_N1_over_N2, FRAC_N1_over_N2, FRAC_SQRT6_over_N2);
        Triad<Apfloat> vB7  = new Triad<>(FRAC_N1_over_N2, FRAC_N1_over_N2, FRAC_SQRT6_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB8  = new Triad<>(FRAC_N1_over_N2, FRAC_N1_over_N2.multiply(NEG_N1), FRAC_SQRT6_over_N2);
        Triad<Apfloat> vB9  = new Triad<>(FRAC_N1_over_N2, FRAC_N1_over_N2.multiply(NEG_N1), FRAC_SQRT6_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), FRAC_N1_over_N2, FRAC_SQRT6_over_N2);
        Triad<Apfloat> vB11 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), FRAC_N1_over_N2, FRAC_SQRT6_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1), FRAC_SQRT6_over_N2);
        Triad<Apfloat> vB13 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1), FRAC_SQRT6_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(FRAC_SQRT6_over_N2, FRAC_N1_over_N2, FRAC_N1_over_N2);
        Triad<Apfloat> vB15 = new Triad<>(FRAC_SQRT6_over_N2, FRAC_N1_over_N2, FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>(FRAC_SQRT6_over_N2, FRAC_N1_over_N2.multiply(NEG_N1), FRAC_N1_over_N2);
        Triad<Apfloat> vB17 = new Triad<>(FRAC_SQRT6_over_N2, FRAC_N1_over_N2.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>(FRAC_SQRT6_over_N2.multiply(NEG_N1), FRAC_N1_over_N2, FRAC_N1_over_N2);
        Triad<Apfloat> vB19 = new Triad<>(FRAC_SQRT6_over_N2.multiply(NEG_N1), FRAC_N1_over_N2, FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>(FRAC_SQRT6_over_N2.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1), FRAC_N1_over_N2);
        Triad<Apfloat> vB21 = new Triad<>(FRAC_SQRT6_over_N2.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB22 = new Triad<>(FRAC_N1_over_N2, FRAC_SQRT6_over_N2, FRAC_N1_over_N2);
        Triad<Apfloat> vB23 = new Triad<>(FRAC_N1_over_N2, FRAC_SQRT6_over_N2, FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB24 = new Triad<>(FRAC_N1_over_N2, FRAC_SQRT6_over_N2.multiply(NEG_N1), FRAC_N1_over_N2);
        Triad<Apfloat> vB25 = new Triad<>(FRAC_N1_over_N2, FRAC_SQRT6_over_N2.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB26 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), FRAC_SQRT6_over_N2, FRAC_N1_over_N2);
        Triad<Apfloat> vB27 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), FRAC_SQRT6_over_N2, FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB28 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), FRAC_SQRT6_over_N2.multiply(NEG_N1), FRAC_N1_over_N2);
        Triad<Apfloat> vB29 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), FRAC_SQRT6_over_N2.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1));

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

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0  = new Hexad<>(vC0, vC6,  vC22, vC4,  vC26, vC10);
        Triad<Apfloat> hex0_norm  = VectorMath.normalHex(vC0, vC6,  vC22, vC4,  vC26, vC10, true);
        Hexad<Tuple<Apfloat>> hex1  = new Hexad<>(vC0, vC10, vC18, vC3,  vC20, vC12);
        Triad<Apfloat> hex1_norm  = VectorMath.normalHex(vC0, vC10, vC18, vC3,  vC20, vC12, true);
        Hexad<Tuple<Apfloat>> hex2  = new Hexad<>(vC0, vC12, vC28, vC5,  vC24, vC8);
        Triad<Apfloat> hex2_norm  = VectorMath.normalHex(vC0, vC12, vC28, vC5,  vC24, vC8, true);
        Hexad<Tuple<Apfloat>> hex3  = new Hexad<>(vC0, vC8,  vC16, vC2,  vC14, vC6);
        Triad<Apfloat> hex3_norm  = VectorMath.normalHex(vC0, vC8,  vC16, vC2,  vC14, vC6, true );
        Hexad<Tuple<Apfloat>> hex4  = new Hexad<>(vC1, vC7,  vC15, vC2,  vC17, vC9);
        Triad<Apfloat> hex4_norm  = VectorMath.normalHex(vC1, vC7,  vC15, vC2,  vC17, vC9, true);
        Hexad<Tuple<Apfloat>> hex5  = new Hexad<>(vC1, vC9,  vC25, vC5,  vC29, vC13);
        Triad<Apfloat> hex5_norm  = VectorMath.normalHex(vC1, vC9,  vC25, vC5,  vC29, vC13, true);
        Hexad<Tuple<Apfloat>> hex6  = new Hexad<>(vC1, vC13, vC21, vC3,  vC19, vC11);
        Triad<Apfloat> hex6_norm  = VectorMath.normalHex(vC1, vC13, vC21, vC3,  vC19, vC11, true);
        Hexad<Tuple<Apfloat>> hex7  = new Hexad<>(vC1, vC11, vC27, vC4,  vC23, vC7);
        Triad<Apfloat> hex7_norm  = VectorMath.normalHex(vC1, vC11, vC27, vC4,  vC23, vC7, true);
        Hexad<Tuple<Apfloat>> hex8  = new Hexad<>(vC2, vC15, vC23, vC4,  vC22, vC14);
        Triad<Apfloat> hex8_norm  = VectorMath.normalHex(vC2, vC15, vC23, vC4,  vC22, vC14, true);
        Hexad<Tuple<Apfloat>> hex9  = new Hexad<>(vC2, vC16, vC24, vC5,  vC25, vC17);
        Triad<Apfloat> hex9_norm  = VectorMath.normalHex(vC2, vC16, vC24, vC5,  vC25, vC17, true);
        Hexad<Tuple<Apfloat>> hex10 = new Hexad<>(vC3, vC18, vC26, vC4,  vC27, vC19);
        Triad<Apfloat> hex10_norm = VectorMath.normalHex(vC3, vC18, vC26, vC4,  vC27, vC19, true);
        Hexad<Tuple<Apfloat>> hex11 = new Hexad<>(vC3, vC21, vC29, vC5,  vC28, vC20);
        Triad<Apfloat> hex11_norm = VectorMath.normalHex(vC3, vC21, vC29, vC5,  vC28, vC20, true);
        faces_hex = new Dodecad<>(
                hex0, hex1, hex2, hex3,
                hex4, hex5, hex6, hex7,
                hex8, hex9, hex10, hex11
        );
        face_norms_hex = new Dodecad<>(
                hex0_norm, hex1_norm, hex2_norm, hex3_norm,
                hex4_norm, hex5_norm, hex6_norm, hex7_norm,
                hex8_norm, hex9_norm, hex10_norm, hex11_norm
        );

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0  = new Triad<>(vC6,  vC14, vC22);
        Triad<Apfloat> tri0_norm  = VectorMath.normalTriple(vC6,  vC14, vC22, true);
        Triad<Tuple<Apfloat>> tri1  = new Triad<>(vC7,  vC23, vC15);
        Triad<Apfloat> tri1_norm  = VectorMath.normalTriple(vC7,  vC23, vC15, true);
        Triad<Tuple<Apfloat>> tri2  = new Triad<>(vC8,  vC24, vC16);
        Triad<Apfloat> tri2_norm  = VectorMath.normalTriple(vC8,  vC24, vC16, true);
        Triad<Tuple<Apfloat>> tri3  = new Triad<>(vC9,  vC17, vC25);
        Triad<Apfloat> tri3_norm  = VectorMath.normalTriple(vC9,  vC17, vC25, true);
        Triad<Tuple<Apfloat>> tri4  = new Triad<>(vC10, vC26, vC18);
        Triad<Apfloat> tri4_norm  = VectorMath.normalTriple(vC10, vC26, vC18, true);
        Triad<Tuple<Apfloat>> tri5  = new Triad<>(vC11, vC19, vC27);
        Triad<Apfloat> tri5_norm  = VectorMath.normalTriple(vC11, vC19, vC27, true);
        Triad<Tuple<Apfloat>> tri6  = new Triad<>(vC12, vC20, vC28);
        Triad<Apfloat> tri6_norm  = VectorMath.normalTriple(vC12, vC20, vC28, true);
        Triad<Tuple<Apfloat>> tri7  = new Triad<>(vC13, vC29, vC21);
        Triad<Apfloat> tri7_norm  = VectorMath.normalTriple(vC13, vC29, vC21, true);
        faces_tri = new Octad<>(
                tri0, tri1, tri2, tri3,
                tri4, tri5, tri6, tri7
        );
        face_norms_tri = new Octad<>(
                tri0_norm, tri1_norm, tri2_norm, tri3_norm,
                tri4_norm, tri5_norm, tri6_norm, tri7_norm

        );
    }

    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < faces_hex.fetchSize(); i++){

            // === Check triangular faces ===
            if (i < face_norms_tri.fetchSize()) {
                Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.fetch(i);
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
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.fetch(i);
                Apfloat d = VectorMath.dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }


            // === Check hexagonal faces ===
            // check if each point lies within bounds of each face

            // get vertices of vertA
            Hexad<Tuple<Apfloat>> face = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
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
            Triad<Apfloat> face_norm = (Triad<Apfloat>) face_norms_hex.fetch(i);
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
