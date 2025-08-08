package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/ChamferedCube2.txt>...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class CubeChamferedBiscribed extends Shape {

    // 6 square faces
    private final Hexad<Tuple<Tuple<Apfloat>>> faces_sqr;
    private final Hexad<Tuple<Apfloat>> face_norms_sqr;

    // 12 hexagonal faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_hex;
    private final Dodecad<Tuple<Apfloat>> face_norms_hex;

    // Apfloat constants
    private final Apfloat NEG_N1    = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat SQRT2 = ApfloatMath.sqrt(N2);

    // C0 = sqrt(2) - 1
    private final Apfloat C0 = SQRT2.subtract(N1);

    // C1 = sqrt(2) / 2
    private final Apfloat C1 = SQRT2.divide(N2);

    public CubeChamferedBiscribed(
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
        Triad<Apfloat> vB0 = new Triad<>(C0, C0, N1);
        Triad<Apfloat> vB1 = new Triad<>(C0, C0, NEG_N1);
        Triad<Apfloat> vB2 = new Triad<>(C0, C0.multiply(NEG_N1), N1);
        Triad<Apfloat> vB3 = new Triad<>(C0, C0.multiply(NEG_N1), NEG_N1);
        Triad<Apfloat> vB4 = new Triad<>(C0.multiply(NEG_N1), C0, N1);
        Triad<Apfloat> vB5 = new Triad<>(C0.multiply(NEG_N1), C0, NEG_N1);
        Triad<Apfloat> vB6 = new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), N1);
        Triad<Apfloat> vB7 = new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), NEG_N1);
        Triad<Apfloat> vB8 = new Triad<>(N1, C0, C0);
        Triad<Apfloat> vB9 = new Triad<>(N1, C0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(N1, C0.multiply(NEG_N1), C0);
        Triad<Apfloat> vB11 = new Triad<>(N1, C0.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(NEG_N1, C0, C0);
        Triad<Apfloat> vB13 = new Triad<>(NEG_N1, C0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(NEG_N1, C0.multiply(NEG_N1), C0);
        Triad<Apfloat> vB15 = new Triad<>(NEG_N1, C0.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>(C0, N1, C0);
        Triad<Apfloat> vB17 = new Triad<>(C0, N1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>(C0, NEG_N1, C0);
        Triad<Apfloat> vB19 = new Triad<>(C0, NEG_N1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>(C0.multiply(NEG_N1), N1, C0);
        Triad<Apfloat> vB21 = new Triad<>(C0.multiply(NEG_N1), N1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB22 = new Triad<>(C0.multiply(NEG_N1), NEG_N1, C0);
        Triad<Apfloat> vB23 = new Triad<>(C0.multiply(NEG_N1), NEG_N1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB24 = new Triad<>(C1, C1, C1);
        Triad<Apfloat> vB25 = new Triad<>(C1, C1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB26 = new Triad<>(C1, C1.multiply(NEG_N1), C1);
        Triad<Apfloat> vB27 = new Triad<>(C1, C1.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB28 = new Triad<>(C1.multiply(NEG_N1), C1, C1);
        Triad<Apfloat> vB29 = new Triad<>(C1.multiply(NEG_N1), C1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB30 = new Triad<>(C1.multiply(NEG_N1), C1.multiply(NEG_N1), C1);
        Triad<Apfloat> vB31 = new Triad<>(C1.multiply(NEG_N1), C1.multiply(NEG_N1), C1.multiply(NEG_N1));

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

        // ==== QUADRILATERAL FACES ====
        Tetrad<Tuple<Apfloat>> sqr0 = new Tetrad<>(vC0, vC4, vC6, vC2);
        Triad<Apfloat> sqr0_norm = VectorMath.normalQuad(vC0, vC4, vC6, vC2, true);
        Tetrad<Tuple<Apfloat>> sqr1 = new Tetrad<>(vC1, vC3, vC7, vC5);
        Triad<Apfloat> sqr1_norm = VectorMath.normalQuad(vC1, vC3, vC7, vC5, true);
        Tetrad<Tuple<Apfloat>> sqr2 = new Tetrad<>(vC8, vC10, vC11, vC9);
        Triad<Apfloat> sqr2_norm = VectorMath.normalQuad(vC8, vC10, vC11, vC9, true);
        Tetrad<Tuple<Apfloat>> sqr3 = new Tetrad<>(vC12, vC13, vC15, vC14);
        Triad<Apfloat> sqr3_norm = VectorMath.normalQuad(vC12, vC13, vC15, vC14, true);
        Tetrad<Tuple<Apfloat>> sqr4 = new Tetrad<>(vC16, vC17, vC21, vC20);
        Triad<Apfloat> sqr4_norm = VectorMath.normalQuad(vC16, vC17, vC21, vC20, true);
        Tetrad<Tuple<Apfloat>> sqr5 = new Tetrad<>(vC18, vC22, vC23, vC19);
        Triad<Apfloat> sqr5_norm = VectorMath.normalQuad(vC18, vC22, vC23, vC19, true);
        faces_sqr = new Hexad<>(sqr0,sqr1,sqr2,sqr3,sqr4,sqr5);
        face_norms_sqr = new Hexad<>(sqr0_norm,sqr1_norm,sqr2_norm,sqr3_norm,sqr4_norm,sqr5_norm);

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC24, vC0, vC2, vC26, vC10, vC8);
        Triad<Apfloat> hex0_norm = VectorMath.normalHex(vC24, vC0, vC2, vC26, vC10, vC8, true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC24, vC8, vC9, vC25, vC17, vC16);
        Triad<Apfloat> hex1_norm = VectorMath.normalHex(vC24, vC8, vC9, vC25, vC17, vC16, true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC24, vC16, vC20, vC28, vC4, vC0);
        Triad<Apfloat> hex2_norm = VectorMath.normalHex(vC24, vC16, vC20, vC28, vC4, vC0, true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC27, vC3, vC1, vC25, vC9, vC11);
        Triad<Apfloat> hex3_norm = VectorMath.normalHex(vC27, vC3, vC1, vC25, vC9, vC11, true);
        Hexad<Tuple<Apfloat>> hex4 = new Hexad<>(vC27, vC11, vC10, vC26, vC18, vC19);
        Triad<Apfloat> hex4_norm = VectorMath.normalHex(vC27, vC11, vC10, vC26, vC18, vC19, true);
        Hexad<Tuple<Apfloat>> hex5 = new Hexad<>(vC27, vC19, vC23, vC31, vC7, vC3);
        Triad<Apfloat> hex5_norm = VectorMath.normalHex(vC27, vC19, vC23, vC31, vC7, vC3, true);
        Hexad<Tuple<Apfloat>> hex6 = new Hexad<>(vC29, vC5, vC7, vC31, vC15, vC13);
        Triad<Apfloat> hex6_norm = VectorMath.normalHex(vC29, vC5, vC7, vC31, vC15, vC13, true);
        Hexad<Tuple<Apfloat>> hex7 = new Hexad<>(vC29, vC13, vC12, vC28, vC20, vC21);
        Triad<Apfloat> hex7_norm = VectorMath.normalHex(vC29, vC13, vC12, vC28, vC20, vC21, true);
        Hexad<Tuple<Apfloat>> hex8 = new Hexad<>(vC29, vC21, vC17, vC25, vC1, vC5);
        Triad<Apfloat> hex8_norm = VectorMath.normalHex(vC29, vC21, vC17, vC25, vC1, vC5, true);
        Hexad<Tuple<Apfloat>> hex9 = new Hexad<>(vC30, vC6, vC4, vC28, vC12, vC14);
        Triad<Apfloat> hex9_norm = VectorMath.normalHex(vC30, vC6, vC4, vC28, vC12, vC14, true);
        Hexad<Tuple<Apfloat>> hex10 = new Hexad<>(vC30, vC14, vC15, vC31, vC23, vC22);
        Triad<Apfloat> hex10_norm = VectorMath.normalHex(vC30, vC14, vC15, vC31, vC23, vC22, true);
        Hexad<Tuple<Apfloat>> hex11 = new Hexad<>(vC30, vC22, vC18, vC26, vC2, vC6);
        Triad<Apfloat> hex11_norm = VectorMath.normalHex(vC30, vC22, vC18, vC26, vC2, vC6, true);
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
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 12; i++) {
            if (i < 6){
                // square faces
                Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_sqr.fetch(i);

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
            // Hexagonal  faces
            Hexad<Tuple<Apfloat>> face = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_hex.fetch(i);

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
