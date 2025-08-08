package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/TruncatedOctahedron.txt...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class OctahedronTruncated extends Shape {

    // 6 square faces
    Hexad<Tuple<Tuple<Apfloat>>> faces_sqr;
    Hexad<Tuple<Apfloat>> face_norms_sqr;

    // 8 Hexagonal  faces
    Octad<Tuple<Tuple<Apfloat>>> faces_hex;
    Octad<Tuple<Apfloat>> face_norms_hex;

    // Apfloat constants
    private final Apfloat NEG_N1    = new Apfloat("-1", super.precision);
    private final Apfloat N0   = new Apfloat("0", super.precision);
    private final Apfloat N2   = new Apfloat("2", super.precision);
    private final Apfloat SQRT2   = ApfloatMath.sqrt(N2);

    private final Apfloat C0 = SQRT2.divide(N2);
    private final Apfloat C1 = SQRT2;

    public OctahedronTruncated(
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
        Triad<Apfloat> vB0  = new Triad<>(C0, N0, C1);
        Triad<Apfloat> vB1  = new Triad<>(C0, N0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB2  = new Triad<>(C0.multiply(NEG_N1), N0, C1);
        Triad<Apfloat> vB3  = new Triad<>(C0.multiply(NEG_N1), N0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB4  = new Triad<>(C1, C0, N0);
        Triad<Apfloat> vB5  = new Triad<>(C1, C0.multiply(NEG_N1), N0);
        Triad<Apfloat> vB6  = new Triad<>(C1.multiply(NEG_N1), C0, N0);
        Triad<Apfloat> vB7  = new Triad<>(C1.multiply(NEG_N1), C0.multiply(NEG_N1), N0);
        Triad<Apfloat> vB8  = new Triad<>(N0, C1, C0);
        Triad<Apfloat> vB9  = new Triad<>(N0, C1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(N0, C1.multiply(NEG_N1), C0);
        Triad<Apfloat> vB11 = new Triad<>(N0, C1.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(N0, C0, C1);
        Triad<Apfloat> vB13 = new Triad<>(N0, C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(N0, C0.multiply(NEG_N1), C1);
        Triad<Apfloat> vB15 = new Triad<>(N0, C0.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>(C1, N0, C0);
        Triad<Apfloat> vB17 = new Triad<>(C1, N0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>(C1.multiply(NEG_N1), N0, C0);
        Triad<Apfloat> vB19 = new Triad<>(C1.multiply(NEG_N1), N0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>(C0, C1, N0);
        Triad<Apfloat> vB21 = new Triad<>(C0, C1.multiply(NEG_N1), N0);
        Triad<Apfloat> vB22 = new Triad<>(C0.multiply(NEG_N1), C1, N0);
        Triad<Apfloat> vB23 = new Triad<>(C0.multiply(NEG_N1), C1.multiply(NEG_N1), N0);

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

        // ==== RECTANGULAR FACES ====
        Tetrad<Tuple<Apfloat>> sqr0 = new Tetrad<>(vC0, vC12, vC2, vC14);
        Triad<Apfloat> sqr0_norm = VectorMath.normalQuad(vC0, vC12, vC2, vC14, true);
        Tetrad<Tuple<Apfloat>> sqr1 = new Tetrad<>(vC1, vC13, vC9, vC20);
        Triad<Apfloat> sqr1_norm = VectorMath.normalQuad(vC1, vC13, vC9, vC20, true);
        Tetrad<Tuple<Apfloat>> sqr2 = new Tetrad<>(vC2, vC12, vC8, vC22);
        Triad<Apfloat> sqr2_norm = VectorMath.normalQuad(vC2, vC12, vC8, vC22, true);
        Tetrad<Tuple<Apfloat>> sqr3 = new Tetrad<>(vC3, vC15, vC11, vC23);
        Triad<Apfloat> sqr3_norm = VectorMath.normalQuad(vC3, vC15, vC11, vC23, true);
        Tetrad<Tuple<Apfloat>> sqr4 = new Tetrad<>(vC4, vC20, vC8, vC12);
        Triad<Apfloat> sqr4_norm = VectorMath.normalQuad(vC4, vC20, vC8, vC12, true);
        Tetrad<Tuple<Apfloat>> sqr5 = new Tetrad<>(vC5, vC21, vC11, vC15);
        Triad<Apfloat> sqr5_norm = VectorMath.normalQuad(vC5, vC21, vC11, vC15, true);
        faces_sqr = new Hexad<>(sqr0, sqr1, sqr2, sqr3, sqr4, sqr5);
        face_norms_sqr = new Hexad<>(sqr0_norm, sqr1_norm, sqr2_norm, sqr3_norm, sqr4_norm, sqr5_norm);

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC0, vC14, vC10, vC21, vC5, vC16);
        Triad<Apfloat> hex0_norm = VectorMath.normalHex(vC0, vC14, vC10, vC21, vC5, vC16, true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC1, vC15, vC3, vC13, vC9, vC22);
        Triad<Apfloat> hex1_norm = VectorMath.normalHex(vC1, vC15, vC3, vC13, vC9, vC22, true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC4, vC16, vC5, vC17, vC7, vC18);
        Triad<Apfloat> hex2_norm = VectorMath.normalHex(vC4, vC16, vC5, vC17, vC7, vC18, true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC6, vC19, vC7, vC18, vC2, vC12);
        Triad<Apfloat> hex3_norm = VectorMath.normalHex(vC6, vC19, vC7, vC18, vC2, vC12, true);
        Hexad<Tuple<Apfloat>> hex4 = new Hexad<>(vC8, vC20, vC9, vC22, vC0, vC14);
        Triad<Apfloat> hex4_norm = VectorMath.normalHex(vC8, vC20, vC9, vC22, vC0, vC14, true);
        Hexad<Tuple<Apfloat>> hex5 = new Hexad<>(vC10, vC23, vC11, vC21, vC1, vC13);
        Triad<Apfloat> hex5_norm = VectorMath.normalHex(vC10, vC23, vC11, vC21, vC1, vC13, true);
        Hexad<Tuple<Apfloat>> hex6 = new Hexad<>(vC6, vC19, vC3, vC15, vC5, vC17);
        Triad<Apfloat> hex6_norm = VectorMath.normalHex(vC6, vC19, vC3, vC15, vC5, vC17, true);
        Hexad<Tuple<Apfloat>> hex7 = new Hexad<>(vC2, vC18, vC7, vC23, vC11, vC21);
        Triad<Apfloat> hex7_norm = VectorMath.normalHex(vC2, vC18, vC7, vC23, vC11, vC21, true);
        faces_hex = new Octad<>(hex0, hex1, hex2, hex3, hex4, hex5, hex6, hex7);
        face_norms_hex = new Octad<>(hex0_norm, hex1_norm, hex2_norm, hex3_norm, hex4_norm, hex5_norm, hex6_norm, hex7_norm);
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 8; i++) {
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
