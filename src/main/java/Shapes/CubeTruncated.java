package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/TruncatedCube.txt...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class CubeTruncated extends Shape {

    // 6 octagonal faces
    Hexad<Tuple<Tuple<Apfloat>>> faces_oct;
    Hexad<Tuple<Apfloat>> face_norms_oct;

    // 8 Triangular  faces
    Octad<Tuple<Tuple<Apfloat>>> faces_tri;
    Octad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat SQRT2 = ApfloatMath.sqrt(N2);
    private final Apfloat FRAC_N1_over_N2 = N1.divide(N2);
    private final Apfloat C0 = (N1.add(SQRT2)).divide(N2);

    public CubeTruncated(
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
        Triad<Apfloat> vB0  = new Triad<>(C0, FRAC_N1_over_N2, C0);
        Triad<Apfloat> vB1  = new Triad<>(C0, FRAC_N1_over_N2, C0.multiply(NEG_N1));
        Triad<Apfloat> vB2  = new Triad<>(C0, FRAC_N1_over_N2.multiply(NEG_N1), C0);
        Triad<Apfloat> vB3  = new Triad<>(C0, FRAC_N1_over_N2.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB4  = new Triad<>(C0.multiply(NEG_N1), FRAC_N1_over_N2, C0);
        Triad<Apfloat> vB5  = new Triad<>(C0.multiply(NEG_N1), FRAC_N1_over_N2, C0.multiply(NEG_N1));
        Triad<Apfloat> vB6  = new Triad<>(C0.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1), C0);
        Triad<Apfloat> vB7  = new Triad<>(C0.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB8  = new Triad<>(C0, C0, FRAC_N1_over_N2);
        Triad<Apfloat> vB9  = new Triad<>(C0, C0, FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(C0, C0.multiply(NEG_N1), FRAC_N1_over_N2);
        Triad<Apfloat> vB11 = new Triad<>(C0, C0.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(C0.multiply(NEG_N1), C0, FRAC_N1_over_N2);
        Triad<Apfloat> vB13 = new Triad<>(C0.multiply(NEG_N1), C0, FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), FRAC_N1_over_N2);
        Triad<Apfloat> vB15 = new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), FRAC_N1_over_N2.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>(FRAC_N1_over_N2, C0, C0);
        Triad<Apfloat> vB17 = new Triad<>(FRAC_N1_over_N2, C0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>(FRAC_N1_over_N2, C0.multiply(NEG_N1), C0);
        Triad<Apfloat> vB19 = new Triad<>(FRAC_N1_over_N2, C0.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), C0, C0);
        Triad<Apfloat> vB21 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), C0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB22 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), C0.multiply(NEG_N1), C0);
        Triad<Apfloat> vB23 = new Triad<>(FRAC_N1_over_N2.multiply(NEG_N1), C0.multiply(NEG_N1), C0.multiply(NEG_N1));

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
        Octad<Tuple<Apfloat>> oct0 = new Octad<>(vC0, vC2, vC10, vC11, vC3, vC1, vC9, vC8);
        Triad<Apfloat> oct0_norm = VectorMath.normalOct(vC0, vC2, vC10, vC11, vC3, vC1, vC9, vC8, true);
        Octad<Tuple<Apfloat>> oct1 = new Octad<>(vC0, vC16, vC20, vC4, vC6, vC22, vC18, vC2);
        Triad<Apfloat> oct1_norm = VectorMath.normalOct(vC0, vC16, vC20, vC4, vC6, vC22, vC18, vC2, true);
        Octad<Tuple<Apfloat>> oct2 = new Octad<>(vC12, vC13, vC5, vC7, vC15, vC14, vC6, vC4);
        Triad<Apfloat> oct2_norm = VectorMath.normalOct(vC12, vC13, vC5, vC7, vC15, vC14, vC6, vC4, true);
        Octad<Tuple<Apfloat>> oct3 = new Octad<>(vC12, vC20, vC16, vC8, vC9, vC17, vC21, vC13);
        Triad<Apfloat> oct3_norm = VectorMath.normalOct(vC12, vC20, vC16, vC8, vC9, vC17, vC21, vC13, true);
        Octad<Tuple<Apfloat>> oct4 = new Octad<>(vC19, vC23, vC7, vC5, vC21, vC17, vC1, vC3);
        Triad<Apfloat> oct4_norm = VectorMath.normalOct(vC19, vC23, vC7, vC5, vC21, vC17, vC1, vC3, true);
        Octad<Tuple<Apfloat>> oct5 = new Octad<>(vC19, vC11, vC10, vC18, vC22, vC14, vC15, vC23);
        Triad<Apfloat> oct5_norm = VectorMath.normalOct(vC19, vC11, vC10, vC18, vC22, vC14, vC15, vC23, true);
        faces_oct = new Hexad<>(oct0, oct1, oct2, oct3, oct4, oct5);
        face_norms_oct = new Hexad<>(oct0_norm, oct1_norm, oct2_norm, oct3_norm, oct4_norm, oct5_norm);

        // ==== Triangular FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC0, vC8, vC16);
        Triad<Apfloat> tri0_norm = VectorMath.normalTriple(vC0, vC8, vC16, true);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC1, vC17, vC9);
        Triad<Apfloat> tri1_norm = VectorMath.normalTriple(vC1, vC17, vC9, true);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC2, vC18, vC10);
        Triad<Apfloat> tri2_norm = VectorMath.normalTriple(vC2, vC18, vC10, true);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC3, vC11, vC19);
        Triad<Apfloat> tri3_norm = VectorMath.normalTriple(vC3, vC11, vC19, true);
        Triad<Tuple<Apfloat>> tri4 = new Triad<>(vC4, vC20, vC12);
        Triad<Apfloat> tri4_norm = VectorMath.normalTriple(vC4, vC20, vC12, true);
        Triad<Tuple<Apfloat>> tri5 = new Triad<>(vC5, vC13, vC21);
        Triad<Apfloat> tri5_norm = VectorMath.normalTriple(vC5, vC13, vC21, true);
        Triad<Tuple<Apfloat>> tri6 = new Triad<>(vC6, vC14, vC22);
        Triad<Apfloat> tri6_norm = VectorMath.normalTriple(vC6, vC14, vC22, true);
        Triad<Tuple<Apfloat>> tri7 = new Triad<>(vC7, vC23, vC15);
        Triad<Apfloat> tri7_norm = VectorMath.normalTriple(vC7, vC23, vC15, true);
        faces_tri = new Octad<>(tri0, tri1, tri2, tri3, tri4, tri5, tri6, tri7);
        face_norms_tri = new Octad<>(tri0_norm, tri1_norm, tri2_norm, tri3_norm, tri4_norm, tri5_norm, tri6_norm, tri7_norm);
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 8; i++) {
            if (i < 6){
                // Octagonal faces
                Octad<Tuple<Apfloat>> face = (Octad<Tuple<Apfloat>>) faces_oct.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_oct.fetch(i);

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
            // Triangular  faces
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
        return true;
    }
}
