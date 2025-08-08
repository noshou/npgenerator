package Shapes;

import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;
import Lattice.*;
import Atom.*;
import Utilities.VectorMath;

/*
 * https://dmccooey.com/polyhedra/Dodecahedron.txt
 */
@SuppressWarnings("FieldCanBeLocal")
public class Dodecahedron extends Shape {

    // 12 pentagonal faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces;
    private final Dodecad<Tuple<Apfloat>> face_norms;

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(new Apfloat("5", super.precision));
    private final Apfloat FRAC_1_over_2 = N1.divide(N2);

    // C0 = (1 + sqrt(5)) / 4
    private final Apfloat C0 = (N1.add(SQRT5)).divide(N4);

    // C1 = (3 + sqrt(5)) / 4
    private final Apfloat C1 = (N3.add(SQRT5)).divide(N4);


    public Dodecahedron(
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

        Triad<Apfloat> vB0 = new Triad<>(N0, FRAC_1_over_2, C1);
        Triad<Apfloat> vB1 = new Triad<>(N0, FRAC_1_over_2, C1.multiply(NEG_N1));
        Triad<Apfloat> vB2 = new Triad<>(N0, FRAC_1_over_2.multiply(NEG_N1), C1);
        Triad<Apfloat> vB3 = new Triad<>(N0, FRAC_1_over_2.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB4 = new Triad<>(C1, N0, FRAC_1_over_2);
        Triad<Apfloat> vB5 = new Triad<>(C1, N0, FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB6 = new Triad<>(C1.multiply(NEG_N1), N0, FRAC_1_over_2);
        Triad<Apfloat> vB7 = new Triad<>(C1.multiply(NEG_N1), N0, FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB8 = new Triad<>(FRAC_1_over_2, C1, N0);
        Triad<Apfloat> vB9 = new Triad<>(FRAC_1_over_2, C1.multiply(NEG_N1), N0);
        Triad<Apfloat> vB10 = new Triad<>(FRAC_1_over_2.multiply(NEG_N1), C1, N0);
        Triad<Apfloat> vB11 = new Triad<>(FRAC_1_over_2.multiply(NEG_N1), C1.multiply(NEG_N1), N0);
        Triad<Apfloat> vB12 = new Triad<>(C0, C0, C0);
        Triad<Apfloat> vB13 = new Triad<>(C0, C0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(C0, C0.multiply(NEG_N1), C0);
        Triad<Apfloat> vB15 = new Triad<>(C0, C0.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>(C0.multiply(NEG_N1), C0, C0);
        Triad<Apfloat> vB17 = new Triad<>(C0.multiply(NEG_N1), C0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), C0);
        Triad<Apfloat> vB19 = new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), C0.multiply(NEG_N1));

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

        // ==== PENTAGONAL FACES ====
        Pentad<Tuple<Apfloat>> pnt0 = new Pentad<>(vC0, vC2, vC14, vC4, vC12);
        Triad<Apfloat> pnt0_norm = VectorMath.normalPent(vC0, vC2, vC14, vC4, vC12, true);
        Pentad<Tuple<Apfloat>> pnt1 = new Pentad<>(vC0, vC12, vC8, vC10, vC16);
        Triad<Apfloat> pnt1_norm = VectorMath.normalPent(vC0, vC12, vC8, vC10, vC16, true);
        Pentad<Tuple<Apfloat>> pnt2 = new Pentad<>(vC0, vC16, vC6, vC18, vC2);
        Triad<Apfloat> pnt2_norm = VectorMath.normalPent(vC0, vC16, vC6, vC18, vC2, true);
        Pentad<Tuple<Apfloat>> pnt3 = new Pentad<>(vC7, vC6, vC16, vC10, vC17);
        Triad<Apfloat> pnt3_norm = VectorMath.normalPent(vC7, vC6, vC16, vC10, vC17, true);
        Pentad<Tuple<Apfloat>> pnt4 = new Pentad<>(vC7, vC17, vC1, vC3, vC19);
        Triad<Apfloat> pnt4_norm = VectorMath.normalPent(vC7, vC17, vC1, vC3, vC19, true);
        Pentad<Tuple<Apfloat>> pnt5 = new Pentad<>(vC7, vC19, vC11, vC18, vC6);
        Triad<Apfloat> pnt5_norm = VectorMath.normalPent(vC7, vC19, vC11, vC18, vC6, true);
        Pentad<Tuple<Apfloat>> pnt6 = new Pentad<>(vC9, vC11, vC19, vC3, vC15);
        Triad<Apfloat> pnt6_norm = VectorMath.normalPent(vC9, vC11, vC19, vC3, vC15, true);
        Pentad<Tuple<Apfloat>> pnt7 = new Pentad<>(vC9, vC15, vC5, vC4, vC14);
        Triad<Apfloat> pnt7_norm = VectorMath.normalPent(vC9, vC15, vC5, vC4, vC14, true);
        Pentad<Tuple<Apfloat>> pnt8 = new Pentad<>(vC9, vC14, vC2, vC18, vC11);
        Triad<Apfloat> pnt8_norm = VectorMath.normalPent(vC9, vC14, vC2, vC18, vC11, true);
        Pentad<Tuple<Apfloat>> pnt9 = new Pentad<>(vC13, vC1, vC17, vC10, vC8);
        Triad<Apfloat> pnt9_norm = VectorMath.normalPent(vC13, vC1, vC17, vC10, vC8, true);
        Pentad<Tuple<Apfloat>> pnt10 = new Pentad<>(vC13, vC8, vC12, vC4, vC5);
        Triad<Apfloat> pnt10_norm = VectorMath.normalPent(vC13, vC8, vC12, vC4, vC5, true);
        Pentad<Tuple<Apfloat>> pnt11 = new Pentad<>(vC13, vC5, vC15, vC3, vC1);
        Triad<Apfloat> pnt11_norm = VectorMath.normalPent(vC13, vC5, vC15, vC3, vC1, true);
        faces = new Dodecad<>(
                pnt0, pnt1, pnt2,
                pnt3, pnt4, pnt5,
                pnt6, pnt7, pnt8,
                pnt9, pnt10, pnt11
        );
        face_norms = new Dodecad<>(
                pnt0_norm, pnt1_norm, pnt2_norm,
                pnt3_norm, pnt4_norm, pnt5_norm,
                pnt6_norm, pnt7_norm, pnt8_norm,
                pnt9_norm, pnt10_norm, pnt11_norm
        );
    }

    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        // check if each point lies within bounds of each face
        for (int i = 0; i < faces.fetchSize(); i++) {

            // get vertices of vertA
            Pentad<Tuple<Apfloat>> face = (Pentad<Tuple<Apfloat>>) faces.fetch(i);
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
            Triad<Apfloat> face_norm = (Triad<Apfloat>) face_norms.fetch(i);
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
