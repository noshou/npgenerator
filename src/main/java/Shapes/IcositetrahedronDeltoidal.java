package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.Polyad;
import com.oson.tuple.Tetrad;
import com.oson.tuple.Triad;
import com.oson.tuple.Tuple;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;

/*
* @see https://dmccooey.com/polyhedra/DeltoidalIcositetrahedron.txt
*/
public class IcositetrahedronDeltoidal extends Shape{

    // kite faces
    private final ArrayList<Tetrad<Tuple<Apfloat>>> faces_kte = new ArrayList<>();
    private final ArrayList<Triad<Apfloat>> face_norms_kte = new ArrayList<>();

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N7 = new Apfloat("7", super.precision);
    private final Apfloat SQRT2 = ApfloatMath.sqrt(new Apfloat("2", super.precision)    );

    // C0 = (4 + sqrt(2)) / 7
    final Apfloat C0 = (N4.add(SQRT2)).divide(N7);

    // NEG_C0 = -((4 + sqrt(2)) / 7)
    final Apfloat NEG_C0 = C0.multiply(NEG_N1);

    // C1 = sqrt(2)
    final Apfloat C1 = SQRT2;

    // NEG_C1 = -sqrt(2)
    final Apfloat NEG_C1 = C1.multiply(NEG_N1);

    public IcositetrahedronDeltoidal(
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
        Triad<Apfloat> vB0 = new Triad<>(N0, N0, C1);
        Triad<Apfloat> vB1 = new Triad<>(N0, N0, NEG_C1);
        Triad<Apfloat> vB2 = new Triad<>(C1, N0, N0);
        Triad<Apfloat> vB3 = new Triad<>(NEG_C1, N0, N0);
        Triad<Apfloat> vB4 = new Triad<>(N0, C1, N0);
        Triad<Apfloat> vB5 = new Triad<>(N0, NEG_C1, N0);
        Triad<Apfloat> vB6 = new Triad<>(N1, N0, N1);
        Triad<Apfloat> vB7 = new Triad<>(N1, N0, NEG_N1);
        Triad<Apfloat> vB8 = new Triad<>(NEG_N1, N0, N1);
        Triad<Apfloat> vB9 = new Triad<>(NEG_N1, N0, NEG_N1);
        Triad<Apfloat> vB10 = new Triad<>(N1, N1, N0);
        Triad<Apfloat> vB11 = new Triad<>(N1, NEG_N1, N0);
        Triad<Apfloat> vB12 = new Triad<>(NEG_N1, N1, N0);
        Triad<Apfloat> vB13 = new Triad<>(NEG_N1, NEG_N1, N0);
        Triad<Apfloat> vB14 = new Triad<>(N0, N1, N1);
        Triad<Apfloat> vB15 = new Triad<>(N0, N1, NEG_N1);
        Triad<Apfloat> vB16 = new Triad<>(N0, NEG_N1, N1);
        Triad<Apfloat> vB17 = new Triad<>(N0, NEG_N1, NEG_N1);
        Triad<Apfloat> vB18 = new Triad<>(C0, C0, C0);
        Triad<Apfloat> vB19 = new Triad<>(C0, C0, NEG_C0);
        Triad<Apfloat> vB20 = new Triad<>(C0, NEG_C0, C0);
        Triad<Apfloat> vB21 = new Triad<>(C0, NEG_C0, NEG_C0);
        Triad<Apfloat> vB22 = new Triad<>(NEG_C0, C0, C0);
        Triad<Apfloat> vB23 = new Triad<>(NEG_C0, C0, NEG_C0);
        Triad<Apfloat> vB24 = new Triad<>(NEG_C0, NEG_C0, C0);
        Triad<Apfloat> vB25 = new Triad<>(NEG_C0, NEG_C0, NEG_C0);
        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = VectorMath.mult(VectorMath.normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1  = VectorMath.mult(VectorMath.normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2  = VectorMath.mult(VectorMath.normalize(vB2), super.getRadius().toString());
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

        // ==== KITE FACES ====
        Tetrad<Tuple<Apfloat>> kte0 = new Tetrad<>(vC0, vC6, vC18, vC14);
        Triad<Apfloat> kte0_norm = VectorMath.normalQuad(vC0, vC6, vC18, vC14,true);
        faces_kte.add(kte0);
        face_norms_kte.add(kte0_norm);
        Tetrad<Tuple<Apfloat>> kte1 = new Tetrad<>(vC0, vC14, vC22, vC8);
        Triad<Apfloat> kte1_norm = VectorMath.normalQuad(vC0, vC14, vC22, vC8, true);
        faces_kte.add(kte1);
        face_norms_kte.add(kte1_norm);
        Tetrad<Tuple<Apfloat>> kte2 = new Tetrad<>(vC0, vC8, vC24, vC16);
        Triad<Apfloat> kte2_norm = VectorMath.normalQuad(vC0, vC8, vC24, vC16, true);
        faces_kte.add(kte2);
        face_norms_kte.add(kte2_norm);
        Tetrad<Tuple<Apfloat>> kte3 = new Tetrad<>(vC0, vC16, vC20, vC6);
        Triad<Apfloat> kte3_norm = VectorMath.normalQuad(vC0, vC16, vC20, vC6, true);
        faces_kte.add(kte3);
        face_norms_kte.add(kte3_norm);
        Tetrad<Tuple<Apfloat>> kte4 = new Tetrad<>(vC1, vC9, vC23, vC15);
        Triad<Apfloat> kte4_norm = VectorMath.normalQuad(vC1, vC9, vC23, vC15, true);
        faces_kte.add(kte4);
        face_norms_kte.add(kte4_norm);
        Tetrad<Tuple<Apfloat>> kte5 = new Tetrad<>(vC1, vC15, vC19, vC7);
        Triad<Apfloat> kte5_norm = VectorMath.normalQuad(vC1, vC15, vC19, vC7, true);
        faces_kte.add(kte5);
        face_norms_kte.add(kte5_norm);
        Tetrad<Tuple<Apfloat>> kte6 = new Tetrad<>(vC1, vC7, vC21, vC17);
        Triad<Apfloat> kte6_norm = VectorMath.normalQuad(vC1, vC7, vC21, vC17, true);
        faces_kte.add(kte6);
        face_norms_kte.add(kte6_norm);
        Tetrad<Tuple<Apfloat>> kte7 = new Tetrad<>(vC1, vC17, vC25, vC9);
        Triad<Apfloat> kte7_norm = VectorMath.normalQuad(vC1, vC17, vC25, vC9, true);
        faces_kte.add(kte7);
        face_norms_kte.add(kte7_norm);
        Tetrad<Tuple<Apfloat>> kte8 = new Tetrad<>(vC2, vC7, vC19, vC10);
        Triad<Apfloat> kte8_norm = VectorMath.normalQuad(vC2, vC7, vC19, vC10, true);
        faces_kte.add(kte8);
        face_norms_kte.add(kte8_norm);
        Tetrad<Tuple<Apfloat>> kte9 = new Tetrad<>(vC2, vC10, vC18, vC6);
        Triad<Apfloat> kte9_norm = VectorMath.normalQuad(vC2, vC10, vC18, vC6, true);
        faces_kte.add(kte9);
        face_norms_kte.add(kte9_norm);
        Tetrad<Tuple<Apfloat>> kte10 = new Tetrad<>(vC2, vC6, vC20, vC11);
        Triad<Apfloat> kte10_norm = VectorMath.normalQuad(vC2, vC6, vC20, vC11, true);
        faces_kte.add(kte10);
        face_norms_kte.add(kte10_norm);
        Tetrad<Tuple<Apfloat>> kte11 = new Tetrad<>(vC2, vC11, vC21, vC7);
        Triad<Apfloat> kte11_norm = VectorMath.normalQuad(vC2, vC11, vC21, vC7, true);
        faces_kte.add(kte11);
        face_norms_kte.add(kte11_norm);
        Tetrad<Tuple<Apfloat>> kte12 = new Tetrad<>(vC3, vC8, vC22, vC12);
        Triad<Apfloat> kte12_norm = VectorMath.normalQuad(vC3, vC8, vC22, vC12, true);
        faces_kte.add(kte12);
        face_norms_kte.add(kte12_norm);
        Tetrad<Tuple<Apfloat>> kte13 = new Tetrad<>(vC3, vC12, vC23, vC9);
        Triad<Apfloat> kte13_norm = VectorMath.normalQuad(vC3, vC12, vC23, vC9, true);
        faces_kte.add(kte13);
        face_norms_kte.add(kte13_norm);
        Tetrad<Tuple<Apfloat>> kte14 = new Tetrad<>(vC3, vC9, vC25, vC13);
        Triad<Apfloat> kte14_norm = VectorMath.normalQuad(vC3, vC9, vC25, vC13, true);
        faces_kte.add(kte14);
        face_norms_kte.add(kte14_norm);
        Tetrad<Tuple<Apfloat>> kte15 = new Tetrad<>(vC3, vC13, vC24, vC8);
        Triad<Apfloat> kte15_norm = VectorMath.normalQuad(vC3, vC13, vC24, vC8, true);
        faces_kte.add(kte15);
        face_norms_kte.add(kte15_norm);
        Tetrad<Tuple<Apfloat>> kte16 = new Tetrad<>(vC4, vC10, vC19, vC15);
        Triad<Apfloat> kte16_norm = VectorMath.normalQuad(vC4, vC10, vC19, vC15, true);
        faces_kte.add(kte16);
        face_norms_kte.add(kte16_norm);
        Tetrad<Tuple<Apfloat>> kte17 = new Tetrad<>(vC4, vC15, vC23, vC12);
        Triad<Apfloat> kte17_norm = VectorMath.normalQuad(vC4, vC15, vC23, vC12, true);
        faces_kte.add(kte17);
        face_norms_kte.add(kte17_norm);
        Tetrad<Tuple<Apfloat>> kte18 = new Tetrad<>(vC4, vC12, vC22, vC14);
        Triad<Apfloat> kte18_norm = VectorMath.normalQuad(vC4, vC12, vC22, vC14, true);
        faces_kte.add(kte18);
        face_norms_kte.add(kte18_norm);
        Tetrad<Tuple<Apfloat>> kte19 = new Tetrad<>(vC4, vC14, vC18, vC10);
        Triad<Apfloat> kte19_norm = VectorMath.normalQuad(vC4, vC14, vC18, vC10, true);
        faces_kte.add(kte19);
        face_norms_kte.add(kte19_norm);
        Tetrad<Tuple<Apfloat>> kte20 = new Tetrad<>(vC5, vC11, vC20, vC16);
        Triad<Apfloat> kte20_norm = VectorMath.normalQuad(vC5, vC11, vC20, vC16, true);
        faces_kte.add(kte20);
        face_norms_kte.add(kte20_norm);
        Tetrad<Tuple<Apfloat>> kte21 = new Tetrad<>(vC5, vC16, vC24, vC13);
        Triad<Apfloat> kte21_norm = VectorMath.normalQuad(vC5, vC16, vC24, vC13, true);
        faces_kte.add(kte21);
        face_norms_kte.add(kte21_norm);
        Tetrad<Tuple<Apfloat>> kte22 = new Tetrad<>(vC5, vC13, vC25, vC17);
        Triad<Apfloat> kte22_norm = VectorMath.normalQuad(vC5, vC13, vC25, vC17, true);
        faces_kte.add(kte22);
        face_norms_kte.add(kte22_norm);
        Tetrad<Tuple<Apfloat>> kte23 = new Tetrad<>(vC5, vC17, vC21, vC11);
        Triad<Apfloat> kte23_norm = VectorMath.normalQuad(vC5, vC17, vC21, vC11, true);
        faces_kte.add(kte23);
        face_norms_kte.add(kte23_norm);
    }

    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_kte.size(); i++){
            // === Check square faces ===
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_kte.get(i);
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
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_kte.get(i);
            Apfloat d = VectorMath.dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }

        return true;
    }
}
