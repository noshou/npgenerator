package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/TriakisOctahedron.txt...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class OctahedronTriakis extends Shape {

    // 24 isosceles faces
    ArrayList<Triad<Tuple<Apfloat>>> faces_tri;
    ArrayList<Triad<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N1plusSQRT2 = N1.add(ApfloatMath.sqrt(new Apfloat("2", super.precision)));


    public OctahedronTriakis(
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
        Triad<Apfloat> vB0 = new Triad<>(N0, N0, N1plusSQRT2);
        Triad<Apfloat> vB1 = new Triad<>(N0, N0, N1plusSQRT2.multiply(NEG_N1));
        Triad<Apfloat> vB2 = new Triad<>(N1plusSQRT2, N0, N0);
        Triad<Apfloat> vB3 = new Triad<>(N1plusSQRT2.multiply(NEG_N1), N0, N0);
        Triad<Apfloat> vB4 = new Triad<>(N0, N1plusSQRT2, N0);
        Triad<Apfloat> vB5 = new Triad<>(N0, N1plusSQRT2.multiply(NEG_N1), N0);
        Triad<Apfloat> vB6 = new Triad<>(N1, N1, N1);
        Triad<Apfloat> vB7 = new Triad<>(N1, N1, NEG_N1);
        Triad<Apfloat> vB8 = new Triad<>(N1, NEG_N1, N1);
        Triad<Apfloat> vB9 = new Triad<>(N1, NEG_N1, NEG_N1);
        Triad<Apfloat> vB10 = new Triad<>(NEG_N1, N1, N1);
        Triad<Apfloat> vB11 = new Triad<>(NEG_N1, N1, NEG_N1);
        Triad<Apfloat> vB12 = new Triad<>(NEG_N1, NEG_N1, N1);
        Triad<Apfloat> vB13 = new Triad<>(NEG_N1, NEG_N1, NEG_N1);

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

        // ==== TRIANGULAR FACES ====
        faces_tri = new ArrayList<>();
        face_norms_tri = new ArrayList<>();

        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC6, vC0, vC2);
        Triad<Apfloat> tri0_norm = VectorMath.normalTriple(vC6, vC0, vC2, true);
        faces_tri.add(tri0);
        face_norms_tri.add(tri0_norm);

        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC6, vC2, vC4);
        Triad<Apfloat> tri1_norm = VectorMath.normalTriple(vC6, vC2, vC4, true);
        faces_tri.add(tri1);
        face_norms_tri.add(tri1_norm);

        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC6, vC4, vC0);
        Triad<Apfloat> tri2_norm = VectorMath.normalTriple(vC6, vC4, vC0, true);
        faces_tri.add(tri2);
        face_norms_tri.add(tri2_norm);

        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC7, vC1, vC4);
        Triad<Apfloat> tri3_norm = VectorMath.normalTriple(vC7, vC1, vC4, true);
        faces_tri.add(tri3);
        face_norms_tri.add(tri3_norm);

        Triad<Tuple<Apfloat>> tri4 = new Triad<>(vC7, vC4, vC2);
        Triad<Apfloat> tri4_norm = VectorMath.normalTriple(vC7, vC4, vC2, true);
        faces_tri.add(tri4);
        face_norms_tri.add(tri4_norm);

        Triad<Tuple<Apfloat>> tri5 = new Triad<>(vC7, vC2, vC1);
        Triad<Apfloat> tri5_norm = VectorMath.normalTriple(vC7, vC2, vC1, true);
        faces_tri.add(tri5);
        face_norms_tri.add(tri5_norm);

        Triad<Tuple<Apfloat>> tri6 = new Triad<>(vC8, vC0, vC5);
        Triad<Apfloat> tri6_norm = VectorMath.normalTriple(vC8, vC0, vC5, true);
        faces_tri.add(tri6);
        face_norms_tri.add(tri6_norm);

        Triad<Tuple<Apfloat>> tri7 = new Triad<>(vC8, vC5, vC2);
        Triad<Apfloat> tri7_norm = VectorMath.normalTriple(vC8, vC5, vC2, true);
        faces_tri.add(tri7);
        face_norms_tri.add(tri7_norm);

        Triad<Tuple<Apfloat>> tri8 = new Triad<>(vC8, vC2, vC0);
        Triad<Apfloat> tri8_norm = VectorMath.normalTriple(vC8, vC2, vC0, true);
        faces_tri.add(tri8);
        face_norms_tri.add(tri8_norm);

        Triad<Tuple<Apfloat>> tri9 = new Triad<>(vC9, vC1, vC2);
        Triad<Apfloat> tri9_norm = VectorMath.normalTriple(vC9, vC1, vC2, true);
        faces_tri.add(tri9);
        face_norms_tri.add(tri9_norm);

        Triad<Tuple<Apfloat>> tri10 = new Triad<>(vC9, vC2, vC5);
        Triad<Apfloat> tri10_norm = VectorMath.normalTriple(vC9, vC2, vC5, true);
        faces_tri.add(tri10);
        face_norms_tri.add(tri10_norm);

        Triad<Tuple<Apfloat>> tri11 = new Triad<>(vC9, vC5, vC1);
        Triad<Apfloat> tri11_norm = VectorMath.normalTriple(vC9, vC5, vC1, true);
        faces_tri.add(tri11);
        face_norms_tri.add(tri11_norm);

        Triad<Tuple<Apfloat>> tri12 = new Triad<>(vC10, vC0, vC4);
        Triad<Apfloat> tri12_norm = VectorMath.normalTriple(vC10, vC0, vC4, true);
        faces_tri.add(tri12);
        face_norms_tri.add(tri12_norm);

        Triad<Tuple<Apfloat>> tri13 = new Triad<>(vC10, vC4, vC3);
        Triad<Apfloat> tri13_norm = VectorMath.normalTriple(vC10, vC4, vC3, true);
        faces_tri.add(tri13);
        face_norms_tri.add(tri13_norm);

        Triad<Tuple<Apfloat>> tri14 = new Triad<>(vC10, vC3, vC0);
        Triad<Apfloat> tri14_norm = VectorMath.normalTriple(vC10, vC3, vC0, true);
        faces_tri.add(tri14);
        face_norms_tri.add(tri14_norm);

        Triad<Tuple<Apfloat>> tri15 = new Triad<>(vC11, vC1, vC3);
        Triad<Apfloat> tri15_norm = VectorMath.normalTriple(vC11, vC1, vC3, true);
        faces_tri.add(tri15);
        face_norms_tri.add(tri15_norm);

        Triad<Tuple<Apfloat>> tri16 = new Triad<>(vC11, vC3, vC4);
        Triad<Apfloat> tri16_norm = VectorMath.normalTriple(vC11, vC3, vC4, true);
        faces_tri.add(tri16);
        face_norms_tri.add(tri16_norm);

        Triad<Tuple<Apfloat>> tri17 = new Triad<>(vC11, vC4, vC1);
        Triad<Apfloat> tri17_norm = VectorMath.normalTriple(vC11, vC4, vC1, true);
        faces_tri.add(tri17);
        face_norms_tri.add(tri17_norm);

        Triad<Tuple<Apfloat>> tri18 = new Triad<>(vC12, vC0, vC3);
        Triad<Apfloat> tri18_norm = VectorMath.normalTriple(vC12, vC0, vC3, true);
        faces_tri.add(tri18);
        face_norms_tri.add(tri18_norm);

        Triad<Tuple<Apfloat>> tri19 = new Triad<>(vC12, vC3, vC5);
        Triad<Apfloat> tri19_norm = VectorMath.normalTriple(vC12, vC3, vC5, true);
        faces_tri.add(tri19);
        face_norms_tri.add(tri19_norm);

        Triad<Tuple<Apfloat>> tri20 = new Triad<>(vC12, vC5, vC0);
        Triad<Apfloat> tri20_norm = VectorMath.normalTriple(vC12, vC5, vC0, true);
        faces_tri.add(tri20);
        face_norms_tri.add(tri20_norm);

        Triad<Tuple<Apfloat>> tri21 = new Triad<>(vC13, vC1, vC5);
        Triad<Apfloat> tri21_norm = VectorMath.normalTriple(vC13, vC1, vC5, true);
        faces_tri.add(tri21);
        face_norms_tri.add(tri21_norm);

        Triad<Tuple<Apfloat>> tri22 = new Triad<>(vC13, vC5, vC3);
        Triad<Apfloat> tri22_norm = VectorMath.normalTriple(vC13, vC5, vC3, true);
        faces_tri.add(tri22);
        face_norms_tri.add(tri22_norm);

        Triad<Tuple<Apfloat>> tri23 = new Triad<>(vC13, vC3, vC1);
        Triad<Apfloat> tri23_norm = VectorMath.normalTriple(vC13, vC3, vC1, true);
        faces_tri.add(tri23);
        face_norms_tri.add(tri23_norm);
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_tri.size(); i++) {

        // triangular faces
        Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.get(i);
        Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
        Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.get(i);

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
