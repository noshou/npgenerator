package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/TruncatedTetrahedron.txt...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class TetrahedronTruncated extends Shape {

    // 4 Hexagonal faces
    Tetrad<Tuple<Tuple<Apfloat>>> faces_hex;
    Tetrad<Tuple<Apfloat>> face_norms_hex;

    // 4 triangular faces
    Tetrad <Tuple<Tuple<Apfloat>>> faces_tri;
    Tetrad <Tuple<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1    = new Apfloat("-1", super.precision);
    private final Apfloat N0   = new Apfloat("0", super.precision);
    private final Apfloat N3  = new Apfloat("3", super.precision);
    private final Apfloat N4   = new Apfloat("4", super.precision);
    private final Apfloat SQRT2   = ApfloatMath.sqrt(new Apfloat("2", super.precision));

    private final Apfloat C0 = SQRT2.divide(N4);
    private final Apfloat C1 = N3.multiply(SQRT2).divide(N4);

    public TetrahedronTruncated(
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
        Triad<Apfloat> vB0 = new Triad<>(C0, C0.multiply(NEG_N1), C1);
        Triad<Apfloat> vB1 = new Triad<>(C0, C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB2 = new Triad<>(C0.multiply(NEG_N1), C0, C1);
        Triad<Apfloat> vB3 = new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB4 = new Triad<>(C1, C0.multiply(NEG_N1), C0);
        Triad<Apfloat> vB5 = new Triad<>(C1, C0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB6 = new Triad<>(C1.multiply(NEG_N1), C0, C0);
        Triad<Apfloat> vB7 = new Triad<>(C1.multiply(NEG_N1), C0.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB8 = new Triad<>(C0, C1.multiply(NEG_N1), C0);
        Triad<Apfloat> vB9 = new Triad<>(C0, C1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(C0.multiply(NEG_N1), C1, C0);
        Triad<Apfloat> vB11 = new Triad<>(C0.multiply(NEG_N1), C1.multiply(NEG_N1), C0.multiply(NEG_N1));


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

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC0, vC4, vC5, vC9, vC10, vC2);
        Triad<Apfloat> hex0_norm = VectorMath.normalHex(vC0, vC4, vC5, vC9, vC10, vC2, true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC1, vC5, vC4, vC8, vC11, vC3);
        Triad<Apfloat> hex1_norm = VectorMath.normalHex(vC1, vC5, vC4, vC8, vC11, vC3, true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC2, vC6, vC7, vC11, vC8, vC0);
        Triad<Apfloat> hex2_norm = VectorMath.normalHex(vC2, vC6, vC7, vC11, vC8, vC0, true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC3, vC7, vC6, vC10, vC9, vC1);
        Triad<Apfloat> hex3_norm = VectorMath.normalHex(vC3, vC7, vC6, vC10, vC9, vC1, true);
        faces_hex = new Tetrad<>(hex0, hex1, hex2, hex3);
        face_norms_hex = new Tetrad<>(hex0_norm, hex1_norm, hex2_norm, hex3_norm);

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC0, vC8, vC4);
        Triad<Apfloat> tri0_norm = VectorMath.normalTriple(vC0, vC8, vC4, true);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC1, vC9, vC5);
        Triad<Apfloat> tri1_norm = VectorMath.normalTriple(vC1, vC9, vC5, true);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC2, vC10, vC6);
        Triad<Apfloat> tri2_norm = VectorMath.normalTriple(vC2, vC10, vC6, true);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC3, vC11, vC7);
        Triad<Apfloat> tri3_norm = VectorMath.normalTriple(vC3, vC11, vC7, true);
        faces_tri = new Tetrad<>(tri0, tri1, tri2, tri3);
        face_norms_tri = new Tetrad<>(tri0_norm, tri1_norm, tri2_norm, tri3_norm);
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 4; i++) {
            // hexagonal faces
            Hexad<Tuple<Apfloat>> face_h = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
            Triad<Apfloat> vertA_h = (Triad<Apfloat>) face_h.fetch(0);
            Triad<Apfloat> norm_h = (Triad<Apfloat>) face_norms_hex.fetch(i);

            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
            Triad<Apfloat> m_h = VectorMath.subs(point_cart, vertA_h);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Apfloat d_h = VectorMath.dot_prod(norm_h, m_h);

            // point outside if dot product > 0
            if (d_h.compareTo(N0) > 0) {
                return false;
            }

            // triangular faces
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
