package Shapes;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.*;
import Lattice.*;
import Atom.*;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/Tetrahedron.txt...</a>}
 */
public class Tetrahedron extends Shape{

    // 4 triangular faces
    Tetrad<Tuple<Tuple<Apfloat>>> faces_tri;
    Tetrad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat constant
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat SQRT2 = ApfloatMath.sqrt(new Apfloat("2", super.precision));
    private final Apfloat C0 = SQRT2.divide(N4);

    public Tetrahedron(
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
        Triad<Apfloat> vB0 = new Triad<>(C0, C0.multiply(NEG_N1), C0);
        Triad<Apfloat> vB1 = new Triad<>(C0, C0, C0.multiply(NEG_N1));
        Triad<Apfloat> vB2 = new Triad<>(C0.multiply(NEG_N1), C0, C0);
        Triad<Apfloat> vB3 = new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), C0.multiply(NEG_N1));

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = VectorMath.mult(VectorMath.normalize(vB0),  super.getRadius().toString());
        Triad<Apfloat> vC1  = VectorMath.mult(VectorMath.normalize(vB1),  super.getRadius().toString());
        Triad<Apfloat> vC2  = VectorMath.mult(VectorMath.normalize(vB2),  super.getRadius().toString());
        Triad<Apfloat> vC3  = VectorMath.mult(VectorMath.normalize(vB3),  super.getRadius().toString());

        // ==== TRIANGULAR FACES VERTICES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC0, vC1, vC2);
        Triad<Apfloat> tri0_norm = VectorMath.normalTriple(vC0, vC1, vC2, true);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC1, vC0, vC3);
        Triad<Apfloat> tri1_norm = VectorMath.normalTriple(vC1, vC0, vC3, true);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC2, vC3, vC0);
        Triad<Apfloat> tri2_norm = VectorMath.normalTriple(vC2, vC3, vC0, true);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC3, vC2, vC1);
        Triad<Apfloat> tri3_norm = VectorMath.normalTriple(vC3, vC2, vC1, true);
        faces_tri = new Tetrad<>(tri0,tri1,tri2,tri3);
        face_norms_tri = new Tetrad<>(tri0_norm,tri1_norm,tri2_norm,tri3_norm);
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 4; i++) {
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
