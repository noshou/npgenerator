package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/TriangularOrthobicupola.txt...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class OrthobicupolaTriangular extends Shape {

    // 6 square faces
    Hexad<Tuple<Tuple<Apfloat>>> faces_sqr;
    Hexad<Tuple<Apfloat>> face_norms_sqr;

    // 8 triangular faces
    Octad<Tuple<Tuple<Apfloat>>> faces_tri;
    Octad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1    = new Apfloat("-1", super.precision);
    private final Apfloat N0   = new Apfloat("0", super.precision);
    private final Apfloat N2   = new Apfloat("2", super.precision);
    private final Apfloat N3  = new Apfloat("3", super.precision);
    private final Apfloat N6  = new Apfloat("6", super.precision);
    private final Apfloat SQRT2   = ApfloatMath.sqrt(N2);

    // golden ratio based constants
    private final Apfloat C0 = SQRT2.divide(N6);
    private final Apfloat C1 = N3.multiply(SQRT2).divide(N2);
    private final Apfloat C2 = N2.multiply(SQRT2).divide(N3);

    public OrthobicupolaTriangular(
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
        Triad<Apfloat> vB0  = new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), C2.multiply(NEG_N1));
        Triad<Apfloat> vB1  = new Triad<>(C0.multiply(NEG_N1), C2.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB2  = new Triad<>(C2.multiply(NEG_N1), C0.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB3  = new Triad<>(C1, N0, C1);
        Triad<Apfloat> vB4  = new Triad<>(C1, N0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB5  = new Triad<>(C1.multiply(NEG_N1), N0, C1);
        Triad<Apfloat> vB6  = new Triad<>(C1, C1, N0);
        Triad<Apfloat> vB7  = new Triad<>(C1, C1.multiply(NEG_N1), N0);
        Triad<Apfloat> vB8  = new Triad<>(C1.multiply(NEG_N1), C1, N0);
        Triad<Apfloat> vB9  = new Triad<>(N0, C1, C1);
        Triad<Apfloat> vB10 = new Triad<>(N0, C1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB11 = new Triad<>(N0, C1.multiply(NEG_N1), C1);


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

        // ==== SQUARE FACES ====
        Tetrad<Tuple<Apfloat>> sqr0 = new Tetrad<>(vC10, vC0, vC2, vC8);
        Triad<Apfloat> sqr0_norm = VectorMath.normalQuad(vC10, vC0, vC2, vC8, true);
        Tetrad<Tuple<Apfloat>> sqr1 = new Tetrad<>(vC5, vC2, vC1, vC11);
        Triad<Apfloat> sqr1_norm = VectorMath.normalQuad(vC5, vC2, vC1, vC11, true);
        Tetrad<Tuple<Apfloat>> sqr2 = new Tetrad<>(vC7, vC1, vC0, vC4);
        Triad<Apfloat> sqr2_norm = VectorMath.normalQuad(vC7, vC1, vC0, vC4, true);
        Tetrad<Tuple<Apfloat>> sqr3 = new Tetrad<>(vC10, vC8, vC9, vC6);
        Triad<Apfloat> sqr3_norm = VectorMath.normalQuad(vC10, vC8, vC9, vC6, true);
        Tetrad<Tuple<Apfloat>> sqr4 = new Tetrad<>(vC5, vC11, vC3, vC9);
        Triad<Apfloat> sqr4_norm = VectorMath.normalQuad(vC5, vC11, vC3, vC9, true);
        Tetrad<Tuple<Apfloat>> sqr5 = new Tetrad<>(vC7, vC4, vC6, vC3);
        Triad<Apfloat> sqr5_norm = VectorMath.normalQuad(vC7, vC4, vC6, vC3, true);
        faces_sqr = new Hexad<>(sqr0, sqr1, sqr2,sqr3,sqr4,sqr5);
        face_norms_sqr = new Hexad<>(sqr0_norm,sqr1_norm,sqr2_norm,sqr3_norm,sqr4_norm,sqr5_norm);

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC10, vC6, vC4);
        Triad<Apfloat> tri0_norm = VectorMath.normalTriple(vC10, vC6, vC4, true);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC5, vC9, vC8);
        Triad<Apfloat> tri1_norm = VectorMath.normalTriple(vC5, vC9, vC8, true);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC7, vC3, vC11);
        Triad<Apfloat> tri2_norm = VectorMath.normalTriple(vC7, vC3, vC11, true);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC10, vC4, vC0);
        Triad<Apfloat> tri3_norm = VectorMath.normalTriple(vC10, vC4, vC0, true);
        Triad<Tuple<Apfloat>> tri4 = new Triad<>(vC5, vC8, vC2);
        Triad<Apfloat> tri4_norm = VectorMath.normalTriple(vC5, vC8, vC2, true);
        Triad<Tuple<Apfloat>> tri5 = new Triad<>(vC7, vC11, vC1);
        Triad<Apfloat> tri5_norm = VectorMath.normalTriple(vC7, vC11, vC1, true);
        Triad<Tuple<Apfloat>> tri6 = new Triad<>(vC0, vC1, vC2);
        Triad<Apfloat> tri6_norm = VectorMath.normalTriple(vC0, vC1, vC2, true);
        Triad<Tuple<Apfloat>> tri7 = new Triad<>(vC3, vC6, vC9);
        Triad<Apfloat> tri7_norm = VectorMath.normalTriple(vC3, vC6, vC9, true);
        faces_tri = new Octad<>(tri0,tri1,tri2,tri3,tri4,tri5,tri6,tri7);
        face_norms_tri = new Octad<>(tri0_norm,tri1_norm,tri2_norm,tri3_norm,tri4_norm,tri5_norm,tri6_norm,tri7_norm);

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
