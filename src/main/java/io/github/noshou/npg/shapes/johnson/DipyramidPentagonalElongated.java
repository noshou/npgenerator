package io.github.noshou.npg.shapes.johnson;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.*;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents an <b>Elongated Pentagonal Dipyramid</b>.
 * <p> This Johnson solid (J14) is formed by elongating a pentagonal dipyramid with a pentagonal prism inserted between
 * its two halves. It consists of 15 faces: 10 equilateral triangles and 5 squares, with 25 edges and 12 vertices.
 * The shape exhibits D5h symmetry and is convex.
 * @see <a href="https://dmccooey.com/polyhedra/ElongatedPentagonalDipyramid.html">
 *      Elongated Pentagonal Dipyramid (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class DipyramidPentagonalElongated extends Shape {

    // 5 square face
    Pentad<Tuple<Tuple<Apfloat>>> faces_sqr;
    Pentad<Tuple<Apfloat>> face_norms_sqr;

    // 10 Triangular  faces
    Decad<Tuple<Tuple<Apfloat>>> faces_tri;
    Decad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(N5);
    private final Apfloat N10 = new Apfloat("10", super.precision);
    private final Apfloat N11 = new Apfloat("11", super.precision);
    private final Apfloat N40 = new Apfloat("40", super.precision);
    private final Apfloat N185 = new Apfloat("185", super.precision);

    //  C0 = 0.309016994374947424102293417183 = (sqrt(5) - 1) / 4
    private final Apfloat C0 = (SQRT5.subtract(N1)).divide(N4);

    //  C1 = sqrt(2 * (5 - sqrt(5))) / 4
    private final Apfloat C1 = ApfloatMath.sqrt(N2.multiply(N5.subtract(SQRT5))).divide(N4);

    //  C2 = (1 + sqrt(5)) / 4
    private final Apfloat C2 = (N1.add(SQRT5)).divide(N4);

    //  C3 = sqrt(2 * (5 + sqrt(5))) / 4
    private final Apfloat C3 = ApfloatMath.sqrt(N2.multiply(N5.add(SQRT5))).divide(N4);

    //  C4 = sqrt(10 * (185 + 11 * sqrt(5))) / 40
    private final Apfloat C4 = ApfloatMath.sqrt(N10.multiply(N185.add(N11.multiply(SQRT5)))).divide(N40);


    public DipyramidPentagonalElongated(
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
        Triad<Apfloat> vB0  = new Triad<>(C3, C0, C1);
        Triad<Apfloat> vB1  = new Triad<>(C3, C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB2  = new Triad<>(C3.multiply(NEG_N1), C0, C1);
        Triad<Apfloat> vB3  = new Triad<>(C3.multiply(NEG_N1), C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB4  = new Triad<>(C1, C2.multiply(NEG_N1), C1);
        Triad<Apfloat> vB5  = new Triad<>(C1, C2.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB6  = new Triad<>(C1.multiply(NEG_N1), C2.multiply(NEG_N1), C1);
        Triad<Apfloat> vB7  = new Triad<>(C1.multiply(NEG_N1), C2.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB8  = new Triad<>(N0, N1, C1);
        Triad<Apfloat> vB9  = new Triad<>(N0, N1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(N0, N0, C4);
        Triad<Apfloat> vB11 = new Triad<>(N0, N0, C4.multiply(NEG_N1));

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = mult(normalize(vB0),  super.getRadius().toString());
        Triad<Apfloat> vC1  = mult(normalize(vB1),  super.getRadius().toString());
        Triad<Apfloat> vC2  = mult(normalize(vB2),  super.getRadius().toString());
        Triad<Apfloat> vC3  = mult(normalize(vB3),  super.getRadius().toString());
        Triad<Apfloat> vC4  = mult(normalize(vB4),  super.getRadius().toString());
        Triad<Apfloat> vC5  = mult(normalize(vB5),  super.getRadius().toString());
        Triad<Apfloat> vC6  = mult(normalize(vB6),  super.getRadius().toString());
        Triad<Apfloat> vC7  = mult(normalize(vB7),  super.getRadius().toString());
        Triad<Apfloat> vC8  = mult(normalize(vB8),  super.getRadius().toString());
        Triad<Apfloat> vC9  = mult(normalize(vB9),  super.getRadius().toString());
        Triad<Apfloat> vC10 = mult(normalize(vB10), super.getRadius().toString());
        Triad<Apfloat> vC11 = mult(normalize(vB11), super.getRadius().toString());

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0  = new Triad<>(vC10, vC0, vC8);
        Triad<Apfloat> tri0_norm  = normalTriple(vC10, vC0, vC8, true);
        Triad<Tuple<Apfloat>> tri1  = new Triad<>(vC10, vC8, vC2);
        Triad<Apfloat> tri1_norm  = normalTriple(vC10, vC8, vC2, true);
        Triad<Tuple<Apfloat>> tri2  = new Triad<>(vC10, vC2, vC6);
        Triad<Apfloat> tri2_norm  = normalTriple(vC10, vC2, vC6, true);
        Triad<Tuple<Apfloat>> tri3  = new Triad<>(vC10, vC6, vC4);
        Triad<Apfloat> tri3_norm  = normalTriple(vC10, vC6, vC4, true);
        Triad<Tuple<Apfloat>> tri4  = new Triad<>(vC10, vC4, vC0);
        Triad<Apfloat> tri4_norm  = normalTriple(vC10, vC4, vC0, true);
        Triad<Tuple<Apfloat>> tri5  = new Triad<>(vC11, vC1, vC5);
        Triad<Apfloat> tri5_norm  = normalTriple(vC11, vC1, vC5, true);
        Triad<Tuple<Apfloat>> tri6  = new Triad<>(vC11, vC5, vC7);
        Triad<Apfloat> tri6_norm  = normalTriple(vC11, vC5, vC7, true);
        Triad<Tuple<Apfloat>> tri7  = new Triad<>(vC11, vC7, vC3);
        Triad<Apfloat> tri7_norm  = normalTriple(vC11, vC7, vC3, true);
        Triad<Tuple<Apfloat>> tri8  = new Triad<>(vC11, vC3, vC9);
        Triad<Apfloat> tri8_norm  = normalTriple(vC11, vC3, vC9, true);
        Triad<Tuple<Apfloat>> tri9  = new Triad<>(vC11, vC9, vC1);
        Triad<Apfloat> tri9_norm  = normalTriple(vC11, vC9, vC1, true);
        faces_tri = new Decad<>(
                tri0, tri1, tri2, tri3, tri4,
                tri5, tri6, tri7, tri8, tri9
        );
        face_norms_tri = new Decad<>(
                tri0_norm, tri1_norm, tri2_norm, tri3_norm, tri4_norm,
                tri5_norm, tri6_norm, tri7_norm, tri8_norm, tri9_norm

        );

        // ==== SQUARE FACES ====
        Tetrad<Tuple<Apfloat>> sqr0 = new Tetrad<>(vC0, vC1, vC9, vC8);
        Triad<Apfloat> sqr0_norm = normalQuad(vC0, vC1, vC9, vC8, true);
        Tetrad<Tuple<Apfloat>> sqr1 = new Tetrad<>(vC8, vC9, vC3, vC2);
        Triad<Apfloat> sqr1_norm = normalQuad(vC8, vC9, vC3, vC2, true);
        Tetrad<Tuple<Apfloat>> sqr2 = new Tetrad<>(vC2, vC3, vC7, vC6);
        Triad<Apfloat> sqr2_norm = normalQuad(vC2, vC3, vC7, vC6, true);
        Tetrad<Tuple<Apfloat>> sqr3 = new Tetrad<>(vC6, vC7, vC5, vC4);
        Triad<Apfloat> sqr3_norm = normalQuad(vC6, vC7, vC5, vC4, true);
        Tetrad<Tuple<Apfloat>> sqr4 = new Tetrad<>(vC4, vC5, vC1, vC0);
        Triad<Apfloat> sqr4_norm = normalQuad(vC4, vC5, vC1, vC0, true);
        faces_sqr = new Pentad<>(sqr0, sqr1, sqr2, sqr3, sqr4);
        face_norms_sqr = new Pentad<>(sqr0_norm, sqr1_norm,sqr2_norm,sqr3_norm,sqr4_norm);
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 10; i++) {

            // Square faces
            if (i < 5){
                Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_sqr.fetch(i);

                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Apfloat d = dot_prod(norm, m);

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
            Triad<Apfloat> m = subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}

