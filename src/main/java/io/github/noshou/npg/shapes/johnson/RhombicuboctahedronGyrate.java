package io.github.noshou.npg.shapes.johnson;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.npg.shapes.archimedean.*;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * <p>Represents a <b>Gyrate Rhombicuboctahedron</b> (aka Elongated Square Gyrobicupola).</p>
 * <p>This solid is formed by rotating one of the square copula of the
 * {@link Rhombicuboctahedron} by 45°.</p>
 * <p>It has the following properties:</p>
 * <ul>
 *   <li>Faces: 18 squares and 8 triangles</li>
 *   <li>Vertices: 24</li>
 *   <li>Edges: 48</li>
 * </ul>
 *  @see <a href="https://dmccooey.com/polyhedra/ElongatedSquareGyrobicupola.html">
 *      Gyrate Rhombicuboctahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class RhombicuboctahedronGyrate extends Shape {

    // 18 square faces
    private final Octakaidecad<Tuple<Tuple<Apfloat>>> faces_sqr;
    private final Octakaidecad<Tuple<Apfloat>> face_norms_sqr;

    // 8 triangular faces
    private final Octad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Octad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat SQRT2 = ApfloatMath.sqrt(N2);
    private final Apfloat HALF = N1.divide(N2);
    private final Apfloat NEG_HALF = HALF.multiply(NEG_N1);

    // C0 = sqrt(2) / 2
    private final Apfloat C0 = SQRT2.divide(N2);
    private final Apfloat NEG_C0 = C0.multiply(NEG_N1);

    // C1 = (1 + sqrt(2)) / 2
    private final Apfloat C1 = (N1.add(SQRT2)).divide(N2);
    private final Apfloat NEG_C1 = C1.multiply(NEG_N1);

    public RhombicuboctahedronGyrate(
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
        Triad<Apfloat> vB0 = new Triad<>(C0, N0, C1);
        Triad<Apfloat> vB1 = new Triad<>(NEG_C0, N0, C1);
        Triad<Apfloat> vB2 = new Triad<>(N0, C0, C1);
        Triad<Apfloat> vB3 = new Triad<>(N0, NEG_C0, C1);
        Triad<Apfloat> vB4 = new Triad<>(HALF, HALF, NEG_C1);
        Triad<Apfloat> vB5 = new Triad<>(HALF, NEG_HALF, NEG_C1);
        Triad<Apfloat> vB6 = new Triad<>(NEG_HALF, HALF, NEG_C1);
        Triad<Apfloat> vB7 = new Triad<>(NEG_HALF, NEG_HALF, NEG_C1);
        Triad<Apfloat> vB8 = new Triad<>(C1, HALF, HALF);
        Triad<Apfloat> vB9 = new Triad<>(C1, HALF, NEG_HALF);
        Triad<Apfloat> vB10 = new Triad<>(C1, NEG_HALF, HALF);
        Triad<Apfloat> vB11 = new Triad<>(C1, NEG_HALF, NEG_HALF);
        Triad<Apfloat> vB12 = new Triad<>(NEG_C1, HALF, HALF);
        Triad<Apfloat> vB13 = new Triad<>(NEG_C1, HALF, NEG_HALF);
        Triad<Apfloat> vB14 = new Triad<>(NEG_C1, NEG_HALF, HALF);
        Triad<Apfloat> vB15 = new Triad<>(NEG_C1, NEG_HALF, NEG_HALF);
        Triad<Apfloat> vB16 = new Triad<>(HALF, C1, HALF);
        Triad<Apfloat> vB17 = new Triad<>(HALF, C1, NEG_HALF);
        Triad<Apfloat> vB18 = new Triad<>(HALF, NEG_C1, HALF);
        Triad<Apfloat> vB19 = new Triad<>(HALF, NEG_C1, NEG_HALF);
        Triad<Apfloat> vB20 = new Triad<>(NEG_HALF, C1, HALF);
        Triad<Apfloat> vB21 = new Triad<>(NEG_HALF, C1, NEG_HALF);
        Triad<Apfloat> vB22 = new Triad<>(NEG_HALF, NEG_C1, HALF);
        Triad<Apfloat> vB23 = new Triad<>(NEG_HALF, NEG_C1, NEG_HALF);

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = mult(normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1  = mult(normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2  = mult(normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3  = mult(normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4  = mult(normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5  = mult(normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6  = mult(normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7  = mult(normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8  = mult(normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9  = mult(normalize(vB9), super.getRadius().toString());
        Triad<Apfloat> vC10 = mult(normalize(vB10), super.getRadius().toString());
        Triad<Apfloat> vC11 = mult(normalize(vB11), super.getRadius().toString());
        Triad<Apfloat> vC12 = mult(normalize(vB12), super.getRadius().toString());
        Triad<Apfloat> vC13 = mult(normalize(vB13), super.getRadius().toString());
        Triad<Apfloat> vC14 = mult(normalize(vB14), super.getRadius().toString());
        Triad<Apfloat> vC15 = mult(normalize(vB15), super.getRadius().toString());
        Triad<Apfloat> vC16 = mult(normalize(vB16), super.getRadius().toString());
        Triad<Apfloat> vC17 = mult(normalize(vB17), super.getRadius().toString());
        Triad<Apfloat> vC18 = mult(normalize(vB18), super.getRadius().toString());
        Triad<Apfloat> vC19 = mult(normalize(vB19), super.getRadius().toString());
        Triad<Apfloat> vC20 = mult(normalize(vB20), super.getRadius().toString());
        Triad<Apfloat> vC21 = mult(normalize(vB21), super.getRadius().toString());
        Triad<Apfloat> vC22 = mult(normalize(vB22), super.getRadius().toString());
        Triad<Apfloat> vC23 = mult(normalize(vB23), super.getRadius().toString());

        // ==== QUADRILATERAL FACES ====
        Tetrad<Tuple<Apfloat>> sqr0  = new Tetrad<>(vC0, vC3, vC18, vC10);
        Triad<Apfloat> sqr0_norm    = normalQuad(vC0, vC3, vC18, vC10, true);

        Tetrad<Tuple<Apfloat>> sqr1  = new Tetrad<>(vC1, vC2, vC20, vC12);
        Triad<Apfloat> sqr1_norm    = normalQuad(vC1, vC2, vC20, vC12, true);

        Tetrad<Tuple<Apfloat>> sqr2  = new Tetrad<>(vC2, vC0, vC8, vC16);
        Triad<Apfloat> sqr2_norm    = normalQuad(vC2, vC0, vC8, vC16, true);

        Tetrad<Tuple<Apfloat>> sqr3  = new Tetrad<>(vC3, vC1, vC14, vC22);
        Triad<Apfloat> sqr3_norm    = normalQuad(vC3, vC1, vC14, vC22, true);

        Tetrad<Tuple<Apfloat>> sqr4  = new Tetrad<>(vC4, vC6, vC21, vC17);
        Triad<Apfloat> sqr4_norm    = normalQuad(vC4, vC6, vC21, vC17, true);

        Tetrad<Tuple<Apfloat>> sqr5  = new Tetrad<>(vC5, vC4, vC9, vC11);
        Triad<Apfloat> sqr5_norm    = normalQuad(vC5, vC4, vC9, vC11, true);

        Tetrad<Tuple<Apfloat>> sqr6  = new Tetrad<>(vC6, vC7, vC15, vC13);
        Triad<Apfloat> sqr6_norm    = normalQuad(vC6, vC7, vC15, vC13, true);

        Tetrad<Tuple<Apfloat>> sqr7  = new Tetrad<>(vC7, vC5, vC19, vC23);
        Triad<Apfloat> sqr7_norm    = normalQuad(vC7, vC5, vC19, vC23, true);

        Tetrad<Tuple<Apfloat>> sqr8  = new Tetrad<>(vC16, vC17, vC21, vC20);
        Triad<Apfloat> sqr8_norm    = normalQuad(vC16, vC17, vC21, vC20, true);

        Tetrad<Tuple<Apfloat>> sqr9  = new Tetrad<>(vC17, vC16, vC8, vC9);
        Triad<Apfloat> sqr9_norm    = normalQuad(vC17, vC16, vC8, vC9, true);

        Tetrad<Tuple<Apfloat>> sqr10 = new Tetrad<>(vC18, vC19, vC11, vC10);
        Triad<Apfloat> sqr10_norm   = normalQuad(vC18, vC19, vC11, vC10, true);

        Tetrad<Tuple<Apfloat>> sqr11 = new Tetrad<>(vC19, vC18, vC22, vC23);
        Triad<Apfloat> sqr11_norm   = normalQuad(vC19, vC18, vC22, vC23, true);

        Tetrad<Tuple<Apfloat>> sqr12 = new Tetrad<>(vC8, vC10, vC11, vC9);
        Triad<Apfloat> sqr12_norm   = normalQuad(vC8, vC10, vC11, vC9, true);

        Tetrad<Tuple<Apfloat>> sqr13 = new Tetrad<>(vC15, vC14, vC12, vC13);
        Triad<Apfloat> sqr13_norm   = normalQuad(vC15, vC14, vC12, vC13, true);

        Tetrad<Tuple<Apfloat>> sqr14 = new Tetrad<>(vC20, vC21, vC13, vC12);
        Triad<Apfloat> sqr14_norm   = normalQuad(vC20, vC21, vC13, vC12, true);

        Tetrad<Tuple<Apfloat>> sqr15 = new Tetrad<>(vC23, vC22, vC14, vC15);
        Triad<Apfloat> sqr15_norm   = normalQuad(vC23, vC22, vC14, vC15, true);

        Tetrad<Tuple<Apfloat>> sqr16 = new Tetrad<>(vC0, vC2, vC1, vC3);
        Triad<Apfloat> sqr16_norm   = normalQuad(vC0, vC2, vC1, vC3, true);

        Tetrad<Tuple<Apfloat>> sqr17 = new Tetrad<>(vC7, vC6, vC4, vC5);
        Triad<Apfloat> sqr17_norm   = normalQuad(vC7, vC6, vC4, vC5, true);

        faces_sqr = new Octakaidecad<>(
                sqr0, sqr1, sqr2, sqr3, sqr4, sqr5,
                sqr6, sqr7, sqr8, sqr9, sqr10, sqr11,
                sqr12, sqr13, sqr14, sqr15, sqr16, sqr17
        );
        face_norms_sqr = new Octakaidecad<>(
                sqr0_norm, sqr1_norm, sqr2_norm, sqr3_norm, sqr4_norm, sqr5_norm,
                sqr6_norm, sqr7_norm, sqr8_norm, sqr9_norm, sqr10_norm, sqr11_norm,
                sqr12_norm, sqr13_norm, sqr14_norm, sqr15_norm, sqr16_norm, sqr17_norm
        );

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC0, vC8, vC16);
        Triad<Apfloat> tri0_norm = normalTriple(vC0, vC8, vC16, true);

        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC1, vC12, vC14);
        Triad<Apfloat> tri1_norm = normalTriple(vC1, vC12, vC14, true);

        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC2, vC16, vC20);
        Triad<Apfloat> tri2_norm = normalTriple(vC2, vC16, vC20, true);

        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC3, vC22, vC18);
        Triad<Apfloat> tri3_norm = normalTriple(vC3, vC22, vC18, true);

        Triad<Tuple<Apfloat>> tri4 = new Triad<>(vC4, vC17, vC9);
        Triad<Apfloat> tri4_norm = normalTriple(vC4, vC17, vC9, true);

        Triad<Tuple<Apfloat>> tri5 = new Triad<>(vC5, vC11, vC19);
        Triad<Apfloat> tri5_norm = normalTriple(vC5, vC11, vC19, true);

        Triad<Tuple<Apfloat>> tri6 = new Triad<>(vC6, vC13, vC21);
        Triad<Apfloat> tri6_norm = normalTriple(vC6, vC13, vC21, true);

        Triad<Tuple<Apfloat>> tri7 = new Triad<>(vC7, vC23, vC15);
        Triad<Apfloat> tri7_norm = normalTriple(vC7, vC23, vC15, true);

        faces_tri = new Octad<>(
                tri0, tri1, tri2, tri3,
                tri4, tri5, tri6, tri7
        );
        face_norms_tri = new Octad<>(
                tri0_norm, tri1_norm, tri2_norm, tri3_norm,
                tri4_norm, tri5_norm, tri6_norm, tri7_norm
        );
    }

    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_sqr.fetchSize(); i++){

            // === Check triangular faces ===
            if (i < 8) {
                Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.fetch(i);
                Apfloat d = dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check square faces ===
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
            Triad<Apfloat> m = subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_sqr.fetch(i);
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
