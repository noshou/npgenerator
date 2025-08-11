package io.github.noshou.npg.shapes.johnson;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b>Bilunabirotunda</b>
 * <p> A Johnson solid characterized by 14 vertices, 16 faces, and 30 edges.
 * <p> The Bilunabirotunda is composed of a combination of polygons, including pentagons and
 * triangles, arranged in a specific manner that forms this unique convex polyhedron.
 * It belongs to the family of Johnson solids, which are strictly convex polyhedra
 * with regular faces but are not uniform like Platonic or Archimedean solids.
 * @see <a href="https://dmccooey.com/polyhedra/Bilunabirotunda.html">
 *      Bilunabirotunda (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class Bilunabirotunda extends Shape {

    // 4 pentagonal faces
    Tetrad<Tuple<Tuple<Apfloat>>> faces_pnt;
    Tetrad<Tuple<Apfloat>> face_norms_pnt;

    // 2 square faces
    Dyad<Tuple<Tuple<Apfloat>>> faces_sqr;
    Dyad<Tuple<Apfloat>> face_norms_sqr;

    // 8 triangular faces
    Octad<Tuple<Tuple<Apfloat>>> faces_tri;
    Octad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1    = new Apfloat("-1", super.precision);
    private final Apfloat N0   = new Apfloat("0", super.precision);
    private final Apfloat N1    = new Apfloat("1", super.precision);
    private final Apfloat N3  = new Apfloat("3", super.precision);
    private final Apfloat N4   = new Apfloat("4", super.precision);
    private final Apfloat N5   = new Apfloat("5", super.precision);
    private final Apfloat FRAC_1_over_2   = new Apfloat("0.5", super.precision);

    // golden ratio based constants
    private final Apfloat C0 = ApfloatMath.sqrt(N5).add(N1).divide(N4);                      // (1 + sqrt(5)) / 4
    private final Apfloat C1 = ApfloatMath.sqrt(N5).add(N3).divide(N4);                    // (3 + sqrt(5)) / 4

    public Bilunabirotunda(
            @NotNull String radius,
            @NotNull String radius_type,
            @NotNull LatticeType lattice_type,
            int precision,
            Polyad<Atom> basis,
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
        Triad<Apfloat> vB0  = new Triad<>( FRAC_1_over_2, C0, FRAC_1_over_2);
        Triad<Apfloat> vB1  = new Triad<>( FRAC_1_over_2, C0, FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB2  = new Triad<>( FRAC_1_over_2, C0.multiply(NEG_N1), FRAC_1_over_2);
        Triad<Apfloat> vB3  = new Triad<>( FRAC_1_over_2, C0.multiply(NEG_N1), FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB4  = new Triad<>( FRAC_1_over_2.multiply(NEG_N1),C0, FRAC_1_over_2);
        Triad<Apfloat> vB5  = new Triad<>( FRAC_1_over_2.multiply(NEG_N1), C0,  FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB6  = new Triad<>( FRAC_1_over_2.multiply(NEG_N1), C0.multiply(NEG_N1), FRAC_1_over_2);
        Triad<Apfloat> vB7  = new Triad<>( FRAC_1_over_2.multiply(NEG_N1), C0.multiply(NEG_N1), FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB8  = new Triad<>( C1, FRAC_1_over_2, N0);
        Triad<Apfloat> vB9  = new Triad<>( C1, FRAC_1_over_2.multiply(NEG_N1), N0);
        Triad<Apfloat> vB10 = new Triad<>( C1.multiply(NEG_N1), FRAC_1_over_2,N0);
        Triad<Apfloat> vB11 = new Triad<>( C1.multiply(NEG_N1), FRAC_1_over_2.multiply(NEG_N1), N0);
        Triad<Apfloat> vB12 = new Triad<>( N0, N0, C0);
        Triad<Apfloat> vB13 = new Triad<>( N0, N0, C0.multiply(NEG_N1));

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
        Triad<Apfloat> vC12 = mult(normalize(vB12), super.getRadius().toString());
        Triad<Apfloat> vC13 = mult(normalize(vB13), super.getRadius().toString());

        // ==== PENTAGONAL FACES ====
        Pentad<Tuple<Apfloat>> pnt0 = new Pentad<>(vC12, vC2,  vC9,  vC8,  vC0);
        Triad<Apfloat> pnt0_norm = normalPent(vC12, vC2, vC9, vC8, vC0, true);
        Pentad<Tuple<Apfloat>> pnt1 = new Pentad<>(vC12, vC4,  vC10, vC11, vC6);
        Triad<Apfloat> pnt1_norm = normalPent(vC12, vC4,  vC10, vC11, vC6, true);
        Pentad<Tuple<Apfloat>> pnt2 = new Pentad<>(vC13, vC1,  vC8,  vC9,  vC3);
        Triad<Apfloat> pnt2_norm = normalPent(vC13, vC1,  vC8,  vC9,  vC3, true);
        Pentad<Tuple<Apfloat>> pnt3 = new Pentad<>(vC13, vC7,  vC11, vC10, vC5);
        Triad<Apfloat> pnt3_norm = normalPent(vC13, vC7,  vC11, vC10, vC5, true);
        faces_pnt = new Tetrad<>(pnt0, pnt1, pnt2, pnt3);
        face_norms_pnt = new Tetrad<>(pnt0_norm, pnt1_norm, pnt2_norm, pnt3_norm);

        // ==== QUADRILATERAL FACES ====
        Tetrad<Tuple<Apfloat>> sqr0 = new Tetrad<>(vC0, vC1, vC5, vC4);
        Triad<Apfloat> sqr0_norm = normalQuad(vC0, vC1, vC5, vC4, true);
        Tetrad<Tuple<Apfloat>> sqr1 = new Tetrad<>(vC3, vC2, vC6, vC7);
        Triad<Apfloat> sqr1_norm = normalQuad(vC3, vC2, vC6, vC7, true);
        faces_sqr = new Dyad<>(sqr0, sqr1);
        face_norms_sqr = new Dyad<>(sqr0_norm, sqr1_norm);

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC8,  vC1,  vC0);
        Triad<Apfloat> tri0_norm  = normalTriple(vC8, vC1, vC0, true);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC9,  vC2,  vC3);
        Triad<Apfloat> tri1_norm  = normalTriple(vC9,  vC2,  vC3, true);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC10, vC4,  vC5);
        Triad<Apfloat> tri2_norm  = normalTriple(vC10, vC4,  vC5, true);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC11, vC7,  vC6);
        Triad<Apfloat> tri3_norm  = normalTriple(vC11, vC7,  vC6, true);
        Triad<Tuple<Apfloat>> tri4 = new Triad<>(vC12, vC0,  vC4);
        Triad<Apfloat> tri4_norm  = normalTriple(vC12, vC0,  vC4, true);
        Triad<Tuple<Apfloat>> tri5 = new Triad<>(vC12, vC6,  vC2);
        Triad<Apfloat> tri5_norm  = normalTriple(vC12, vC6,  vC2, true);
        Triad<Tuple<Apfloat>> tri6 = new Triad<>(vC13, vC3,  vC7);
        Triad<Apfloat> tri6_norm  = normalTriple(vC13, vC3,  vC7, true);
        Triad<Tuple<Apfloat>> tri7 = new Triad<>(vC13, vC5,  vC1);
        Triad<Apfloat> tri7_norm  = normalTriple(vC13, vC5,  vC1, true);

        faces_tri = new Octad<>(tri0, tri1, tri2, tri3, tri4, tri5, tri6, tri7);
        face_norms_tri = new Octad<>(tri0_norm, tri1_norm, tri2_norm, tri3_norm, tri4_norm, tri5_norm, tri6_norm, tri7_norm);
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 8; i++) {

            // square faces
            if (i == 0 || i == 1) {
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

            // pentagonal faces
            if (i < 3) {
                Pentad<Tuple<Apfloat>> face = (Pentad<Tuple<Apfloat>>) faces_pnt.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_pnt.fetch(i);

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

            // triangular faces
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
