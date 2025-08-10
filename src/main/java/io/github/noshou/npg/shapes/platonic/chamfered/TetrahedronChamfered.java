package io.github.noshou.npg.shapes.platonic.chamfered;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.npg.shapes.platonic.Tetrahedron;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b>Chamfered Tetrahedron</b>.
 * <p>The chamfered tetrahedron is a convex polyhedron constructed by beveling (chamfering) the edges
 * of a {@link Tetrahedron}. Each edge is replaced by a hexagonal face, and the original
 * vertices become smaller triangles. The result is a polyhedron with 24 edges,
 * 12 vertices, and 8 faces (4 hexagonal and 4 triangular).
 * @see <a href="https://dmccooey.com/polyhedra/ChamferedTetrahedron3.html">
 *     Chamfered Tetrahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class TetrahedronChamfered extends Shape {

    //  6 hexagonal faces
    private final Hexad<Tuple<Tuple<Apfloat>>> faces_hex;
    private final Hexad<Tuple<Apfloat>> face_norms_hex;

    // 4 triangular faces
    private final Tetrad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Tetrad<Tuple<Apfloat>> face_norms_tri;


    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat FRAC_SQRT2_over_N2 = ApfloatMath.sqrt(N2).divide(N2);
    private final Apfloat FRAC_N2minSQRT2_over_N2 = (N2.subtract(ApfloatMath.sqrt(N2))).divide(N2);

    public TetrahedronChamfered(
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
        Triad<Apfloat> vB0 = new Triad<>(
                FRAC_SQRT2_over_N2,
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2
        );
        Triad<Apfloat> vB1 = new Triad<>(
                FRAC_SQRT2_over_N2,
                FRAC_SQRT2_over_N2,
                FRAC_SQRT2_over_N2.multiply(NEG_N1)
        );
        Triad<Apfloat> vB2 = new Triad<>(
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2,
                FRAC_SQRT2_over_N2
        );
        Triad<Apfloat> vB3 = new Triad<>(
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2.multiply(NEG_N1)
        );
        Triad<Apfloat> vB4 = new Triad<>(
                FRAC_SQRT2_over_N2,
                FRAC_N2minSQRT2_over_N2,
                FRAC_SQRT2_over_N2
        );
        Triad<Apfloat> vB5 = new Triad<>(
                FRAC_SQRT2_over_N2,
                FRAC_N2minSQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2.multiply(NEG_N1)
        );
        Triad<Apfloat> vB6 = new Triad<>(
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_N2minSQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2
        );
        Triad<Apfloat> vB7 = new Triad<>(
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_N2minSQRT2_over_N2,
                FRAC_SQRT2_over_N2.multiply(NEG_N1)
        );
        Triad<Apfloat> vB8 = new Triad<>(
                FRAC_SQRT2_over_N2,
                FRAC_SQRT2_over_N2,
                FRAC_N2minSQRT2_over_N2
        );
        Triad<Apfloat> vB9 = new Triad<>(
                FRAC_SQRT2_over_N2,
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_N2minSQRT2_over_N2.multiply(NEG_N1)
        );
        Triad<Apfloat> vB10 = new Triad<>(
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_N2minSQRT2_over_N2
        );
        Triad<Apfloat> vB11 = new Triad<>(
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2,
                FRAC_N2minSQRT2_over_N2.multiply(NEG_N1)
        );
        Triad<Apfloat> vB12 = new Triad<>(
                FRAC_N2minSQRT2_over_N2,
                FRAC_SQRT2_over_N2,
                FRAC_SQRT2_over_N2
        );
        Triad<Apfloat> vB13 = new Triad<>(
                FRAC_N2minSQRT2_over_N2,
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2.multiply(NEG_N1)
        );
        Triad<Apfloat> vB14 = new Triad<>(
                FRAC_N2minSQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2
        );
        Triad<Apfloat> vB15 = new Triad<>(
                FRAC_N2minSQRT2_over_N2.multiply(NEG_N1),
                FRAC_SQRT2_over_N2,
                FRAC_SQRT2_over_N2.multiply(NEG_N1)
        );

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0 = mult(normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1 = mult(normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2 = mult(normalize(vB2), super.getRadius().toString());
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

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC0, vC9, vC5, vC1, vC8, vC4);
        Triad<Apfloat> hex0_norm = normalHex(vC0, vC9, vC5, vC1, vC8, vC4, true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC0, vC14, vC10, vC3, vC13, vC9);
        Triad<Apfloat> hex1_norm = normalHex(vC0, vC14, vC10, vC3, vC13, vC9, true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC1, vC5, vC13, vC3, vC7, vC15);
        Triad<Apfloat> hex2_norm = normalHex(vC1, vC5, vC13, vC3, vC7, vC15, true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC1, vC15, vC11, vC2, vC12, vC8);
        Triad<Apfloat> hex3_norm = normalHex(vC1, vC15, vC11, vC2, vC12, vC8, true);
        Hexad<Tuple<Apfloat>> hex4 = new Hexad<>(vC2, vC6, vC14, vC0, vC4, vC12);
        Triad<Apfloat> hex4_norm = normalHex(vC2, vC6, vC14, vC0, vC4, vC12, true);
        Hexad<Tuple<Apfloat>> hex5 = new Hexad<>(vC2, vC11, vC7, vC3, vC10, vC6);
        Triad<Apfloat> hex5_norm = normalHex(vC2, vC11, vC7, vC3, vC10, vC6, true);
        faces_hex = new Hexad<>(hex0,hex1,hex2,hex3,hex4,hex5);
        face_norms_hex = new Hexad<>(hex0_norm,hex1_norm,hex2_norm,hex3_norm,hex4_norm,hex5_norm);

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC4, vC8, vC12);
        Triad<Apfloat> tri0_norm = normalTriple(vC4,vC8,vC12,true);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC5, vC9, vC13);
        Triad<Apfloat> tri1_norm = normalTriple(vC5, vC9, vC13,true);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC6, vC10, vC14);
        Triad<Apfloat> tri2_norm = normalTriple(vC6, vC10, vC14, true);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC7, vC11, vC15);
        Triad<Apfloat> tri3_norm = normalTriple(vC7, vC11, vC15,true);
        faces_tri = new Tetrad<>(tri0,tri1,tri2,tri3);
        face_norms_tri = new Tetrad<>(tri0_norm,tri1_norm,tri2_norm,tri3_norm);
    }

    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < faces_hex.fetchSize(); i++){

            // === Check triangular faces ===
            if (i < face_norms_tri.fetchSize()) {
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


            // === Check hexagonal faces ===
            // check if each point lies within bounds of each face

            // get vertices of vertA
            Hexad<Tuple<Apfloat>> face = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x - vertA_z)
            Triad<Apfloat> m = subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> face_norm = (Triad<Apfloat>) face_norms_hex.fetch(i);
            Apfloat d = dot_prod(face_norm, m);

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
