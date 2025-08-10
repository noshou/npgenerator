package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.npg.shapes.platonic.Tetrahedron;
import io.github.noshou.npg.shapes.archimedean.TetrahedronTruncated;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b>Triakis Tetrahedron</b>.
 * <p> This polyhedron is the Catalan solid dual of the {@link TetrahedronTruncated}.
 * It is formed by attaching a {@link Tetrahedron} to each face
 * of a regular tetrahedron, producing a shape with 12 isosceles triangular faces,
 * 8 vertices, and 18 edges.
 * @see <a href="https://dmccooey.com/polyhedra/TriakisTetrahedron.html">
 *      Triakis Tetrahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class TetrahedronTriakis extends Shape {

    // 12 triangular faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Dodecad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat Constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N9 = new Apfloat("9", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N20 = new Apfloat("20", super.precision);
    private final Apfloat SQRT2 = ApfloatMath.sqrt(new Apfloat("2", super.precision));

    // C0 = 9 * sqrt(2) / 20
    private final Apfloat C0 = N9.multiply(SQRT2).divide(N20);

    // C1 3 * sqrt(2) / 4
    private final Apfloat C1 = N3.multiply(SQRT2).divide(N4);

    public TetrahedronTriakis (
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

        // ==== BASIS VERTICES ====
        Triad<Apfloat> vB0 = new Triad<>(C1, C1, C1);
        Triad<Apfloat> vB1 = new Triad<>( C1,  NEG_N1.multiply(C1),  NEG_N1.multiply(C1));
        Triad<Apfloat> vB2 = new Triad<>(NEG_N1.multiply(C1),  NEG_N1.multiply(C1),   C1);
        Triad<Apfloat> vB3 = new Triad<>(NEG_N1.multiply(C1),   C1,  NEG_N1.multiply(C1));
        Triad<Apfloat> vB4 = new Triad<>(   C0,  NEG_N1.multiply(C0),    C0);
        Triad<Apfloat> vB5 = new Triad<>(   C0,    C0,  NEG_N1.multiply(C0));
        Triad<Apfloat> vB6 = new Triad<>(NEG_N1.multiply(C0),    C0,    C0);
        Triad<Apfloat> vB7 = new Triad<>(NEG_N1.multiply(C0),  NEG_N1.multiply(C0),  NEG_N1.multiply(C0));

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = mult(normalize(vB0),  super.getRadius().toString());
        Triad<Apfloat> vC1  = mult(normalize(vB1),  super.getRadius().toString());
        Triad<Apfloat> vC2  = mult(normalize(vB2),  super.getRadius().toString());
        Triad<Apfloat> vC3  = mult(normalize(vB3),  super.getRadius().toString());
        Triad<Apfloat> vC4  = mult(normalize(vB4),  super.getRadius().toString());
        Triad<Apfloat> vC5  = mult(normalize(vB5),  super.getRadius().toString());
        Triad<Apfloat> vC6  = mult(normalize(vB6),  super.getRadius().toString());
        Triad<Apfloat> vC7  = mult(normalize(vB7),  super.getRadius().toString());

        // ==== TRIANGULAR VERTICES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC4, vC0, vC2);
        Triad<Apfloat> tri0_norm = normalTriple(vC4, vC0, vC2, true);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC4, vC2, vC1);
        Triad<Apfloat> tri1_norm = normalTriple(vC4, vC2, vC1, true);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC4, vC1, vC0);
        Triad<Apfloat> tri2_norm = normalTriple(vC4, vC1, vC0, true);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC5, vC0, vC1);
        Triad<Apfloat> tri3_norm = normalTriple(vC5, vC0, vC1, true);
        Triad<Tuple<Apfloat>> tri4 = new Triad<>(vC5, vC1, vC3);
        Triad<Apfloat> tri4_norm = normalTriple(vC5, vC1, vC3, true);
        Triad<Tuple<Apfloat>> tri5 = new Triad<>(vC5, vC3, vC0);
        Triad<Apfloat> tri5_norm = normalTriple(vC5, vC3, vC0, true);
        Triad<Tuple<Apfloat>> tri6 = new Triad<>(vC6, vC0, vC3);
        Triad<Apfloat> tri6_norm = normalTriple(vC6, vC0, vC3, true);
        Triad<Tuple<Apfloat>> tri7 = new Triad<>(vC6, vC3, vC2);
        Triad<Apfloat> tri7_norm = normalTriple(vC6, vC3, vC2, true);
        Triad<Tuple<Apfloat>> tri8 = new Triad<>(vC6, vC2, vC0);
        Triad<Apfloat> tri8_norm = normalTriple(vC6, vC2, vC0, true);
        Triad<Tuple<Apfloat>> tri9 = new Triad<>(vC7, vC1, vC2);
        Triad<Apfloat> tri9_norm = normalTriple(vC7, vC1, vC2, true);
        Triad<Tuple<Apfloat>> tri10 = new Triad<>(vC7, vC2, vC3);
        Triad<Apfloat> tri10_norm = normalTriple(vC7, vC2, vC3, true);
        Triad<Tuple<Apfloat>> tri11 = new Triad<>(vC7, vC3, vC1);
        Triad<Apfloat> tri11_norm = normalTriple(vC7, vC3, vC1, true);
        faces_tri = new Dodecad<>(
                tri0, tri1, tri2, tri3,
                tri4, tri5, tri6, tri7,
                tri8, tri9, tri10, tri11
        );
        face_norms_tri = new Dodecad<>(
                tri0_norm, tri1_norm, tri2_norm, tri3_norm,
                tri4_norm, tri5_norm, tri6_norm, tri7_norm,
                tri8_norm, tri9_norm, tri10_norm, tri11_norm
        );
    }
    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 12; i++) {
            // triangular faces  faces
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
