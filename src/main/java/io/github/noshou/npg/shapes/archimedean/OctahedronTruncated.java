package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.npg.shapes.platonic.Octahedron;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;
import static io.github.noshou.npg.nputil.VectorMath.*;
import io.github.noshou.npg.shapes.catalan.*;

/**
 * Represents a <b>Truncated Octahedron</b>.
 * <p> An Archimedean solid formed by truncating the vertices of a
 * {@link Octahedron}, it is a convex polyhedron with 14 faces: 8 regular hexagons
 * and 6 squares. It has 36 edges and 24 vertices, part of the Oh symmetry group.
 * <p> It can be circumscribed ({@link OctahedronTruncatedCanonical}) and
 * biscribed ({@link OctahedronTruncatedBiscribed}).
 * <p> It is the dual of the {@link HexahedronTetrakis}
 */
@SuppressWarnings("FieldCanBeLocal")
public abstract class OctahedronTruncated extends Shape {

    // 6 square faces
    private Hexad<Tuple<Tuple<Apfloat>>> faces_sqr;
    private Hexad<Tuple<Apfloat>> face_norms_sqr;

    // 8 Hexagonal  faces
    private Octad<Tuple<Tuple<Apfloat>>> faces_hex;
    private Octad<Tuple<Apfloat>> face_norms_hex;

    // Apfloat constants
    protected final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    protected final Apfloat N0 = new Apfloat("0", super.precision);
    protected final Apfloat N1 = new Apfloat("1", super.precision);
    protected final Apfloat N2 = new Apfloat("2", super.precision);
    protected final Apfloat SQRT2 = ApfloatMath.sqrt(N2);
    protected final Apfloat SQRT3 = ApfloatMath.sqrt(new Apfloat("3", super.precision));

    // basis vertices
    private ArrayList<Triad<Apfloat>> vertices;

    /*
     * Sets the basis vertices.
     * child classes must implement this method in their constructors.
     */
    protected void setVertices(
        Triad<Apfloat> vB0,
        Triad<Apfloat> vB1,
        Triad<Apfloat> vB2,
        Triad<Apfloat> vB3,
        Triad<Apfloat> vB4,
        Triad<Apfloat> vB5,
        Triad<Apfloat> vB6,
        Triad<Apfloat> vB7,
        Triad<Apfloat> vB8,
        Triad<Apfloat> vB9,
        Triad<Apfloat> vB10,
        Triad<Apfloat> vB11,
        Triad<Apfloat> vB12,
        Triad<Apfloat> vB13,
        Triad<Apfloat> vB14,
        Triad<Apfloat> vB15,
        Triad<Apfloat> vB16,
        Triad<Apfloat> vB17,
        Triad<Apfloat> vB18,
        Triad<Apfloat> vB19,
        Triad<Apfloat> vB20,
        Triad<Apfloat> vB21,
        Triad<Apfloat> vB22,
        Triad<Apfloat> vB23
    ) {
        ArrayList<Triad<Apfloat>> verts = new ArrayList<>();

        // ==== SCALED VERTICES ====
        verts.add(mult(normalize(vB0),  super.getRadius().toString()));
        verts.add(mult(normalize(vB1),  super.getRadius().toString()));
        verts.add(mult(normalize(vB2),  super.getRadius().toString()));
        verts.add(mult(normalize(vB3),  super.getRadius().toString()));
        verts.add(mult(normalize(vB4),  super.getRadius().toString()));
        verts.add(mult(normalize(vB5),  super.getRadius().toString()));
        verts.add(mult(normalize(vB6),  super.getRadius().toString()));
        verts.add(mult(normalize(vB7),  super.getRadius().toString()));
        verts.add(mult(normalize(vB8),  super.getRadius().toString()));
        verts.add(mult(normalize(vB9),  super.getRadius().toString()));
        verts.add(mult(normalize(vB10), super.getRadius().toString()));
        verts.add(mult(normalize(vB11), super.getRadius().toString()));
        verts.add(mult(normalize(vB12), super.getRadius().toString()));
        verts.add(mult(normalize(vB13), super.getRadius().toString()));
        verts.add(mult(normalize(vB14), super.getRadius().toString()));
        verts.add(mult(normalize(vB15), super.getRadius().toString()));
        verts.add(mult(normalize(vB16), super.getRadius().toString()));
        verts.add(mult(normalize(vB17), super.getRadius().toString()));
        verts.add(mult(normalize(vB18), super.getRadius().toString()));
        verts.add(mult(normalize(vB19), super.getRadius().toString()));
        verts.add(mult(normalize(vB20), super.getRadius().toString()));
        verts.add(mult(normalize(vB21), super.getRadius().toString()));
        verts.add(mult(normalize(vB22), super.getRadius().toString()));
        verts.add(mult(normalize(vB23), super.getRadius().toString()));
        this.vertices = verts;

        // ==== RECTANGULAR FACES ====
        Tetrad<Tuple<Apfloat>> sqr0 = new Tetrad<>(vertices.get(0), vertices.get(12), vertices.get(2), vertices.get(14));
        Triad<Apfloat> sqr0_norm = normalQuad(vertices.get(0), vertices.get(12), vertices.get(2), vertices.get(14), true);
        Tetrad<Tuple<Apfloat>> sqr1 = new Tetrad<>(vertices.get(1), vertices.get(13), vertices.get(9), vertices.get(20));
        Triad<Apfloat> sqr1_norm = normalQuad(vertices.get(1), vertices.get(13), vertices.get(9), vertices.get(20), true);
        Tetrad<Tuple<Apfloat>> sqr2 = new Tetrad<>(vertices.get(2), vertices.get(12), vertices.get(8), vertices.get(22));
        Triad<Apfloat> sqr2_norm = normalQuad(vertices.get(2), vertices.get(12), vertices.get(8), vertices.get(22), true);
        Tetrad<Tuple<Apfloat>> sqr3 = new Tetrad<>(vertices.get(3), vertices.get(15), vertices.get(11), vertices.get(23));
        Triad<Apfloat> sqr3_norm = normalQuad(vertices.get(3), vertices.get(15), vertices.get(11), vertices.get(23), true);
        Tetrad<Tuple<Apfloat>> sqr4 = new Tetrad<>(vertices.get(4), vertices.get(20), vertices.get(8), vertices.get(12));
        Triad<Apfloat> sqr4_norm = normalQuad(vertices.get(4), vertices.get(20), vertices.get(8), vertices.get(12), true);
        Tetrad<Tuple<Apfloat>> sqr5 = new Tetrad<>(vertices.get(5), vertices.get(21), vertices.get(11), vertices.get(15));
        Triad<Apfloat> sqr5_norm = normalQuad(vertices.get(5), vertices.get(21), vertices.get(11), vertices.get(15), true);
        faces_sqr = new Hexad<>(sqr0, sqr1, sqr2, sqr3, sqr4, sqr5);
        face_norms_sqr = new Hexad<>(sqr0_norm, sqr1_norm, sqr2_norm, sqr3_norm, sqr4_norm, sqr5_norm);

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vertices.get(0), vertices.get(14), vertices.get(10), vertices.get(21), vertices.get(5), vertices.get(16));
        Triad<Apfloat> hex0_norm = normalHex(vertices.get(0), vertices.get(14), vertices.get(10), vertices.get(21), vertices.get(5), vertices.get(16), true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vertices.get(1), vertices.get(15), vertices.get(3), vertices.get(13), vertices.get(9), vertices.get(22));
        Triad<Apfloat> hex1_norm = normalHex(vertices.get(1), vertices.get(15), vertices.get(3), vertices.get(13), vertices.get(9), vertices.get(22), true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vertices.get(4), vertices.get(16), vertices.get(5), vertices.get(17), vertices.get(7), vertices.get(18));
        Triad<Apfloat> hex2_norm = normalHex(vertices.get(4), vertices.get(16), vertices.get(5), vertices.get(17), vertices.get(7), vertices.get(18), true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vertices.get(6), vertices.get(19), vertices.get(7), vertices.get(18), vertices.get(2), vertices.get(12));
        Triad<Apfloat> hex3_norm = normalHex(vertices.get(6), vertices.get(19), vertices.get(7), vertices.get(18), vertices.get(2), vertices.get(12), true);
        Hexad<Tuple<Apfloat>> hex4 = new Hexad<>(vertices.get(8), vertices.get(20), vertices.get(9), vertices.get(22), vertices.get(0), vertices.get(14));
        Triad<Apfloat> hex4_norm = normalHex(vertices.get(8), vertices.get(20), vertices.get(9), vertices.get(22), vertices.get(0), vertices.get(14), true);
        Hexad<Tuple<Apfloat>> hex5 = new Hexad<>(vertices.get(10), vertices.get(23), vertices.get(11), vertices.get(21), vertices.get(1), vertices.get(13));
        Triad<Apfloat> hex5_norm = normalHex(vertices.get(10), vertices.get(23), vertices.get(11), vertices.get(21), vertices.get(1), vertices.get(13), true);
        Hexad<Tuple<Apfloat>> hex6 = new Hexad<>(vertices.get(6), vertices.get(19), vertices.get(3), vertices.get(15), vertices.get(5), vertices.get(17));
        Triad<Apfloat> hex6_norm = normalHex(vertices.get(6), vertices.get(19), vertices.get(3), vertices.get(15), vertices.get(5), vertices.get(17), true);
        Hexad<Tuple<Apfloat>> hex7 = new Hexad<>(vertices.get(2), vertices.get(18), vertices.get(7), vertices.get(23), vertices.get(11), vertices.get(21));
        Triad<Apfloat> hex7_norm = normalHex(vertices.get(2), vertices.get(18), vertices.get(7), vertices.get(23), vertices.get(11), vertices.get(21), true);
        faces_hex = new Octad<>(hex0, hex1, hex2, hex3, hex4, hex5, hex6, hex7);
        face_norms_hex = new Octad<>(hex0_norm, hex1_norm, hex2_norm, hex3_norm, hex4_norm, hex5_norm, hex6_norm, hex7_norm);

    }

    public OctahedronTruncated(
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
            // Hexagonal  faces
            Hexad<Tuple<Apfloat>> face = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_hex.fetch(i);

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
