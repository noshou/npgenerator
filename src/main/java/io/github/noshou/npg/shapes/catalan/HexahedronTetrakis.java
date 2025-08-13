package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.npg.shapes.platonic.Cube;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;
import static io.github.noshou.npg.nputil.VectorMath.*;
import io.github.noshou.npg.shapes.archimedean.*;

/**
 * Represents a <b>Tetrakis Hexahedron</b>.
 * <p> A Catalan solid formed by augmenting each face of a
 * {@link Cube} (hexahedron) with a square pyramid, it is a convex polyhedron
 * with 24 isosceles triangular faces, 36 edges, and 14 vertices, belonging
 * to the Oh symmetry group.
 * <p> It can be circumscribed ({@link HexahedronTetrakisCanonical}) and
 * biscribed ({@link HexahedronTetrakisBiscribed}).
 * <p> It is the dual of the {@link OctahedronTruncatedBiscribed}.
 */
public abstract class HexahedronTetrakis extends Shape {

    // 24 triangular faces
    private final ArrayList<Triad<Tuple<Apfloat>>> faces_tri;
    private final ArrayList<Triad<Apfloat>> face_norms_tri;

    // basis vertices
    private ArrayList<Triad<Apfloat>> vertices;

    // Apfloat constants
    protected final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    protected final Apfloat N0 = new Apfloat("0", super.precision);
    protected final Apfloat N1 = new Apfloat("1", super.precision);
    protected final Apfloat N2 = new Apfloat("2", super.precision);
    protected final Apfloat N3 = new Apfloat("3", super.precision);
    protected final Apfloat N4 = new Apfloat("4", super.precision);
    protected final Apfloat N8 = new Apfloat("8", super.precision);
    protected final Apfloat N9 = new Apfloat("9", super.precision);
    protected final Apfloat SQRT2 = ApfloatMath.sqrt(N2);
    protected final Apfloat SQRT3 = ApfloatMath.sqrt(N3);


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
            Triad<Apfloat> vB13
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
        this.vertices = verts;
    }

    public HexahedronTetrakis (
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

        // ==== TRIANGULAR FACES ====

        faces_tri = new ArrayList<>();
        face_norms_tri = new ArrayList<>();

        faces_tri.add(new Triad<>(vertices.get(0), vertices.get(10), vertices.get(12)));
        face_norms_tri.add(normalTriple(vertices.get(0), vertices.get(10), vertices.get(12), true));

        faces_tri.add(new Triad<>(vertices.get(0), vertices.get(12), vertices.get(8)));
        face_norms_tri.add(normalTriple(vertices.get(0), vertices.get(12), vertices.get(8), true));

        faces_tri.add(new Triad<>(vertices.get(0), vertices.get(8), vertices.get(6)));
        face_norms_tri.add(normalTriple(vertices.get(0), vertices.get(8), vertices.get(6), true));

        faces_tri.add(new Triad<>(vertices.get(1), vertices.get(7), vertices.get(9)));
        face_norms_tri.add(normalTriple(vertices.get(1), vertices.get(7), vertices.get(9), true));

        faces_tri.add(new Triad<>(vertices.get(1), vertices.get(9), vertices.get(13)));
        face_norms_tri.add(normalTriple(vertices.get(1), vertices.get(9), vertices.get(13), true));

        faces_tri.add(new Triad<>(vertices.get(1), vertices.get(13), vertices.get(11)));
        face_norms_tri.add(normalTriple(vertices.get(1), vertices.get(13), vertices.get(11), true));

        faces_tri.add(new Triad<>(vertices.get(1), vertices.get(11), vertices.get(7)));
        face_norms_tri.add(normalTriple(vertices.get(1), vertices.get(11), vertices.get(7), true));

        faces_tri.add(new Triad<>(vertices.get(2), vertices.get(6), vertices.get(8)));
        face_norms_tri.add(normalTriple(vertices.get(2), vertices.get(6), vertices.get(8), true));

        faces_tri.add(new Triad<>(vertices.get(2), vertices.get(8), vertices.get(9)));
        face_norms_tri.add(normalTriple(vertices.get(2), vertices.get(8), vertices.get(9), true));

        faces_tri.add(new Triad<>(vertices.get(2), vertices.get(9), vertices.get(7)));
        face_norms_tri.add(normalTriple(vertices.get(2), vertices.get(9), vertices.get(7), true));

        faces_tri.add(new Triad<>(vertices.get(2), vertices.get(7), vertices.get(6)));
        face_norms_tri.add(normalTriple(vertices.get(2), vertices.get(7), vertices.get(6), true));

        faces_tri.add(new Triad<>(vertices.get(3), vertices.get(10), vertices.get(11)));
        face_norms_tri.add(normalTriple(vertices.get(3), vertices.get(10), vertices.get(11), true));

        faces_tri.add(new Triad<>(vertices.get(3), vertices.get(11), vertices.get(13)));
        face_norms_tri.add(normalTriple(vertices.get(3), vertices.get(11), vertices.get(13), true));

        faces_tri.add(new Triad<>(vertices.get(3), vertices.get(13), vertices.get(12)));
        face_norms_tri.add(normalTriple(vertices.get(3), vertices.get(13), vertices.get(12), true));

        faces_tri.add(new Triad<>(vertices.get(3), vertices.get(12), vertices.get(10)));
        face_norms_tri.add(normalTriple(vertices.get(3), vertices.get(12), vertices.get(10), true));

        faces_tri.add(new Triad<>(vertices.get(4), vertices.get(6), vertices.get(7)));
        face_norms_tri.add(normalTriple(vertices.get(4), vertices.get(6), vertices.get(7), true));

        faces_tri.add(new Triad<>(vertices.get(4), vertices.get(7), vertices.get(11)));
        face_norms_tri.add(normalTriple(vertices.get(4), vertices.get(7), vertices.get(11), true));

        faces_tri.add(new Triad<>(vertices.get(4), vertices.get(11), vertices.get(10)));
        face_norms_tri.add(normalTriple(vertices.get(4), vertices.get(11), vertices.get(10), true));

        faces_tri.add(new Triad<>(vertices.get(4), vertices.get(10), vertices.get(6)));
        face_norms_tri.add(normalTriple(vertices.get(4), vertices.get(10), vertices.get(6), true));

        faces_tri.add(new Triad<>(vertices.get(5), vertices.get(8), vertices.get(12)));
        face_norms_tri.add(normalTriple(vertices.get(5), vertices.get(8), vertices.get(12), true));

        faces_tri.add(new Triad<>(vertices.get(5), vertices.get(12), vertices.get(13)));
        face_norms_tri.add(normalTriple(vertices.get(5), vertices.get(12), vertices.get(13), true));

        faces_tri.add(new Triad<>(vertices.get(5), vertices.get(13), vertices.get(9)));
        face_norms_tri.add(normalTriple(vertices.get(5), vertices.get(13), vertices.get(9), true));

        faces_tri.add(new Triad<>(vertices.get(5), vertices.get(9), vertices.get(8)));
        face_norms_tri.add(normalTriple(vertices.get(5), vertices.get(9), vertices.get(8), true));
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < faces_tri.size(); i++) {

            // Triangular  faces
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.get(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.get(i);

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
