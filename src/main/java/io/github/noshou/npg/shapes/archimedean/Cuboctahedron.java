package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.npg.shapes.platonic.*;
import io.github.noshou.npg.shapes.catalan.*;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.*;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b>Cuboctahedron</b>
 * <p> An Archimedean solid with 12 vertices, 24 edges,
 * and 14 faces consisting of 8 equilateral triangles and 6 squares.
 * The Cuboctahedron can be constructed by truncating the vertices of a {@link Cube} or
 * {@link Octahedron}, resulting in a polyhedron where all vertices are equivalent, and all faces are regular polygons.
 * <p> It is the dual of the {@link DodecahedronRhombic}
 * @see <a href="https://dmccooey.com/polyhedra/Cuboctahedron.html">
 *     Cuboctahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class Cuboctahedron extends Shape {

    // the 6 square faces
    private final Hexad<Tuple<Tuple<Apfloat>>> faces_sqr;
    private final Hexad<Tuple<Apfloat>> face_norms_sqr;

    // the 8 triangular faces
    private final Octad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Octad<Tuple<Apfloat>> face_norms_tri;

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);

    public Cuboctahedron(
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

        // (±1,±1,0)
        Triad<Apfloat> vB0 = new Triad<>(N1, N1, N0);    // (+1, +1, 0)
        Triad<Apfloat> vB1 = new Triad<>(N1, NEG_N1, N0);    // (+1, -1, 0)
        Triad<Apfloat> vB2 = new Triad<>(NEG_N1, N1, N0);    // (-1, +1, 0)
        Triad<Apfloat> vB3 = new Triad<>(NEG_N1, NEG_N1, N0);    // (-1, -1, 0)

        // (±1,0,±1)
        Triad<Apfloat> vB4 = new Triad<>(N1, N0, N1);    // (+1, 0, +1)
        Triad<Apfloat> vB5 = new Triad<>(N1, N0, NEG_N1);    // (+1, 0, -1)
        Triad<Apfloat> vB6 = new Triad<>(NEG_N1, N0, N1);    // (-1, 0, +1)
        Triad<Apfloat> vB7 = new Triad<>(NEG_N1, N0, NEG_N1);    // (-1, 0, -1)

        // (0,±1,±1)
        Triad<Apfloat> vB8 = new Triad<>(N0, N1, N1);    // (0, +1, +1)
        Triad<Apfloat> vB9 = new Triad<>(N0, N1, NEG_N1);    // (0, +1, -1)
        Triad<Apfloat> vB10= new Triad<>(N0, NEG_N1, N1);    // (0, -1, +1)
        Triad<Apfloat> vB11= new Triad<>(N0, NEG_N1, NEG_N1);    // (0, -1, -1)

        // ==== SCALING FACTOR ====
        // dist from origin = sqrt(x²+y²+z²) = sqrt(1²+1²) = sqrt(2)
        Apfloat dist = ApfloatMath.sqrt(N1.add(N1));
        Apfloat r = super.getRadius();
        String r_dist = r.divide(dist).toString();

        // ==== SCALED CARTESIAN VERTICES ====

        // (±(r/sqrt(2)),±(r/sqrt(2)),0)
        Triad<Apfloat> vC0 = mult(vB0, r_dist);    // (+(r/sqrt(2)), +(r/sqrt(2)), 0)
        Triad<Apfloat> vC1 = mult(vB1, r_dist);    // (+(r/sqrt(2)), -(r/sqrt(2)), 0)
        Triad<Apfloat> vC2 = mult(vB2, r_dist);    // (-(r/sqrt(2)), +(r/sqrt(2)), 0)
        Triad<Apfloat> vC3 = mult(vB3, r_dist);    // (-(r/sqrt(2)), -(r/sqrt(2)), 0)

        // (±(r/sqrt(2)),0,±(r/sqrt(2)))
        Triad<Apfloat> vC4 = mult(vB4, r_dist);    // (+(r/sqrt(2)), 0, +(r/sqrt(2)))
        Triad<Apfloat> vC5 = mult(vB5, r_dist);    // (+(r/sqrt(2)), 0, -(r/sqrt(2)))
        Triad<Apfloat> vC6 = mult(vB6, r_dist);    // (-(r/sqrt(2)), 0, +(r/sqrt(2)))
        Triad<Apfloat> vC7 = mult(vB7, r_dist);     // (-(r/sqrt(2)), 0, -(r/sqrt(2)))

        // (0,±(r/sqrt(2)),±(r/sqrt(2)))
        Triad<Apfloat> vC8 = mult(vB8, r_dist);    // (0, +(r/sqrt(2)), +(r/sqrt(2)))
        Triad<Apfloat> vC9 = mult(vB9, r_dist);    // (0, +(r/sqrt(2)), -(r/sqrt(2)))
        Triad<Apfloat> vC10= mult(vB10,r_dist);    // (0, -(r/sqrt(2)), +(r/sqrt(2)))
        Triad<Apfloat> vC11= mult(vB11,r_dist);    // (0, -(r/sqrt(2)), -(r/sqrt(2)))

        // ==== QUADRILATERAL FACE DEFINITIONS ====

        // Square 1: X = +r/√2 plane (vertices with x-coordinate = +r/√2)
        // Contains: vC0(+,+,0), vC1(+,-,0), vC4(+,0,+), vC5(+,0,-)
        Tetrad<Tuple<Apfloat>> posX = new Tetrad<>(vC0, vC5, vC1, vC4);
        Triad<Apfloat> posX_norm = normalQuad(vC0, vC5, vC1, vC4, true);

        // Square 2: X = -r/√2 plane (vertices with x-coordinate = -r/√2)
        // Contains: vC2(-,+,0), vC3(-,-,0), vC6(-,0,+), vC7(-,0,-)
        Tetrad<Tuple<Apfloat>> negX = new Tetrad<>(vC2, vC7, vC3, vC6);
        Triad<Apfloat> negX_norm = normalQuad(vC2, vC7, vC3, vC6, true);

        // Square 3: Y = +r/√2 plane (vertices with y-coordinate = +r/√2)
        // Contains: vC0(+,+,0), vC2(-,+,0), vC8(0,+,+), vC9(0,+,-)
        Tetrad<Tuple<Apfloat>> posY = new Tetrad<>(vC0, vC9, vC2, vC8);
        Triad<Apfloat> posY_norm = normalQuad(vC0, vC9, vC2, vC8, true);

        // Square 4: Y = -r/√2 plane (vertices with y-coordinate = -r/√2)
        // Contains: vC1(+,-,0), vC3(-,-,0), vC10(0,-,+), vC11(0,-,-)
        Tetrad<Tuple<Apfloat>> negY = new Tetrad<>(vC1, vC10, vC3, vC11);
        Triad<Apfloat> negY_norm = normalQuad(vC1, vC10, vC3, vC11, true);

        // Square 5: Z = +r/√2 plane (vertices with z-coordinate = +r/√2)
        // Contains: vC4(+,0,+), vC6(-,0,+), vC8(0,+,+), vC10(0,-,+)
        Tetrad<Tuple<Apfloat>> posZ = new Tetrad<>(vC4, vC8, vC6, vC10);
        Triad<Apfloat> posZ_norm = normalQuad(vC4, vC8, vC6, vC10, true);

        // Square 6: Z = -r/√2 plane (vertices with z-coordinate = -r/√2)
        // Contains: vC5(+,0,-), vC7(-,0,-), vC9(0,+,-), vC11(0,-,-)
        Tetrad<Tuple<Apfloat>> negZ = new Tetrad<>(vC5, vC11, vC7, vC9);
        Triad<Apfloat> negZ_norm = normalQuad(vC5, vC11, vC7, vC9, true);

        faces_sqr = new Hexad<>(
                posX,
                negX,
                posY,
                negY,
                posZ,
                negZ
        );
        face_norms_sqr = new Hexad<>(
                posX_norm,
                negX_norm,
                posY_norm,
                negY_norm,
                posZ_norm,
                negZ_norm
        );

        // ==== TRIANGULAR FACE DEFINITIONS ====

        // Triangle 1: (+x,+y,+z) octant
        // Vertices: vC0(+,+,0), vC4(+,0,+), vC8(0,+,+)
        Triad<Tuple<Apfloat>> posX_posY_posZ = new Triad<>(vC0, vC4, vC8);
        Triad<Apfloat> posX_posY_posZ_norm = normalTriple(vC0, vC4, vC8, true);

        // Triangle 2: (+x,+y,-z) octant
        // Vertices: vC0(+,+,0), vC5(+,0,-), vC9(0,+,-)
        Triad<Tuple<Apfloat>> posX_posY_negZ = new Triad<>(vC0, vC9, vC5);
        Triad<Apfloat> posX_posY_negZ_norm = normalTriple(vC0, vC9, vC5, true);

        // Triangle 3: (-x,+y,-z) octant
        // Vertices: vC2(-,+,0), vC7(-,0,-), vC9(0,+,-)
        Triad<Tuple<Apfloat>> negX_posY_negZ = new Triad<>(vC2, vC9, vC7);
        Triad<Apfloat> negX_posY_negZ_norm = normalTriple(vC2, vC9, vC7, true);

        // Triangle 4: (+x,-y,-z) octant
        // Vertices: vC1(+,-,0), vC5(+,0,-), vC11(0,-,-)
        Triad<Tuple<Apfloat>> posX_negY_negZ = new Triad<>(vC1, vC11, vC5);
        Triad<Apfloat> posX_negY_negZ_norm = normalTriple(vC1, vC11, vC5, true);

        // Triangle 5: (-x,-y,-z) octant
        // Vertices: vC3(-,-,0), vC7(-,0,-), vC11(0,-,-)
        Triad<Tuple<Apfloat>> negX_negY_negZ = new Triad<>(vC3, vC11, vC7);
        Triad<Apfloat> negX_negY_negZ_norm = normalTriple(vC3, vC11, vC7, true);

        // Triangle 6: (-x,-y,+z) octant
        // Vertices: vC3(-,-,0), vC6(-,0,+), vC10(0,-,+)
        Triad<Tuple<Apfloat>> negX_negY_posZ = new Triad<>(vC3, vC10, vC6);
        Triad<Apfloat> negX_negY_posZ_norm = normalTriple(vC3, vC10, vC6, true);

        // Triangle 7: (-x,+y,+z) octant
        // Vertices: vC2(-,+,0), vC6(-,0,+), vC8(0,+,+)
        Triad<Tuple<Apfloat>> negX_posY_posZ = new Triad<>(vC2, vC8, vC6);
        Triad<Apfloat> negX_posY_posZ_norm = normalTriple(vC2, vC8, vC6, true);

        // Triangle 8: (+x,-y,+z) octant
        // Vertices: vC1(+,-,0), vC4(+,0,+), vC10(0,-,+)
        Triad<Tuple<Apfloat>> posX_negY_posZ = new Triad<>(vC1, vC10, vC4);
        Triad<Apfloat> posX_negY_posZ_norm = normalTriple(vC1, vC10, vC4, true);

        faces_tri = new Octad<>(
                posX_posY_posZ,
                posX_posY_negZ,
                negX_posY_negZ,
                posX_negY_negZ,
                negX_negY_negZ,
                negX_negY_posZ,
                negX_posY_posZ,
                posX_negY_posZ
                );

        face_norms_tri = new Octad<>(
                posX_posY_posZ_norm,
                posX_posY_negZ_norm,
                negX_posY_negZ_norm,
                posX_negY_negZ_norm,
                negX_negY_negZ_norm,
                negX_negY_posZ_norm,
                negX_posY_posZ_norm,
                posX_negY_posZ_norm
        );
    }

    /**
     * Tests whether the given Cartesian point lies inside or on the boundary
     * of the Cuboctahedron.
     *
     * <p>The test is performed by checking dot products against all
     * face normals (squares and triangles).
     *
     * @param point_cart the Cartesian coordinate point
     * @return {@code true} if inside or on surface, {@code false} if outside
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {


        for (int i = 0; i < faces_tri.fetchSize(); i++) {

            // === Check square faces ===
            if (i < faces_sqr.fetchSize()) {
                Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.fetch(i);
                Triad<Apfloat> vert0 = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> vert1 = (Triad<Apfloat>) face.fetch(1);
                Triad<Apfloat> vert2 = (Triad<Apfloat>) face.fetch(2);
                Triad<Apfloat> vert3 = (Triad<Apfloat>) face.fetch(3);

                // centroid of square
                Triad<Apfloat> centroid = new Triad<>(
                        vert0.fetch(0).add(vert1.fetch(0)).add(vert2.fetch(0)).add(vert3.fetch(0)).divide(N4),
                        vert0.fetch(1).add(vert1.fetch(1)).add(vert2.fetch(1)).add(vert3.fetch(1)).divide(N4),
                        vert0.fetch(2).add(vert1.fetch(2)).add(vert2.fetch(2)).add(vert3.fetch(2)).divide(N4)
                );

                // given point p and vertA, calculate vector from centroid -> p:
                // m = p - centroid = (p_x - centroid_x, p_x - centroid_y, p_x -centroid_z)
                Triad<Apfloat> m = subs(point_cart, centroid);

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

            // === Check triangular faces ===
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.fetch(i);
            Triad<Apfloat> vert0 = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> vert1 = (Triad<Apfloat>) face.fetch(1);
            Triad<Apfloat> vert2 = (Triad<Apfloat>) face.fetch(2);

            // centroid of triangle
            Triad<Apfloat> centroid = new Triad<>(
                    vert0.fetch(0).add(vert1.fetch(0)).add(vert2.fetch(0)).divide(N3),
                    vert0.fetch(1).add(vert1.fetch(1)).add(vert2.fetch(1)).divide(N3),
                    vert0.fetch(2).add(vert1.fetch(2)).add(vert2.fetch(2)).divide(N3)
            );

            // given point p and vertA, calculate vector from centroid -> p:
            // m = p - centroid = (p_x - centroid_x, p_x - centroid_y, p_x -centroid_z)
            Triad<Apfloat> m = subs(point_cart, centroid);

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
        return true;
    }
}
