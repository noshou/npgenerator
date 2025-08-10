package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.*;
import static io.github.noshou.npg.nputil.VectorMath.*;
import java.util.ArrayList;

/**
 * Represents a <b>Snub Cuboctahedron (aka: Snub Cube)</b>
 * <p>
 * An Archimedean solid with 60 vertices, 92 faces
 * (comprised of 32 triangles and 6 squares), and 150 edges.
 * <p>
 * It is a chiral polyhedron, meaning it exists in two mirror-image forms:
 * the <i>levo</i> (left-handed) and <i>dextro</i> (right-handed) variants.
 * These enantiomorphs are not superimposable on each other.
 *
 */
@SuppressWarnings("FieldCanBeLocal")
public abstract class CuboctahedronSnub extends Shape {

    // children must fill in these faces!

    // 32 Equilateral triangles
    private ArrayList<Triad<Tuple<Apfloat>>> faces_tri = new ArrayList<>();
    private final ArrayList<Triad<Apfloat>> face_norms_tri= new ArrayList<>();

    // the 6 square faces
    private ArrayList<Tetrad<Tuple<Apfloat>>> faces_sqr = new ArrayList<>();
    private final ArrayList<Triad<Apfloat>> face_norms_sqr = new ArrayList<>();

    // 24 vertices
    protected ArrayList<Triad<Apfloat>> vertices = new ArrayList<>();

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N6 = new Apfloat("6", super.precision);
    private final Apfloat N17 = new Apfloat("17", super.precision);
    private final Apfloat N33 = new Apfloat("33", super.precision);
    private final Apfloat N199 = new Apfloat("199", super.precision);
    private final Apfloat N3mulSQRT33 = ApfloatMath.sqrt(N3.multiply(N33));

    /**
     * C0 = sqrt(3 * (4 - cbrt(17 + 3*sqrt(33)) - cbrt(17 - 3*sqrt(33)))) / 6
     */
    public Apfloat C0() {
        return ApfloatMath.sqrt(
                N3.multiply(N4.subtract(ApfloatMath.cbrt(N17.add(N3mulSQRT33))))
                .subtract(ApfloatMath.cbrt(N17.subtract(N3mulSQRT33)))
                ).divide(N6);
    }

    /**
     * NEG_C0 = -(sqrt(3 * (4 - cbrt(17 + 3*sqrt(33)) - cbrt(17 - 3*sqrt(33)))) / 6)
     */
    public Apfloat NEG_C0() {
        return C0().multiply(NEG_N1);
    }

    /**
     * C1 =sqrt(3 * (2 + cbrt(17 + 3*sqrt(33)) + cbrt(17 - 3*sqrt(33)))) / 6
     */
    public Apfloat C1() {
        return ApfloatMath.sqrt(
                N3.multiply(N2.add(ApfloatMath.cbrt(N17.add(N3mulSQRT33))))
                .add(ApfloatMath.cbrt(N17.subtract(N3mulSQRT33)))
                ).divide(N6);
    }

    /**
     * NEG_C1 = -(sqrt(3 * (2 + cbrt(17 + 3*sqrt(33)) + cbrt(17 - 3*sqrt(33)))) / 6)
     */
    public Apfloat NEG_C1() {
        return C1().multiply(NEG_N1);
    }

    /*
     * C2 = sqrt(3 * (4 + cbrt(199 + 3*sqrt(33)) + cbrt(199 - 3*sqrt(33)))) / 6
     */
    public Apfloat C2(){
        return ApfloatMath.sqrt(
                N3.multiply(N4.add(ApfloatMath.cbrt(N199.add(N3mulSQRT33))))
                .add(ApfloatMath.cbrt(N199.subtract(N3mulSQRT33)))
                ).divide(N6);
    }

    /*
     * NEG_C2 = -(sqrt(3 * (4 + cbrt(199 + 3*sqrt(33)) + cbrt(199 - 3*sqrt(33)))) / 6)
     */
    public Apfloat NEG_C2() {
        return C2().multiply(NEG_N1);
    }

    public CuboctahedronSnub(
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

    /** fills scaled vertices */
    protected void scaledVert(ArrayList<Triad<Apfloat>> basis) {
        for (Triad<Apfloat> vB_i: basis) {
            vertices.add(mult(normalize(vB_i), super.getRadius().toString()));
        }
    }

    /** fills faces_tri and face_norms_tri vertices */
    protected void triFaces(ArrayList<Triad<Tuple<Apfloat>>> faces) {
        faces_tri = faces;
        for (Triad<Tuple<Apfloat>> tri_i: faces) {
            face_norms_tri.add(normalTriple(
                    (Triad<Apfloat>) tri_i.fetch(0),
                    (Triad<Apfloat>) tri_i.fetch(1),
                    (Triad<Apfloat>) tri_i.fetch(2),
                    true
                    )
            );
        }
    }

    /** fills faces_sqr and face_norms_sqr vertices */
    protected void sqrFaces(ArrayList<Tetrad<Tuple<Apfloat>>> faces) {
        faces_sqr = faces;
        for (Tetrad<Tuple<Apfloat>> tri_i: faces) {
            face_norms_sqr.add(normalQuad(
                    (Triad<Apfloat>) tri_i.fetch(0),
                    (Triad<Apfloat>) tri_i.fetch(1),
                    (Triad<Apfloat>) tri_i.fetch(2),
                    (Triad<Apfloat>) tri_i.fetch(3),
                    true
                    )
            );
        }
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
        for (int i = 0; i < faces_tri.size(); i++) {

            // === Check square faces ===
            if (i < faces_sqr.size()) {
                Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.get(i);
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
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_sqr.get(i);
                Apfloat d = dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check triangular faces ===
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.get(i);
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
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.get(i);
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
