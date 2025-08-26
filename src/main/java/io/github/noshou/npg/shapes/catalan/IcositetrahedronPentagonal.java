package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.npg.shapes.archimedean.*;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b>Pentagonal Icositetrahedron</b>
 * <p> A Catalan solid with 38 vertices, 24 faces
 * (composed entirely of irregular pentagons), and 60 edges.
 * <p> It is the dual polyhedron of the {@link CuboctahedronSnub},
 * and therefore also exists in two chiral (mirror-image) enantiomorphs,
 * the right-handed {@link IcositetrahedronPentagonalDextro} and the
 * left-handed {@link IcositetrahedronPentagonalLevo}.
 */
@SuppressWarnings("FieldCanBeLocal")
public abstract class IcositetrahedronPentagonal extends Shape {


    // children must fill in these faces!

    // 24 Pentagonal Faces
    private ArrayList<Pentad<Tuple<Apfloat>>> faces_pnt;
    private final ArrayList<Triad<Apfloat>> face_norms_pnt= new ArrayList<>();

    // 38 vertices
    protected ArrayList<Triad<Apfloat>> vertices = new ArrayList<>();

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N6 = new Apfloat("6", super.precision);
    private final Apfloat N9 = new Apfloat("9", super.precision);
    private final Apfloat N12 = new Apfloat("12", super.precision);
    private final Apfloat N18 = new Apfloat("18", super.precision);
    private final Apfloat N33 = new Apfloat("33", super.precision);
    private final Apfloat N14 = new Apfloat("14", super.precision);
    private final Apfloat N1777 = new Apfloat("1777", super.precision);
    private final Apfloat SQRT33 = ApfloatMath.sqrt(N33);

    // cbrt(6*(9 + sqrt(33))
    private final Apfloat CBRT_6_mul_9_add_SQRT33 = ApfloatMath.cbrt(N6.multiply(N9.add(SQRT33)));

    // cbrt(6*(9 - sqrt(33))
    private final Apfloat CBRT_6_mul_9_sub_SQRT33 = ApfloatMath.cbrt(N6.multiply(N9.subtract(SQRT33)));

    // 33 * SQRT33
    private final Apfloat N33_mul_SQRT33 = N33.multiply(SQRT33);

    // cbrt(2 * (1777 + 33 * sqrt(33)))
    private final Apfloat CBRT_N2_mul_N1777_add_N33_mul_SQRT33 = ApfloatMath.cbrt(
            N2.multiply(N1777.add(N33_mul_SQRT33))
    );

    // cbrt(2 * (1777 - 33 * sqrt(33)))
    private final Apfloat CBRT_N2_mul_N1777_sub_N33_mul_SQRT33 = ApfloatMath.cbrt(
            N2.multiply(N1777.subtract(N33_mul_SQRT33))
    );

    /** @return 0 */
    public Apfloat N0() {return N0;}

    /**
     * @return sqrt(6 * (cbrt(6*(9 + sqrt(33))) + cbrt(6*(9 - sqrt(33))) - 6)) / 12
     */
    public Apfloat C0() {
        return ApfloatMath.sqrt(
                    N6.multiply(
                        CBRT_6_mul_9_add_SQRT33
                        .add(CBRT_6_mul_9_sub_SQRT33)
                        .subtract(N6)
                        )
                ).divide(N12);
    }

    /**
     * @return -(sqrt(6 * (cbrt(6*(9 + sqrt(33))) + cbrt(6*(9 - sqrt(33))) - 6)) / 12)
     */
    public Apfloat NEG_C0() {
        return C0().multiply(NEG_N1);
    }

    /**
     * @return sqrt(6 * (6 + cbrt(6*(9 + sqrt(33))) + cbrt(6*(9 - sqrt(33))))) / 12
     */
    public Apfloat C1() {
        return ApfloatMath.sqrt(
                    N6.multiply(
                        CBRT_6_mul_9_add_SQRT33
                        .add(CBRT_6_mul_9_sub_SQRT33)
                        .add(N6)
                    )
                ).divide(N12);
    }

    /**
     * @return -(sqrt(6 * (6 + cbrt(6*(9 + sqrt(33))) + cbrt(6*(9 - sqrt(33))))) / 12)
     */
    public Apfloat NEG_C1() {
        return C1().multiply(NEG_N1);
    }

    /**
     * @return sqrt(6 * (18 + cbrt(6*(9 + sqrt(33))) + cbrt(6*(9 - sqrt(33))))) / 12
     */
    public Apfloat C2() {
        return  ApfloatMath.sqrt(
                    N6.multiply(
                        CBRT_6_mul_9_add_SQRT33
                        .add(CBRT_6_mul_9_sub_SQRT33)
                        .add(N18))
                ).divide(N12);
    }

    /**
     * @return -(sqrt(6 * (18 + cbrt(6*(9 + sqrt(33))) + cbrt(6*(9 - sqrt(33))))) / 12)
     */
    public Apfloat NEG_C2() {
        return C2().multiply(NEG_N1);
    }

    /**
     * @return sqrt(6 * (14+cbrt(2*(1777+33*sqrt(33)))+cbrt(2*(1777-33*sqrt(33))))) / 12
     */
    public Apfloat C3() {
        return ApfloatMath.sqrt(
                    N6.multiply(
                        CBRT_N2_mul_N1777_add_N33_mul_SQRT33
                        .add(CBRT_N2_mul_N1777_sub_N33_mul_SQRT33)
                        .add(N14))
            ).divide(N12);
    }

    /**
     * @return -(sqrt(6 * (14+cbrt(2*(1777+33*sqrt(33)))+cbrt(2*(1777-33*sqrt(33))))) / 12)
     */
    public Apfloat NEG_C3() {
        return C3().multiply(NEG_N1);
    }
    // list of vertices
    private ArrayList<Triad<Apfloat>> vBase;

    // ==== BASIS VERTICES ====
    // MUST be filled by children
    protected Triad<Apfloat> vB0;
    protected Triad<Apfloat> vB1;
    protected Triad<Apfloat> vB2;
    protected Triad<Apfloat> vB3;
    protected Triad<Apfloat> vB4;
    protected Triad<Apfloat> vB5;
    protected Triad<Apfloat> vB6;
    protected Triad<Apfloat> vB7;
    protected Triad<Apfloat> vB8;
    protected Triad<Apfloat> vB9;
    protected Triad<Apfloat> vB10;
    protected Triad<Apfloat> vB11;
    protected Triad<Apfloat> vB12;
    protected Triad<Apfloat> vB13;
    protected Triad<Apfloat> vB14;
    protected Triad<Apfloat> vB15;
    protected Triad<Apfloat> vB16;
    protected Triad<Apfloat> vB17;
    protected Triad<Apfloat> vB18;
    protected Triad<Apfloat> vB19;
    protected Triad<Apfloat> vB20;
    protected Triad<Apfloat> vB21;
    protected Triad<Apfloat> vB22;
    protected Triad<Apfloat> vB23;
    protected Triad<Apfloat> vB24;
    protected Triad<Apfloat> vB25;
    protected Triad<Apfloat> vB26;
    protected Triad<Apfloat> vB27;
    protected Triad<Apfloat> vB28;
    protected Triad<Apfloat> vB29;
    protected Triad<Apfloat> vB30;
    protected Triad<Apfloat> vB31;
    protected Triad<Apfloat> vB32;
    protected Triad<Apfloat> vB33;
    protected Triad<Apfloat> vB34;
    protected Triad<Apfloat> vB35;
    protected Triad<Apfloat> vB36;
    protected Triad<Apfloat> vB37;

    /** fills scaled vertices */
    private void scaledVert(ArrayList<Triad<Apfloat>> basis) {
        for (Triad<Apfloat> vB_i: basis) {
            vertices.add(mult(normalize(vB_i), super.getRadius().toString()));
        }
    }

    /** fills faces_pnt and face_norms_pnt vertices */
    private void pntFaces(ArrayList<Pentad<Tuple<Apfloat>>> faces) {
        faces_pnt = faces;
        for (Pentad<Tuple<Apfloat>> pnt_i: faces) {
            face_norms_pnt.add(normalPent(
                            (Triad<Apfloat>) pnt_i.fetch(0),
                            (Triad<Apfloat>) pnt_i.fetch(1),
                            (Triad<Apfloat>) pnt_i.fetch(2),
                            (Triad<Apfloat>) pnt_i.fetch(3),
                            (Triad<Apfloat>) pnt_i.fetch(4),
                            true
                    )
            );
        }
    }

    /** Child classes must implement this method in their constructors */
    protected void setVerts() {
        vBase.add(vB0);
        vBase.add(vB1);
        vBase.add(vB2);
        vBase.add(vB3);
        vBase.add(vB4);
        vBase.add(vB5);
        vBase.add(vB6);
        vBase.add(vB7);
        vBase.add(vB8);
        vBase.add(vB9);
        vBase.add(vB10);
        vBase.add(vB11);
        vBase.add(vB12);
        vBase.add(vB13);
        vBase.add(vB14);
        vBase.add(vB16);
        vBase.add(vB17);
        vBase.add(vB18);
        vBase.add(vB19);
        vBase.add(vB20);
        vBase.add(vB21);
        vBase.add(vB22);
        vBase.add(vB23);
        vBase.add(vB24);
        vBase.add(vB25);
        vBase.add(vB26);
        vBase.add(vB27);
        vBase.add(vB28);
        vBase.add(vB29);
        vBase.add(vB30);
        vBase.add(vB31);
        vBase.add(vB32);
        vBase.add(vB33);
        vBase.add(vB34);
        vBase.add(vB35);
        vBase.add(vB36);
        vBase.add(vB37);

        // ==== SCALED VERTICES ====
        scaledVert(vBase);

        // ==== PENTAGONAL FACES ====
        ArrayList<Pentad<Tuple<Apfloat>>> pnt_faces = new ArrayList<>();
        pnt_faces.add(new Pentad<>(vertices.get(0),vertices.get(6),vertices.get(22),vertices.get(30),vertices.get(18)));
        pnt_faces.add(new Pentad<>(vertices.get(0),vertices.get(18),vertices.get(16),vertices.get(34),vertices.get(8)));
        pnt_faces.add(new Pentad<>(vertices.get(0),vertices.get(8),vertices.get(24),vertices.get(36),vertices.get(20)));
        pnt_faces.add(new Pentad<>(vertices.get(0),vertices.get(20),vertices.get(14),vertices.get(32),vertices.get(6)));
        pnt_faces.add(new Pentad<>(vertices.get(1),vertices.get(7),vertices.get(23),vertices.get(33),vertices.get(19)));
        pnt_faces.add(new Pentad<>(vertices.get(1),vertices.get(19),vertices.get(17),vertices.get(37),vertices.get(9)));
        pnt_faces.add(new Pentad<>(vertices.get(1),vertices.get(9),vertices.get(25),vertices.get(35),vertices.get(21)));
        pnt_faces.add(new Pentad<>(vertices.get(1),vertices.get(21),vertices.get(15),vertices.get(31),vertices.get(7)));
        pnt_faces.add(new Pentad<>(vertices.get(2),vertices.get(10),vertices.get(27),vertices.get(33),vertices.get(23)));
        pnt_faces.add(new Pentad<>(vertices.get(2),vertices.get(23),vertices.get(7),vertices.get(31),vertices.get(11)));
        pnt_faces.add(new Pentad<>(vertices.get(2),vertices.get(11),vertices.get(26),vertices.get(30),vertices.get(22)));
        pnt_faces.add(new Pentad<>(vertices.get(2),vertices.get(22),vertices.get(6),vertices.get(32),vertices.get(10)));
        pnt_faces.add(new Pentad<>(vertices.get(3),vertices.get(12),vertices.get(29),vertices.get(35),vertices.get(25)));
        pnt_faces.add(new Pentad<>(vertices.get(3),vertices.get(25),vertices.get(9),vertices.get(37),vertices.get(13)));
        pnt_faces.add(new Pentad<>(vertices.get(3),vertices.get(13),vertices.get(28),vertices.get(36),vertices.get(24)));
        pnt_faces.add(new Pentad<>(vertices.get(3),vertices.get(24),vertices.get(8),vertices.get(34),vertices.get(12)));
        pnt_faces.add(new Pentad<>(vertices.get(4),vertices.get(15),vertices.get(21),vertices.get(35),vertices.get(29)));
        pnt_faces.add(new Pentad<>(vertices.get(4),vertices.get(29),vertices.get(12),vertices.get(34),vertices.get(16)));
        pnt_faces.add(new Pentad<>(vertices.get(4),vertices.get(16),vertices.get(18),vertices.get(30),vertices.get(26)));
        pnt_faces.add(new Pentad<>(vertices.get(4),vertices.get(26),vertices.get(11),vertices.get(31),vertices.get(15)));
        pnt_faces.add(new Pentad<>(vertices.get(5),vertices.get(14),vertices.get(20),vertices.get(36),vertices.get(28)));
        pnt_faces.add(new Pentad<>(vertices.get(5),vertices.get(28),vertices.get(13),vertices.get(37),vertices.get(17)));
        pnt_faces.add(new Pentad<>(vertices.get(5),vertices.get(17),vertices.get(19),vertices.get(33),vertices.get(27)));
        pnt_faces.add(new Pentad<>(vertices.get(5),vertices.get(27),vertices.get(10),vertices.get(32),vertices.get(14)));
        pntFaces(pnt_faces);
    }
    public IcositetrahedronPentagonal(
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

    /**
     * Tests whether the given Cartesian point lies inside or on the boundary
     * <p>The test is performed by checking dot products against all
     * face normals (squares and triangles).
     * @param point_cart the Cartesian coordinate point
     * @return {@code true} if inside or on surface, {@code false} if outside
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_pnt.size(); i++) {


            // === Check pentagonal faces ===
            Pentad<Tuple<Apfloat>> face = (Pentad<Tuple<Apfloat>>) faces_pnt.get(i);
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
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_pnt.get(i);
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
