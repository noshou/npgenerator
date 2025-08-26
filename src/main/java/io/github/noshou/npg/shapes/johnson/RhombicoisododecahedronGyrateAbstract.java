package io.github.noshou.npg.shapes.johnson;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.npg.shapes.archimedean.*;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * <p>Represents a family of <b>Gyrated Rhombicosidodecahedrons</b>.</p>
 * <p>These are variations of the {@link Rhombicosidodecahedron} in which one,
 * two, or three pentagonal cupola are rotated by 36 degrees.
 * <p> Depending on which cupola are rotated, the solid may be classified as gyrate
 * ({@link RhombicoisododecahedronGyrate}),
 * trigyrate ({@link RhombicoisododecahedronTriGyrate}),
 * paragbiyrate ({@link RhombicoisododecahedronParaBiGyrate}),
 * or metagbiyrate ({@link RhombicoisododecahedronMetaBiGyrate}).</p>
 */
public abstract class RhombicoisododecahedronGyrateAbstract extends Shape{

    // faces must be filled by children in their constructors

    // 30 square faces
    private ArrayList<Tetrad<Tuple<Apfloat>>> faces_sqr;
    private final ArrayList<Triad<Apfloat>> face_norms_sqr = new ArrayList<>();

    // 20 triangular faces
    private ArrayList<Triad<Tuple<Apfloat>>> faces_tri;
    private final ArrayList<Triad<Apfloat>> face_norms_tri = new ArrayList<>();

    // 12 pentagonal faces
    private ArrayList<Pentad<Tuple<Apfloat>>> faces_pnt;
    private final ArrayList<Triad<Apfloat>> face_norms_pnt = new ArrayList<>();

    // 60 vertices
    protected ArrayList<Triad<Apfloat>> vertices = new ArrayList<>();


    // Apfloat Constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N10 = new Apfloat("10", super.precision);
    private final Apfloat N13 = new Apfloat("13", super.precision);
    private final Apfloat N15 = new Apfloat("15", super.precision);
    private final Apfloat N20 = new Apfloat("20", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(N5);

    /** @return 0 */
    public Apfloat N0() {return N0;}

    /** @return ½ */
    public Apfloat HALF() {return N1.divide(N2);}

    /** @return -½ */
    public Apfloat NEG_HALF() {return NEG_N1.divide(N2);}

    // C0  = (5 + sqrt(5)) / 20
    private final Apfloat C0 = (N5.add(SQRT5)).divide(N20);
    /** @return (5 + sqrt(5)) / 20 */
    public Apfloat C0() {return C0;}
    /** @return -(5 + sqrt(5)) / 20 */
    public Apfloat NEG_C0() {return C0.multiply(NEG_N1);}

    // C1  = (1 + sqrt(5)) / 4
    private final Apfloat C1 = (N1.add(SQRT5)).divide(N4);
    /** @return (1 + sqrt(5)) / 4 */
    public Apfloat C1() {return C1;}
    /** @return -(1 + sqrt(5)) / 4 */
    public Apfloat NEG_C1() {return C1.multiply(NEG_N1);}

    // C2  = (15 + sqrt(5)) / 20
    private final Apfloat C2 = (N15.add(SQRT5)).divide(N20);
    /** @return (15 + sqrt(5)) / 20 */
    public Apfloat C2() {return C2;}
    /** @return -(15 + sqrt(5)) / 20 */
    public Apfloat NEG_C2() {return C2.multiply(NEG_N1);}

    // C3  = (3 + sqrt(5)) / 4
    private final Apfloat C3 = (N3.add(SQRT5)).divide(N4);
    /** @return (3 + sqrt(5)) / 4 */
    public Apfloat C3() {return C3;}
    /** @return -(3 + sqrt(5)) / 4 */
    public Apfloat NEG_C3() {return C3.multiply(NEG_N1);}

    // C4  = (5 + 4 * sqrt(5)) / 10
    private final Apfloat C4 = (N5.add(N4.multiply(SQRT5))).divide(N10);
    /** @return (5 + 4 * sqrt(5)) / 10 */
    public Apfloat C4() {return C4;}
    /** @return -(5 + 4 * sqrt(5)) / 10 */
    public Apfloat NEG_C4() {return C4.multiply(NEG_N1);}

    // C5  = (1 + sqrt(5)) / 2
    private final Apfloat C5 = (N1.add(SQRT5)).divide(N2);
    /** @return (1 + sqrt(5)) / 2 */
    public Apfloat C5() {return C5;}
    /** @return -(1 + sqrt(5)) / 2 */
    public Apfloat NEG_C5() {return C5.multiply(NEG_N1);}

    // C6  = (10 + 3 * sqrt(5)) / 10
    private final Apfloat C6 = (N10.add(N3.multiply(SQRT5))).divide(N10);
    /** @return (10 + 3 * sqrt(5)) / 10 */
    public Apfloat C6() {return C6;}
    /** @return -(10 + 3 * sqrt(5)) / 10 */
    public Apfloat NEG_C6() {return C6.multiply(NEG_N1);}

    // C7  = (5 + sqrt(5)) / 4
    private final Apfloat C7 = (N5.add(SQRT5)).divide(N4);
    /** @return (5 + sqrt(5)) / 4 */
    public Apfloat C7() {return C7;}
    /** @return -(5 + sqrt(5)) / 4 */
    public Apfloat NEG_C7() {return C7.multiply(NEG_N1);}

    // C8  = (5 + 2 * sqrt(5)) / 5
    private final Apfloat C8 = (N5.add(N2.multiply(SQRT5))).divide(N5);
    /** @return (5 + 2 * sqrt(5)) / 5 */
    public Apfloat C8() {return C8;}
    /** @return -(5 + 2 * sqrt(5)) / 5 */
    public Apfloat NEG_C8() {return C8.multiply(NEG_N1);}

    // C9  = (2 + sqrt(5)) / 2
    private final Apfloat C9 = (N2.add(SQRT5)).divide(N2);
    /** @return (2 + sqrt(5)) / 2 */
    public Apfloat C9() {return C9;}
    /** @return -(2 + sqrt(5)) / 2 */
    public Apfloat NEG_C9() {return C9.multiply(NEG_N1);}

    // C10 = (15 + 13 * sqrt(5)) / 20
    private final Apfloat C10 = (N15.add(N13.multiply(SQRT5))).divide(N20);
    /** @return (15 + 13 * sqrt(5)) / 20 */
    public Apfloat C10() {return C10;}
    /** @return -(15 + 13 * sqrt(5)) / 20 */
    public Apfloat NEG_C10() {return C10.multiply(NEG_N1);}

    // list of vertices
    private ArrayList<Triad<Apfloat>> vBase;

    // ==== BASIS VERTICES ====
    // MUST be filled by children in constructors
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
    protected Triad<Apfloat> vB38;
    protected Triad<Apfloat> vB39;
    protected Triad<Apfloat> vB40;
    protected Triad<Apfloat> vB41;
    protected Triad<Apfloat> vB42;
    protected Triad<Apfloat> vB43;
    protected Triad<Apfloat> vB44;
    protected Triad<Apfloat> vB45;
    protected Triad<Apfloat> vB46;
    protected Triad<Apfloat> vB47;
    protected Triad<Apfloat> vB48;
    protected Triad<Apfloat> vB49;
    protected Triad<Apfloat> vB50;
    protected Triad<Apfloat> vB51;
    protected Triad<Apfloat> vB52;
    protected Triad<Apfloat> vB53;
    protected Triad<Apfloat> vB54;
    protected Triad<Apfloat> vB55;
    protected Triad<Apfloat> vB56;
    protected Triad<Apfloat> vB57;
    protected Triad<Apfloat> vB58;
    protected Triad<Apfloat> vB59;

    /** fills scaled vertices. must be called in child's constructor. */
    protected void scaledVert() {
        vertices.add(mult(normalize(vB0), super.getRadius().toString()));
        vertices.add(mult(normalize(vB1), super.getRadius().toString()));
        vertices.add(mult(normalize(vB2), super.getRadius().toString()));
        vertices.add(mult(normalize(vB3), super.getRadius().toString()));
        vertices.add(mult(normalize(vB4), super.getRadius().toString()));
        vertices.add(mult(normalize(vB5), super.getRadius().toString()));
        vertices.add(mult(normalize(vB6), super.getRadius().toString()));
        vertices.add(mult(normalize(vB7), super.getRadius().toString()));
        vertices.add(mult(normalize(vB8), super.getRadius().toString()));
        vertices.add(mult(normalize(vB9), super.getRadius().toString()));
        vertices.add(mult(normalize(vB10), super.getRadius().toString()));
        vertices.add(mult(normalize(vB11), super.getRadius().toString()));
        vertices.add(mult(normalize(vB12), super.getRadius().toString()));
        vertices.add(mult(normalize(vB13), super.getRadius().toString()));
        vertices.add(mult(normalize(vB14), super.getRadius().toString()));
        vertices.add(mult(normalize(vB15), super.getRadius().toString()));
        vertices.add(mult(normalize(vB16), super.getRadius().toString()));
        vertices.add(mult(normalize(vB17), super.getRadius().toString()));
        vertices.add(mult(normalize(vB18), super.getRadius().toString()));
        vertices.add(mult(normalize(vB19), super.getRadius().toString()));
        vertices.add(mult(normalize(vB20), super.getRadius().toString()));
        vertices.add(mult(normalize(vB21), super.getRadius().toString()));
        vertices.add(mult(normalize(vB22), super.getRadius().toString()));
        vertices.add(mult(normalize(vB23), super.getRadius().toString()));
        vertices.add(mult(normalize(vB24), super.getRadius().toString()));
        vertices.add(mult(normalize(vB25), super.getRadius().toString()));
        vertices.add(mult(normalize(vB26), super.getRadius().toString()));
        vertices.add(mult(normalize(vB27), super.getRadius().toString()));
        vertices.add(mult(normalize(vB28), super.getRadius().toString()));
        vertices.add(mult(normalize(vB29), super.getRadius().toString()));
        vertices.add(mult(normalize(vB30), super.getRadius().toString()));
        vertices.add(mult(normalize(vB31), super.getRadius().toString()));
        vertices.add(mult(normalize(vB32), super.getRadius().toString()));
        vertices.add(mult(normalize(vB33), super.getRadius().toString()));
        vertices.add(mult(normalize(vB34), super.getRadius().toString()));
        vertices.add(mult(normalize(vB35), super.getRadius().toString()));
        vertices.add(mult(normalize(vB36), super.getRadius().toString()));
        vertices.add(mult(normalize(vB37), super.getRadius().toString()));
        vertices.add(mult(normalize(vB38), super.getRadius().toString()));
        vertices.add(mult(normalize(vB39), super.getRadius().toString()));
        vertices.add(mult(normalize(vB40), super.getRadius().toString()));
        vertices.add(mult(normalize(vB41), super.getRadius().toString()));
        vertices.add(mult(normalize(vB42), super.getRadius().toString()));
        vertices.add(mult(normalize(vB43), super.getRadius().toString()));
        vertices.add(mult(normalize(vB44), super.getRadius().toString()));
        vertices.add(mult(normalize(vB45), super.getRadius().toString()));
        vertices.add(mult(normalize(vB46), super.getRadius().toString()));
        vertices.add(mult(normalize(vB47), super.getRadius().toString()));
        vertices.add(mult(normalize(vB48), super.getRadius().toString()));
        vertices.add(mult(normalize(vB49), super.getRadius().toString()));
        vertices.add(mult(normalize(vB50), super.getRadius().toString()));
        vertices.add(mult(normalize(vB51), super.getRadius().toString()));
        vertices.add(mult(normalize(vB52), super.getRadius().toString()));
        vertices.add(mult(normalize(vB53), super.getRadius().toString()));
        vertices.add(mult(normalize(vB54), super.getRadius().toString()));
        vertices.add(mult(normalize(vB55), super.getRadius().toString()));
        vertices.add(mult(normalize(vB56), super.getRadius().toString()));
        vertices.add(mult(normalize(vB57), super.getRadius().toString()));
        vertices.add(mult(normalize(vB58), super.getRadius().toString()));
        vertices.add(mult(normalize(vB59), super.getRadius().toString()));
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

    /** fills faces_sqr and face_norms_sqr vertices */
    private void sqrFaces(ArrayList<Tetrad<Tuple<Apfloat>>> faces) {
        faces_sqr = faces;
        for (Tetrad<Tuple<Apfloat>> sqr_i: faces) {
            face_norms_sqr.add(normalQuad(
                            (Triad<Apfloat>) sqr_i.fetch(0),
                            (Triad<Apfloat>) sqr_i.fetch(1),
                            (Triad<Apfloat>) sqr_i.fetch(2),
                            (Triad<Apfloat>) sqr_i.fetch(3),
                            true
                    )
            );
        }
    }

    /** fills faces_tri and face_norms_tri vertices */
    private void triFaces(ArrayList<Triad<Tuple<Apfloat>>> faces) {
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

    /** fills faces and face_norms vertices */
    protected void setFaces(
            ArrayList<Pentad<Tuple<Apfloat>>> faces_pnt,
            ArrayList<Tetrad<Tuple<Apfloat>>> faces_sqr,
            ArrayList<Triad<Tuple<Apfloat>>> faces_tri
    ) {
        pntFaces(faces_pnt);
        sqrFaces(faces_sqr);
        triFaces(faces_tri);
    }

    public RhombicoisododecahedronGyrateAbstract (
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

        for (int i = 0; i < 30; i++) {

            // triangular faces
            if (i < 12){
                Triad<Tuple<Apfloat>> face = faces_tri.get(i);
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

            // pentagonal faces
            if (i < 12){
                Pentad<Tuple<Apfloat>> face = faces_pnt.get(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_pnt.get(i);

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
            // Square  faces
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.get(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_sqr.get(i);

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
