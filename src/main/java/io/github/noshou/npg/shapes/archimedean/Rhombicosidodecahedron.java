package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;
import static io.github.noshou.npg.nputil.VectorMath.*;
import io.github.noshou.npg.shapes.catalan.*;

/**
 * Represents a <b>Rhombicosidodecahedron</b>.
 * <p> An Archimedean solid composed of 20 regular triangles, 30 squares,
 * and 12 regular pentagons, for a total of 62 faces, 120 edges, and 60 vertices.
 * <p> It is the dual of the {@link HexecontahedronDeltoidal}.
 * @see <a href="https://dmccooey.com/polyhedra/Rhombicosidodecahedron.html">
 *      Rhombicosidodecahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class Rhombicosidodecahedron extends Shape {


    // 30 square faces
    private final ArrayList<Tetrad<Tuple<Apfloat>>> faces_sqr;
    private final ArrayList<Triad<Apfloat>> face_norms_sqr;

    // 20 triangular faces
    private final ArrayList<Triad<Tuple<Apfloat>>> faces_tri;
    private final ArrayList<Triad<Apfloat>> face_norms_tri;

    // 12 pentagonal faces
    private final ArrayList<Pentad<Tuple<Apfloat>>> faces_pnt;
    private final ArrayList<Triad<Apfloat>> face_norms_pnt;

    // Apfloat Constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(N5);
    private final Apfloat FRAC_SQRT5_over_N2 = SQRT5.divide(N2);
    private final Apfloat FRAC_SQRT5_over_N4 = SQRT5.divide(N4);
    private final Apfloat HALF = N1.divide(N2);

    // C0 = (1 + sqrt(5)) / 4
    private final Apfloat C0 = N1.divide(N4).add(FRAC_SQRT5_over_N4);

    // C1 = (3 + sqrt(5)) / 4
    private final Apfloat C1 = N3.divide(N4).add(FRAC_SQRT5_over_N4);

    // C2 = (1 + sqrt(5)) / 2
    private final Apfloat C2 = HALF.add(FRAC_SQRT5_over_N2);

    // C3 = (5 + sqrt(5)) / 4
    private final Apfloat C3 = N5.divide(N4).add(FRAC_SQRT5_over_N4);

    // C4 = (2 + sqrt(5)) / 2
    private final Apfloat C4 = N1.add(FRAC_SQRT5_over_N2);

    public Rhombicosidodecahedron(
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
        Triad<Apfloat> vB0 = new Triad<>(HALF, HALF, C4);
        Triad<Apfloat> vB1 = new Triad<>(HALF, HALF, C4.multiply(NEG_N1));
        Triad<Apfloat> vB2 = new Triad<>(HALF, HALF.multiply(NEG_N1), C4);
        Triad<Apfloat> vB3  = new Triad<>(HALF, HALF.multiply(NEG_N1), C4.multiply(NEG_N1));
        Triad<Apfloat> vB4  = new Triad<>(HALF.multiply(NEG_N1), HALF, C4);
        Triad<Apfloat> vB5  = new Triad<>(HALF.multiply(NEG_N1), HALF, C4.multiply(NEG_N1));
        Triad<Apfloat> vB6  = new Triad<>(HALF.multiply(NEG_N1), HALF.multiply(NEG_N1), C4);
        Triad<Apfloat> vB7  = new Triad<>(HALF.multiply(NEG_N1), HALF.multiply(NEG_N1), C4.multiply(NEG_N1));
        Triad<Apfloat> vB8  = new Triad<>(C4, HALF, HALF);
        Triad<Apfloat> vB9  = new Triad<>(C4, HALF, HALF.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(C4, HALF.multiply(NEG_N1), HALF);
        Triad<Apfloat> vB11 = new Triad<>(C4, HALF.multiply(NEG_N1), HALF.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(C4.multiply(NEG_N1), HALF, HALF);
        Triad<Apfloat> vB13 = new Triad<>(C4.multiply(NEG_N1), HALF, HALF.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(C4.multiply(NEG_N1), HALF.multiply(NEG_N1), HALF);
        Triad<Apfloat> vB15 = new Triad<>(C4.multiply(NEG_N1), HALF.multiply(NEG_N1), HALF.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>(HALF, C4, HALF);
        Triad<Apfloat> vB17 = new Triad<>(HALF, C4, HALF.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>(HALF, C4.multiply(NEG_N1), HALF);
        Triad<Apfloat> vB19 = new Triad<>(HALF, C4.multiply(NEG_N1), HALF.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>(HALF.multiply(NEG_N1), C4, HALF);
        Triad<Apfloat> vB21 = new Triad<>(HALF.multiply(NEG_N1), C4, HALF.multiply(NEG_N1));
        Triad<Apfloat> vB22 = new Triad<>(HALF.multiply(NEG_N1), C4.multiply(NEG_N1), HALF);
        Triad<Apfloat> vB23 = new Triad<>(HALF.multiply(NEG_N1), C4.multiply(NEG_N1), HALF.multiply(NEG_N1));
        Triad<Apfloat> vB24 = new Triad<>(N0, C1, C3);
        Triad<Apfloat> vB25 = new Triad<>(N0, C1, C3.multiply(NEG_N1));
        Triad<Apfloat> vB26 = new Triad<>(N0, C1.multiply(NEG_N1), C3);
        Triad<Apfloat> vB27 = new Triad<>(N0, C1.multiply(NEG_N1), C3.multiply(NEG_N1));
        Triad<Apfloat> vB28 = new Triad<>(C3, N0, C1);
        Triad<Apfloat> vB29 = new Triad<>(C3, N0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB30 = new Triad<>(C3.multiply(NEG_N1), N0, C1);
        Triad<Apfloat> vB31 = new Triad<>(C3.multiply(NEG_N1), N0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB32 = new Triad<>(C1, C3, N0);
        Triad<Apfloat> vB33 = new Triad<>(C1, C3.multiply(NEG_N1), N0);
        Triad<Apfloat> vB34 = new Triad<>(C1.multiply(NEG_N1), C3, N0);
        Triad<Apfloat> vB35 = new Triad<>(C1.multiply(NEG_N1), C3.multiply(NEG_N1), N0);
        Triad<Apfloat> vB36 = new Triad<>(C1, C0, C2);
        Triad<Apfloat> vB37 = new Triad<>(C1, C0, C2.multiply(NEG_N1));
        Triad<Apfloat> vB38 = new Triad<>(C1, C0.multiply(NEG_N1), C2);
        Triad<Apfloat> vB39 = new Triad<>(C1, C0.multiply(NEG_N1), C2.multiply(NEG_N1));
        Triad<Apfloat> vB40 = new Triad<>(C1.multiply(NEG_N1), C0, C2);
        Triad<Apfloat> vB41 = new Triad<>(C1.multiply(NEG_N1), C0, C2.multiply(NEG_N1));
        Triad<Apfloat> vB42 = new Triad<>(C1.multiply(NEG_N1), C0.multiply(NEG_N1), C2);
        Triad<Apfloat> vB43 = new Triad<>(C1.multiply(NEG_N1), C0.multiply(NEG_N1), C2.multiply(NEG_N1));
        Triad<Apfloat> vB44 = new Triad<>(C2, C1, C0);
        Triad<Apfloat> vB45 = new Triad<>(C2, C1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB46 = new Triad<>(C2, C1.multiply(NEG_N1), C0);
        Triad<Apfloat> vB47 = new Triad<>(C2, C1.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB48 = new Triad<>(C2.multiply(NEG_N1), C1, C0);
        Triad<Apfloat> vB49 = new Triad<>(C2.multiply(NEG_N1), C1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB50 = new Triad<>(C2.multiply(NEG_N1), C1.multiply(NEG_N1), C0);
        Triad<Apfloat> vB51 = new Triad<>(C2.multiply(NEG_N1), C1.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB52 = new Triad<>(C0, C2, C1);
        Triad<Apfloat> vB53 = new Triad<>(C0, C2, C1.multiply(NEG_N1));
        Triad<Apfloat> vB54 = new Triad<>(C0, C2.multiply(NEG_N1), C1);
        Triad<Apfloat> vB55 = new Triad<>(C0, C2.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB56 = new Triad<>(C0.multiply(NEG_N1), C2, C1);
        Triad<Apfloat> vB57 = new Triad<>(C0.multiply(NEG_N1), C2, C1.multiply(NEG_N1));
        Triad<Apfloat> vB58 = new Triad<>(C0.multiply(NEG_N1), C2.multiply(NEG_N1), C1);
        Triad<Apfloat> vB59 = new Triad<>(C0.multiply(NEG_N1), C2.multiply(NEG_N1), C1.multiply(NEG_N1));

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0 = mult(normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1 = mult(normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2 = mult(normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3 = mult(normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4 = mult(normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5 = mult(normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6 = mult(normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7 = mult(normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8 = mult(normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9 = mult(normalize(vB9), super.getRadius().toString());
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
        Triad<Apfloat> vC24 = mult(normalize(vB24), super.getRadius().toString());
        Triad<Apfloat> vC25 = mult(normalize(vB25), super.getRadius().toString());
        Triad<Apfloat> vC26 = mult(normalize(vB26), super.getRadius().toString());
        Triad<Apfloat> vC27 = mult(normalize(vB27), super.getRadius().toString());
        Triad<Apfloat> vC28 = mult(normalize(vB28), super.getRadius().toString());
        Triad<Apfloat> vC29 = mult(normalize(vB29), super.getRadius().toString());
        Triad<Apfloat> vC30 = mult(normalize(vB30), super.getRadius().toString());
        Triad<Apfloat> vC31 = mult(normalize(vB31), super.getRadius().toString());
        Triad<Apfloat> vC32 = mult(normalize(vB32), super.getRadius().toString());
        Triad<Apfloat> vC33 = mult(normalize(vB33), super.getRadius().toString());
        Triad<Apfloat> vC34 = mult(normalize(vB34), super.getRadius().toString());
        Triad<Apfloat> vC35 = mult(normalize(vB35), super.getRadius().toString());
        Triad<Apfloat> vC36 = mult(normalize(vB36), super.getRadius().toString());
        Triad<Apfloat> vC37 = mult(normalize(vB37), super.getRadius().toString());
        Triad<Apfloat> vC38 = mult(normalize(vB38), super.getRadius().toString());
        Triad<Apfloat> vC39 = mult(normalize(vB39), super.getRadius().toString());
        Triad<Apfloat> vC40 = mult(normalize(vB40), super.getRadius().toString());
        Triad<Apfloat> vC41 = mult(normalize(vB41), super.getRadius().toString());
        Triad<Apfloat> vC42 = mult(normalize(vB42), super.getRadius().toString());
        Triad<Apfloat> vC43 = mult(normalize(vB43), super.getRadius().toString());
        Triad<Apfloat> vC44 = mult(normalize(vB44), super.getRadius().toString());
        Triad<Apfloat> vC45 = mult(normalize(vB45), super.getRadius().toString());
        Triad<Apfloat> vC46 = mult(normalize(vB46), super.getRadius().toString());
        Triad<Apfloat> vC47 = mult(normalize(vB47), super.getRadius().toString());
        Triad<Apfloat> vC48 = mult(normalize(vB48), super.getRadius().toString());
        Triad<Apfloat> vC49 = mult(normalize(vB49), super.getRadius().toString());
        Triad<Apfloat> vC50 = mult(normalize(vB50), super.getRadius().toString());
        Triad<Apfloat> vC51 = mult(normalize(vB51), super.getRadius().toString());
        Triad<Apfloat> vC52 = mult(normalize(vB52), super.getRadius().toString());
        Triad<Apfloat> vC53 = mult(normalize(vB53), super.getRadius().toString());
        Triad<Apfloat> vC54 = mult(normalize(vB54), super.getRadius().toString());
        Triad<Apfloat> vC55 = mult(normalize(vB55), super.getRadius().toString());
        Triad<Apfloat> vC56 = mult(normalize(vB56), super.getRadius().toString());
        Triad<Apfloat> vC57 = mult(normalize(vB57), super.getRadius().toString());
        Triad<Apfloat> vC58 = mult(normalize(vB58), super.getRadius().toString());
        Triad<Apfloat> vC59 = mult(normalize(vB59), super.getRadius().toString());

        // ==== SQUARE FACES ====
        faces_sqr = new ArrayList<>();
        face_norms_sqr = new ArrayList<>();

        faces_sqr.add(new Tetrad<>(vC0, vC36, vC52, vC24));
        face_norms_sqr.add(normalQuad(vC0, vC36, vC52, vC24, true));

        faces_sqr.add(new Tetrad<>(vC1, vC25, vC53, vC37));
        face_norms_sqr.add(normalQuad(vC1, vC25, vC53, vC37, true));

        faces_sqr.add(new Tetrad<>(vC2, vC26, vC54, vC38));
        face_norms_sqr.add(normalQuad(vC2, vC26, vC54, vC38, true));

        faces_sqr.add(new Tetrad<>(vC3, vC39, vC55, vC27));
        face_norms_sqr.add(normalQuad(vC3, vC39, vC55, vC27, true));

        faces_sqr.add(new Tetrad<>(vC4, vC24, vC56, vC40));
        face_norms_sqr.add(normalQuad(vC4, vC24, vC56, vC40, true));

        faces_sqr.add(new Tetrad<>(vC5, vC41, vC57, vC25));
        face_norms_sqr.add(normalQuad(vC5, vC41, vC57, vC25, true));

        faces_sqr.add(new Tetrad<>(vC6, vC42, vC58, vC26));
        face_norms_sqr.add(normalQuad(vC6, vC42, vC58, vC26, true));

        faces_sqr.add(new Tetrad<>(vC7, vC27, vC59, vC43));
        face_norms_sqr.add(normalQuad(vC7, vC27, vC59, vC43, true));

        faces_sqr.add(new Tetrad<>(vC8, vC44, vC36, vC28));
        face_norms_sqr.add(normalQuad(vC8, vC44, vC36, vC28, true));

        faces_sqr.add(new Tetrad<>(vC9, vC29, vC37, vC45));
        face_norms_sqr.add(normalQuad(vC9, vC29, vC37, vC45, true));

        faces_sqr.add(new Tetrad<>(vC10, vC28, vC38, vC46));
        face_norms_sqr.add(normalQuad(vC10, vC28, vC38, vC46, true));

        faces_sqr.add(new Tetrad<>(vC11, vC47, vC39, vC29));
        face_norms_sqr.add(normalQuad(vC11, vC47, vC39, vC29, true));

        faces_sqr.add(new Tetrad<>(vC12, vC30, vC40, vC48));
        face_norms_sqr.add(normalQuad(vC12, vC30, vC40, vC48, true));

        faces_sqr.add(new Tetrad<>(vC13, vC49, vC41, vC31));
        face_norms_sqr.add(normalQuad(vC13, vC49, vC41, vC31, true));

        faces_sqr.add(new Tetrad<>(vC14, vC50, vC42, vC30));
        face_norms_sqr.add(normalQuad(vC14, vC50, vC42, vC30, true));

        faces_sqr.add(new Tetrad<>(vC15, vC31, vC43, vC51));
        face_norms_sqr.add(normalQuad(vC15, vC31, vC43, vC51, true));

        faces_sqr.add(new Tetrad<>(vC16, vC52, vC44, vC32));
        face_norms_sqr.add(normalQuad(vC16, vC52, vC44, vC32, true));

        faces_sqr.add(new Tetrad<>(vC17, vC32, vC45, vC53));
        face_norms_sqr.add(normalQuad(vC17, vC32, vC45, vC53, true));

        faces_sqr.add(new Tetrad<>(vC18, vC33, vC46, vC54));
        face_norms_sqr.add(normalQuad(vC18, vC33, vC46, vC54, true));

        faces_sqr.add(new Tetrad<>(vC19, vC55, vC47, vC33));
        face_norms_sqr.add(normalQuad(vC19, vC55, vC47, vC33, true));

        faces_sqr.add(new Tetrad<>(vC20, vC34, vC48, vC56));
        face_norms_sqr.add(normalQuad(vC20, vC34, vC48, vC56, true));

        faces_sqr.add(new Tetrad<>(vC21, vC57, vC49, vC34));
        face_norms_sqr.add(normalQuad(vC21, vC57, vC49, vC34, true));

        faces_sqr.add(new Tetrad<>(vC22, vC58, vC50, vC35));
        face_norms_sqr.add(normalQuad(vC22, vC58, vC50, vC35, true));

        faces_sqr.add(new Tetrad<>(vC23, vC35, vC51, vC59));
        face_norms_sqr.add(normalQuad(vC23, vC35, vC51, vC59, true));

        faces_sqr.add(new Tetrad<>(vC0, vC4, vC6, vC2));
        face_norms_sqr.add(normalQuad(vC0, vC4, vC6, vC2, true));

        faces_sqr.add(new Tetrad<>(vC1, vC3, vC7, vC5));
        face_norms_sqr.add(normalQuad(vC1, vC3, vC7, vC5, true));

        faces_sqr.add(new Tetrad<>(vC8, vC10, vC11, vC9));
        face_norms_sqr.add(normalQuad(vC8, vC10, vC11, vC9, true));

        faces_sqr.add(new Tetrad<>(vC12, vC13, vC15, vC14));
        face_norms_sqr.add(normalQuad(vC12, vC13, vC15, vC14, true));

        faces_sqr.add(new Tetrad<>(vC16, vC17, vC21, vC20));
        face_norms_sqr.add(normalQuad(vC16, vC17, vC21, vC20, true));

        faces_sqr.add(new Tetrad<>(vC18, vC22, vC23, vC19));
        face_norms_sqr.add(normalQuad(vC18, vC22, vC23, vC19, true));

        // ==== TRIANGULAR FACES ====
        faces_tri = new ArrayList<>();
        face_norms_tri = new ArrayList<>();

        faces_tri.add(new Triad<>(vC24, vC4, vC0));
        face_norms_tri.add(normalTriple(vC24, vC4, vC0, true));

        faces_tri.add(new Triad<>(vC25, vC1, vC5));
        face_norms_tri.add(normalTriple(vC25, vC1, vC5, true));

        faces_tri.add(new Triad<>(vC26, vC2, vC6));
        face_norms_tri.add(normalTriple(vC26, vC2, vC6, true));

        faces_tri.add(new Triad<>(vC27, vC7, vC3));
        face_norms_tri.add(normalTriple(vC27, vC7, vC3, true));

        faces_tri.add(new Triad<>(vC28, vC10, vC8));
        face_norms_tri.add(normalTriple(vC28, vC10, vC8, true));

        faces_tri.add(new Triad<>(vC29, vC9, vC11));
        face_norms_tri.add(normalTriple(vC29, vC9, vC11, true));

        faces_tri.add(new Triad<>(vC30, vC12, vC14));
        face_norms_tri.add(normalTriple(vC30, vC12, vC14, true));

        faces_tri.add(new Triad<>(vC31, vC15, vC13));
        face_norms_tri.add(normalTriple(vC31, vC15, vC13, true));

        faces_tri.add(new Triad<>(vC32, vC17, vC16));
        face_norms_tri.add(normalTriple(vC32, vC17, vC16, true));

        faces_tri.add(new Triad<>(vC33, vC18, vC19));
        face_norms_tri.add(normalTriple(vC33, vC18, vC19, true));

        faces_tri.add(new Triad<>(vC34, vC20, vC21));
        face_norms_tri.add(normalTriple(vC34, vC20, vC21, true));

        faces_tri.add(new Triad<>(vC35, vC23, vC22));
        face_norms_tri.add(normalTriple(vC35, vC23, vC22, true));

        faces_tri.add(new Triad<>(vC36, vC44, vC52));
        face_norms_tri.add(normalTriple(vC36, vC44, vC52, true));

        faces_tri.add(new Triad<>(vC37, vC53, vC45));
        face_norms_tri.add(normalTriple(vC37, vC53, vC45, true));

        faces_tri.add(new Triad<>(vC38, vC54, vC46));
        face_norms_tri.add(normalTriple(vC38, vC54, vC46, true));

        faces_tri.add(new Triad<>(vC39, vC47, vC55));
        face_norms_tri.add(normalTriple(vC39, vC47, vC55, true));

        faces_tri.add(new Triad<>(vC40, vC56, vC48));
        face_norms_tri.add(normalTriple(vC40, vC56, vC48, true));

        faces_tri.add(new Triad<>(vC41, vC49, vC57));
        face_norms_tri.add(normalTriple(vC41, vC49, vC57, true));

        faces_tri.add(new Triad<>(vC42, vC50, vC58));
        face_norms_tri.add(normalTriple(vC42, vC50, vC58, true));

        faces_tri.add(new Triad<>(vC43, vC59, vC51));
        face_norms_tri.add(normalTriple(vC43, vC59, vC51, true));

        // ==== PENTAGONAL FACES ====
        faces_pnt = new ArrayList<>();
        face_norms_pnt = new ArrayList<>();

        faces_pnt.add(new Pentad<>(vC24, vC52, vC16, vC20, vC56));
        face_norms_pnt.add(normalPent(vC24, vC52, vC16, vC20, vC56, true));

        faces_pnt.add(new Pentad<>(vC25, vC57, vC21, vC17, vC53));
        face_norms_pnt.add(normalPent(vC25, vC57, vC21, vC17, vC53, true));

        faces_pnt.add(new Pentad<>(vC26, vC58, vC22, vC18, vC54));
        face_norms_pnt.add(normalPent(vC26, vC58, vC22, vC18, vC54, true));

        faces_pnt.add(new Pentad<>(vC27, vC55, vC19, vC23, vC59));
        face_norms_pnt.add(normalPent(vC27, vC55, vC19, vC23, vC59, true));

        faces_pnt.add(new Pentad<>(vC28, vC36, vC0, vC2, vC38));
        face_norms_pnt.add(normalPent(vC28, vC36, vC0, vC2, vC38, true));

        faces_pnt.add(new Pentad<>(vC29, vC39, vC3, vC1, vC37));
        face_norms_pnt.add(normalPent(vC29, vC39, vC3, vC1, vC37, true));

        faces_pnt.add(new Pentad<>(vC30, vC42, vC6, vC4, vC40));
        face_norms_pnt.add(normalPent(vC30, vC42, vC6, vC4, vC40, true));

        faces_pnt.add(new Pentad<>(vC31, vC41, vC5, vC7, vC43));
        face_norms_pnt.add(normalPent(vC31, vC41, vC5, vC7, vC43, true));

        faces_pnt.add(new Pentad<>(vC32, vC44, vC8, vC9, vC45));
        face_norms_pnt.add(normalPent(vC32, vC44, vC8, vC9, vC45, true));

        faces_pnt.add(new Pentad<>(vC33, vC47, vC11, vC10, vC46));
        face_norms_pnt.add(normalPent(vC33, vC47, vC11, vC10, vC46, true));

        faces_pnt.add(new Pentad<>(vC34, vC49, vC13, vC12, vC48));
        face_norms_pnt.add(normalPent(vC34, vC49, vC13, vC12, vC48, true));

        faces_pnt.add(new Pentad<>(vC35, vC50, vC14, vC15, vC51));
        face_norms_pnt.add(normalPent(vC35, vC50, vC14, vC15, vC51, true));
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
