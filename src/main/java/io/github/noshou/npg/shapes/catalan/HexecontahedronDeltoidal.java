package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;
import static io.github.noshou.npg.nputil.VectorMath.*;
import io.github.noshou.npg.shapes.archimedean.*;

/**
 * Represents a <b>Deltoidal Hexecontahedron</b>.
 * <p> A Catalan solid composed of 60 congruent kite-shaped (deltoid) faces,
 *  120 edges, and 62 vertices.
 * <p> It is the dual of the {@link Rhombicosidodecahedron}.
 * @see <a href="https://dmccooey.com/polyhedra/DeltoidalHexecontahedron.html">
 *      Deltoidal Hexecontahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class HexecontahedronDeltoidal extends Shape {

    // 60 kite faces
    private final ArrayList<Tetrad<Tuple<Apfloat>>> faces_kte;
    private final ArrayList<Triad<Apfloat>> face_norms_kte;

    // Apfloat Constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N6 = new Apfloat("6", super.precision);
    private final Apfloat N9 = new Apfloat("9", super.precision);
    private final Apfloat N11 = new Apfloat("11", super.precision);
    private final Apfloat N15 = new Apfloat("15", super.precision);
    private final Apfloat N22 = new Apfloat("22", super.precision);
    private final Apfloat N25 = new Apfloat("25", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(N5);

    // C0 = (5 - sqrt(5)) / 4
    private final Apfloat C0 = (N5.subtract(SQRT5)).divide(N4);

    // C1 = (15 + sqrt(5)) / 22
    private final Apfloat C1 = (N15.add(SQRT5)).divide(N22);

    // C2 = sqrt(5) / 2
    private final Apfloat C2 = SQRT5.divide(N2);

    // C3 = (5 + sqrt(5)) / 6
    private final Apfloat C3 = (N5.add(SQRT5)).divide(N6);

    // C4 = (5 + 4 * sqrt(5)) / 11
    private final Apfloat C4 = (N5.add(N4.multiply(SQRT5))).divide(N11);

    // C5 = (5 + sqrt(5)) / 4
    private final Apfloat C5 = (N5.add(SQRT5)).divide(N4);

    // C6 = (5 + 3 * sqrt(5)) / 6
    private final Apfloat C6 = (N5.add(N3.multiply(SQRT5))).divide(N6);

    // C7 = (25 + 9 * sqrt(5)) / 22
    private final Apfloat C7 = (N25.add(N9.multiply(SQRT5))).divide(N22);

    // C8 = sqrt(5)
    private final Apfloat C8 = SQRT5;

    public HexecontahedronDeltoidal(
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
        Triad<Apfloat> vB0 = new Triad<>(N0, N0, C8);
        Triad<Apfloat> vB1 = new Triad<>(N0, N0, C8.multiply(NEG_N1));
        Triad<Apfloat> vB2  = new Triad<>(C8, N0, N0);
        Triad<Apfloat> vB3  = new Triad<>(C8.multiply(NEG_N1), N0, N0);
        Triad<Apfloat> vB4  = new Triad<>(N0, C8, N0);
        Triad<Apfloat> vB5  = new Triad<>(N0, C8.multiply(NEG_N1), N0);
        Triad<Apfloat> vB6  = new Triad<>(N0, C1, C7);
        Triad<Apfloat> vB7  = new Triad<>(N0, C1, C7.multiply(NEG_N1));
        Triad<Apfloat> vB8  = new Triad<>(N0, C1.multiply(NEG_N1), C7);
        Triad<Apfloat> vB9  = new Triad<>(N0, C1.multiply(NEG_N1), C7.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(C7, N0, C1);
        Triad<Apfloat> vB11 = new Triad<>(C7, N0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(C7.multiply(NEG_N1), N0, C1);
        Triad<Apfloat> vB13 = new Triad<>(C7.multiply(NEG_N1), N0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(C1, C7, N0);
        Triad<Apfloat> vB15 = new Triad<>(C1, C7.multiply(NEG_N1), N0);
        Triad<Apfloat> vB16 = new Triad<>(C1.multiply(NEG_N1), C7, N0);
        Triad<Apfloat> vB17 = new Triad<>(C1.multiply(NEG_N1), C7.multiply(NEG_N1), N0);
        Triad<Apfloat> vB18 = new Triad<>(C3, N0, C6);
        Triad<Apfloat> vB19 = new Triad<>(C3, N0, C6.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>(C3.multiply(NEG_N1), N0, C6);
        Triad<Apfloat> vB21 = new Triad<>(C3.multiply(NEG_N1), N0, C6.multiply(NEG_N1));
        Triad<Apfloat> vB22 = new Triad<>(C6, C3, N0);
        Triad<Apfloat> vB23 = new Triad<>(C6, C3.multiply(NEG_N1), N0);
        Triad<Apfloat> vB24 = new Triad<>(C6.multiply(NEG_N1), C3, N0);
        Triad<Apfloat> vB25 = new Triad<>(C6.multiply(NEG_N1), C3.multiply(NEG_N1), N0);
        Triad<Apfloat> vB26 = new Triad<>(N0, C6, C3);
        Triad<Apfloat> vB27 = new Triad<>(N0, C6, C3.multiply(NEG_N1));
        Triad<Apfloat> vB28 = new Triad<>(N0, C6.multiply(NEG_N1), C3);
        Triad<Apfloat> vB29 = new Triad<>(N0, C6.multiply(NEG_N1), C3.multiply(NEG_N1));
        Triad<Apfloat> vB30 = new Triad<>(C0, C2, C5);
        Triad<Apfloat> vB31 = new Triad<>(C0, C2, C5.multiply(NEG_N1));
        Triad<Apfloat> vB32 = new Triad<>(C0, C2.multiply(NEG_N1), C5);
        Triad<Apfloat> vB33 = new Triad<>(C0, C2.multiply(NEG_N1), C5.multiply(NEG_N1));
        Triad<Apfloat> vB34 = new Triad<>(C0.multiply(NEG_N1), C2, C5);
        Triad<Apfloat> vB35 = new Triad<>(C0.multiply(NEG_N1), C2, C5.multiply(NEG_N1));
        Triad<Apfloat> vB36 = new Triad<>(C0.multiply(NEG_N1), C2.multiply(NEG_N1), C5);
        Triad<Apfloat> vB37 = new Triad<>(C0.multiply(NEG_N1), C2.multiply(NEG_N1), C5.multiply(NEG_N1));
        Triad<Apfloat> vB38 = new Triad<>(C5, C0, C2);
        Triad<Apfloat> vB39 = new Triad<>(C5, C0, C2.multiply(NEG_N1));
        Triad<Apfloat> vB40 = new Triad<>(C5, C0.multiply(NEG_N1), C2);
        Triad<Apfloat> vB41 = new Triad<>(C5, C0.multiply(NEG_N1), C2.multiply(NEG_N1));
        Triad<Apfloat> vB42 = new Triad<>(C5.multiply(NEG_N1), C0, C2);
        Triad<Apfloat> vB43 = new Triad<>(C5.multiply(NEG_N1), C0, C2.multiply(NEG_N1));
        Triad<Apfloat> vB44 = new Triad<>(C5.multiply(NEG_N1), C0.multiply(NEG_N1), C2);
        Triad<Apfloat> vB45 = new Triad<>(C5.multiply(NEG_N1), C0.multiply(NEG_N1), C2.multiply(NEG_N1));
        Triad<Apfloat> vB46 = new Triad<>(C2, C5, C0);
        Triad<Apfloat> vB47 = new Triad<>(C2, C5, C0.multiply(NEG_N1));
        Triad<Apfloat> vB48 = new Triad<>(C2, C5.multiply(NEG_N1), C0);
        Triad<Apfloat> vB49 = new Triad<>(C2, C5.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB50 = new Triad<>(C2.multiply(NEG_N1), C5, C0);
        Triad<Apfloat> vB51 = new Triad<>(C2.multiply(NEG_N1), C5, C0.multiply(NEG_N1));
        Triad<Apfloat> vB52 = new Triad<>(C2.multiply(NEG_N1), C5.multiply(NEG_N1), C0);
        Triad<Apfloat> vB53 = new Triad<>(C2.multiply(NEG_N1), C5.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB54 = new Triad<>(C4, C4, C4);
        Triad<Apfloat> vB55 = new Triad<>(C4, C4, C4.multiply(NEG_N1));
        Triad<Apfloat> vB56 = new Triad<>(C4, C4.multiply(NEG_N1), C4);
        Triad<Apfloat> vB57 = new Triad<>(C4, C4.multiply(NEG_N1), C4.multiply(NEG_N1));
        Triad<Apfloat> vB58 = new Triad<>(C4.multiply(NEG_N1), C4, C4);
        Triad<Apfloat> vB59 = new Triad<>(C4.multiply(NEG_N1), C4, C4.multiply(NEG_N1));
        Triad<Apfloat> vB60 = new Triad<>(C4.multiply(NEG_N1), C4.multiply(NEG_N1), C4);
        Triad<Apfloat> vB61 = new Triad<>(C4.multiply(NEG_N1), C4.multiply(NEG_N1), C4.multiply(NEG_N1));

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
        Triad<Apfloat> vC60 = mult(normalize(vB60), super.getRadius().toString());
        Triad<Apfloat> vC61 = mult(normalize(vB61), super.getRadius().toString());


        // ==== KITE FACES ====
        faces_kte = new ArrayList<>();
        face_norms_kte = new ArrayList<>();

        faces_kte.add(new Tetrad<>(vC18, vC32, vC56, vC40));
        face_norms_kte.add(normalQuad(vC18, vC32, vC56, vC40, true));

        faces_kte.add(new Tetrad<>(vC18, vC0, vC8, vC32));
        face_norms_kte.add(normalQuad(vC18, vC0, vC8, vC32, true));

        faces_kte.add(new Tetrad<>(vC18, vC40, vC10, vC38));
        face_norms_kte.add(normalQuad(vC18, vC40, vC10, vC38, true));

        faces_kte.add(new Tetrad<>(vC18, vC38, vC54, vC30));
        face_norms_kte.add(normalQuad(vC18, vC38, vC54, vC30, true));

        faces_kte.add(new Tetrad<>(vC18, vC30, vC6, vC0));
        face_norms_kte.add(normalQuad(vC18, vC30, vC6, vC0, true));

        faces_kte.add(new Tetrad<>(vC19, vC1, vC7, vC31));
        face_norms_kte.add(normalQuad(vC19, vC1, vC7, vC31, true));

        faces_kte.add(new Tetrad<>(vC19, vC31, vC55, vC39));
        face_norms_kte.add(normalQuad(vC19, vC31, vC55, vC39, true));

        faces_kte.add(new Tetrad<>(vC19, vC39, vC11, vC41));
        face_norms_kte.add(normalQuad(vC19, vC39, vC11, vC41, true));

        faces_kte.add(new Tetrad<>(vC19, vC41, vC57, vC33));
        face_norms_kte.add(normalQuad(vC19, vC41, vC57, vC33, true));

        faces_kte.add(new Tetrad<>(vC19, vC33, vC9, vC1));
        face_norms_kte.add(normalQuad(vC19, vC33, vC9, vC1, true));

        faces_kte.add(new Tetrad<>(vC20, vC0, vC6, vC34));
        face_norms_kte.add(normalQuad(vC20, vC0, vC6, vC34, true));

        faces_kte.add(new Tetrad<>(vC20, vC34, vC58, vC42));
        face_norms_kte.add(normalQuad(vC20, vC34, vC58, vC42, true));

        faces_kte.add(new Tetrad<>(vC20, vC42, vC12, vC44));
        face_norms_kte.add(normalQuad(vC20, vC42, vC12, vC44, true));

        faces_kte.add(new Tetrad<>(vC20, vC44, vC60, vC36));
        face_norms_kte.add(normalQuad(vC20, vC44, vC60, vC36, true));

        faces_kte.add(new Tetrad<>(vC20, vC36, vC8, vC0));
        face_norms_kte.add(normalQuad(vC20, vC36, vC8, vC0, true));

        faces_kte.add(new Tetrad<>(vC21, vC1, vC9, vC37));
        face_norms_kte.add(normalQuad(vC21, vC1, vC9, vC37, true));

        faces_kte.add(new Tetrad<>(vC21, vC37, vC61, vC45));
        face_norms_kte.add(normalQuad(vC21, vC37, vC61, vC45, true));

        faces_kte.add(new Tetrad<>(vC21, vC45, vC13, vC43));
        face_norms_kte.add(normalQuad(vC21, vC45, vC13, vC43, true));

        faces_kte.add(new Tetrad<>(vC21, vC43, vC59, vC35));
        face_norms_kte.add(normalQuad(vC21, vC43, vC59, vC35, true));

        faces_kte.add(new Tetrad<>(vC21, vC35, vC7, vC1));
        face_norms_kte.add(normalQuad(vC21, vC35, vC7, vC1, true));

        faces_kte.add(new Tetrad<>(vC22, vC2, vC11, vC39));
        face_norms_kte.add(normalQuad(vC22, vC2, vC11, vC39, true));

        faces_kte.add(new Tetrad<>(vC22, vC39, vC55, vC47));
        face_norms_kte.add(normalQuad(vC22, vC39, vC55, vC47, true));

        faces_kte.add(new Tetrad<>(vC22, vC47, vC14, vC46));
        face_norms_kte.add(normalQuad(vC22, vC47, vC14, vC46, true));

        faces_kte.add(new Tetrad<>(vC22, vC46, vC54, vC38));
        face_norms_kte.add(normalQuad(vC22, vC46, vC54, vC38, true));

        faces_kte.add(new Tetrad<>(vC22, vC38, vC10, vC2));
        face_norms_kte.add(normalQuad(vC22, vC38, vC10, vC2, true));

        faces_kte.add(new Tetrad<>(vC23, vC2, vC10, vC40));
        face_norms_kte.add(normalQuad(vC23, vC2, vC10, vC40, true));

        faces_kte.add(new Tetrad<>(vC23, vC40, vC56, vC48));
        face_norms_kte.add(normalQuad(vC23, vC40, vC56, vC48, true));

        faces_kte.add(new Tetrad<>(vC23, vC48, vC15, vC49));
        face_norms_kte.add(normalQuad(vC23, vC48, vC15, vC49, true));

        faces_kte.add(new Tetrad<>(vC23, vC49, vC57, vC41));
        face_norms_kte.add(normalQuad(vC23, vC49, vC57, vC41, true));

        faces_kte.add(new Tetrad<>(vC23, vC41, vC11, vC2));
        face_norms_kte.add(normalQuad(vC23, vC41, vC11, vC2, true));

        faces_kte.add(new Tetrad<>(vC24, vC3, vC12, vC42));
        face_norms_kte.add(normalQuad(vC24, vC3, vC12, vC42, true));

        faces_kte.add(new Tetrad<>(vC24, vC42, vC58, vC50));
        face_norms_kte.add(normalQuad(vC24, vC42, vC58, vC50, true));

        faces_kte.add(new Tetrad<>(vC24, vC50, vC16, vC51));
        face_norms_kte.add(normalQuad(vC24, vC50, vC16, vC51, true));

        faces_kte.add(new Tetrad<>(vC24, vC51, vC59, vC43));
        face_norms_kte.add(normalQuad(vC24, vC51, vC59, vC43, true));

        faces_kte.add(new Tetrad<>(vC24, vC43, vC13, vC3));
        face_norms_kte.add(normalQuad(vC24, vC43, vC13, vC3, true));

        faces_kte.add(new Tetrad<>(vC25, vC3, vC13, vC45));
        face_norms_kte.add(normalQuad(vC25, vC3, vC13, vC45, true));

        faces_kte.add(new Tetrad<>(vC25, vC45, vC61, vC53));
        face_norms_kte.add(normalQuad(vC25, vC45, vC61, vC53, true));

        faces_kte.add(new Tetrad<>(vC25, vC53, vC17, vC52));
        face_norms_kte.add(normalQuad(vC25, vC53, vC17, vC52, true));

        faces_kte.add(new Tetrad<>(vC25, vC52, vC60, vC44));
        face_norms_kte.add(normalQuad(vC25, vC52, vC60, vC44, true));

        faces_kte.add(new Tetrad<>(vC25, vC44, vC12, vC3));
        face_norms_kte.add(normalQuad(vC25, vC44, vC12, vC3, true));

        faces_kte.add(new Tetrad<>(vC26, vC4, vC16, vC50));
        face_norms_kte.add(normalQuad(vC26, vC4, vC16, vC50, true));

        faces_kte.add(new Tetrad<>(vC26, vC50, vC58, vC34));
        face_norms_kte.add(normalQuad(vC26, vC50, vC58, vC34, true));

        faces_kte.add(new Tetrad<>(vC26, vC34, vC6, vC30));
        face_norms_kte.add(normalQuad(vC26, vC34, vC6, vC30, true));

        faces_kte.add(new Tetrad<>(vC26, vC30, vC54, vC46));
        face_norms_kte.add(normalQuad(vC26, vC30, vC54, vC46, true));

        faces_kte.add(new Tetrad<>(vC26, vC46, vC14, vC4));
        face_norms_kte.add(normalQuad(vC26, vC46, vC14, vC4, true));

        faces_kte.add(new Tetrad<>(vC27, vC4, vC14, vC47));
        face_norms_kte.add(normalQuad(vC27, vC4, vC14, vC47, true));

        faces_kte.add(new Tetrad<>(vC27, vC47, vC55, vC31));
        face_norms_kte.add(normalQuad(vC27, vC47, vC55, vC31, true));

        faces_kte.add(new Tetrad<>(vC27, vC31, vC7, vC35));
        face_norms_kte.add(normalQuad(vC27, vC31, vC7, vC35, true));

        faces_kte.add(new Tetrad<>(vC27, vC35, vC59, vC51));
        face_norms_kte.add(normalQuad(vC27, vC35, vC59, vC51, true));

        faces_kte.add(new Tetrad<>(vC27, vC51, vC16, vC4));
        face_norms_kte.add(normalQuad(vC27, vC51, vC16, vC4, true));

        faces_kte.add(new Tetrad<>(vC28, vC5, vC15, vC48));
        face_norms_kte.add(normalQuad(vC28, vC5, vC15, vC48, true));

        faces_kte.add(new Tetrad<>(vC28, vC48, vC56, vC32));
        face_norms_kte.add(normalQuad(vC28, vC48, vC56, vC32, true));

        faces_kte.add(new Tetrad<>(vC28, vC32, vC8, vC36));
        face_norms_kte.add(normalQuad(vC28, vC32, vC8, vC36, true));

        faces_kte.add(new Tetrad<>(vC28, vC36, vC60, vC52));
        face_norms_kte.add(normalQuad(vC28, vC36, vC60, vC52, true));

        faces_kte.add(new Tetrad<>(vC28, vC52, vC17, vC5));
        face_norms_kte.add(normalQuad(vC28, vC52, vC17, vC5, true));

        faces_kte.add(new Tetrad<>(vC29, vC5, vC17, vC53));
        face_norms_kte.add(normalQuad(vC29, vC5, vC17, vC53, true));

        faces_kte.add(new Tetrad<>(vC29, vC53, vC61, vC37));
        face_norms_kte.add(normalQuad(vC29, vC53, vC61, vC37, true));

        faces_kte.add(new Tetrad<>(vC29, vC37, vC9, vC33));
        face_norms_kte.add(normalQuad(vC29, vC37, vC9, vC33, true));

        faces_kte.add(new Tetrad<>(vC29, vC33, vC57, vC49));
        face_norms_kte.add(normalQuad(vC29, vC33, vC57, vC49, true));

        faces_kte.add(new Tetrad<>(vC29, vC49, vC15, vC5));
        face_norms_kte.add(normalQuad(vC29, vC49, vC15, vC5, true));
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < faces_kte.size(); i++) {

            // Kite  faces
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_kte.get(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_kte.get(i);

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
