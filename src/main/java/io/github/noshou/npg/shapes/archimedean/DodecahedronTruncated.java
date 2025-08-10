package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b>Truncated Dodecahedron</b>.
 * <p>
 * The truncated dodecahedron is one of the 13 Archimedean solids, formed by truncating (cutting off)
 * the vertices of a regular dodecahedron. This truncation replaces the original pentagonal faces with
 * decagonal faces and the original vertices with triangular faces, resulting in a polyhedron with 32 faces:
 * 12 regular decagons and 20 equilateral triangles. It has 90 edges and 60 vertices, exhibiting icosahedral symmetry.
 * @see <a href="https://dmccooey.com/polyhedra/TruncatedDodecahedron.html">
 *      Truncated Dodecahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class DodecahedronTruncated extends Shape {

    // 20 triangular faces
    private final Icosad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Icosad<Tuple<Apfloat>> face_norms_tri;

    // 12 decagonal faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_dec;
    private final Dodecad<Tuple<Apfloat>> face_norms_dec;

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(N5);
    private final Apfloat FRAC_1_over_2 = N1.divide(N2);

    // C0 = (3 + sqrt(5)) / 4
    private final Apfloat C0 = (N3.add(SQRT5)).divide(N4);

    // C1 = (1 + sqrt(5)) / 2
    private final Apfloat C1 = (N1.add(SQRT5)).divide(N2);

    // C2 = (2 + sqrt(5)) / 2
    private final Apfloat C2 = (N2.add(SQRT5)).divide(N2);

    // C3 = (3 + sqrt(5)) / 2
    private final Apfloat C3 = (N3.add(SQRT5)).divide(N2);

    // C4 = (5 + 3 * sqrt(5)) / 4
    private final Apfloat C4 = (N5.add(N3.multiply(SQRT5))).divide(N4);

    public DodecahedronTruncated(
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

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vB0 = new Triad<>(N0,FRAC_1_over_2, C4);
        Triad<Apfloat> vB1 = new Triad<>(N0, FRAC_1_over_2, C4.multiply(NEG_N1));
        Triad<Apfloat> vB2 = new Triad<>(N0, FRAC_1_over_2.multiply(NEG_N1), C4);
        Triad<Apfloat> vB3 = new Triad<>(N0, FRAC_1_over_2.multiply(NEG_N1), C4.multiply(NEG_N1));
        Triad<Apfloat> vB4 = new Triad<>(C4, N0, FRAC_1_over_2);
        Triad<Apfloat> vB5 = new Triad<>(C4, N0, FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB6 = new Triad<>(C4.multiply(NEG_N1), N0, FRAC_1_over_2);
        Triad<Apfloat> vB7 = new Triad<>(C4.multiply(NEG_N1), N0, FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB8 = new Triad<>(FRAC_1_over_2, C4, N0);
        Triad<Apfloat> vB9 = new Triad<>(FRAC_1_over_2, C4.multiply(NEG_N1), N0);
        Triad<Apfloat> vB10 = new Triad<>(FRAC_1_over_2.multiply(NEG_N1), C4, N0);
        Triad<Apfloat> vB11 = new Triad<>(FRAC_1_over_2.multiply(NEG_N1), C4.multiply(NEG_N1), N0);
        Triad<Apfloat> vB12 = new Triad<>(FRAC_1_over_2, C0, C3);
        Triad<Apfloat> vB13 = new Triad<>(FRAC_1_over_2, C0, C3.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(FRAC_1_over_2, C0.multiply(NEG_N1), C3);
        Triad<Apfloat> vB15 = new Triad<>(FRAC_1_over_2, C0.multiply(NEG_N1), C3.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>(FRAC_1_over_2.multiply(NEG_N1), C0, C3);
        Triad<Apfloat> vB17 = new Triad<>(FRAC_1_over_2.multiply(NEG_N1), C0, C3.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>(FRAC_1_over_2.multiply(NEG_N1), C0.multiply(NEG_N1), C3);
        Triad<Apfloat> vB19 = new Triad<>(FRAC_1_over_2.multiply(NEG_N1), C0.multiply(NEG_N1), C3.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>(C3, FRAC_1_over_2, C0);
        Triad<Apfloat> vB21 = new Triad<>(C3, FRAC_1_over_2, C0.multiply(NEG_N1));
        Triad<Apfloat> vB22 = new Triad<>(C3, FRAC_1_over_2.multiply(NEG_N1), C0);
        Triad<Apfloat> vB23 = new Triad<>(C3, FRAC_1_over_2.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB24 = new Triad<>(C3.multiply(NEG_N1), FRAC_1_over_2, C0);
        Triad<Apfloat> vB25 = new Triad<>(C3.multiply(NEG_N1), FRAC_1_over_2, C0.multiply(NEG_N1));
        Triad<Apfloat> vB26 = new Triad<>(C3.multiply(NEG_N1), FRAC_1_over_2.multiply(NEG_N1), C0);
        Triad<Apfloat> vB27 = new Triad<>(C3.multiply(NEG_N1), FRAC_1_over_2.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB28 = new Triad<>(C0, C3, FRAC_1_over_2);
        Triad<Apfloat> vB29 = new Triad<>(C0, C3, FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB30 = new Triad<>(C0, C3.multiply(NEG_N1), FRAC_1_over_2);
        Triad<Apfloat> vB31 = new Triad<>(C0, C3.multiply(NEG_N1), FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB32 = new Triad<>(C0.multiply(NEG_N1), C3, FRAC_1_over_2);
        Triad<Apfloat> vB33 = new Triad<>(C0.multiply(NEG_N1), C3, FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB34 = new Triad<>(C0.multiply(NEG_N1), C3.multiply(NEG_N1), FRAC_1_over_2);
        Triad<Apfloat> vB35 = new Triad<>(C0.multiply(NEG_N1), C3.multiply(NEG_N1), FRAC_1_over_2.multiply(NEG_N1));
        Triad<Apfloat> vB36 = new Triad<>(C0, C1, C2);
        Triad<Apfloat> vB37 = new Triad<>(C0, C1, C2.multiply(NEG_N1));
        Triad<Apfloat> vB38 = new Triad<>(C0, C1.multiply(NEG_N1), C2);
        Triad<Apfloat> vB39 = new Triad<>(C0, C1.multiply(NEG_N1), C2.multiply(NEG_N1));
        Triad<Apfloat> vB40 = new Triad<>(C0.multiply(NEG_N1), C1, C2);
        Triad<Apfloat> vB41 = new Triad<>(C0.multiply(NEG_N1), C1, C2.multiply(NEG_N1));
        Triad<Apfloat> vB42 = new Triad<>(C0.multiply(NEG_N1), C1.multiply(NEG_N1), C2);
        Triad<Apfloat> vB43 = new Triad<>(C0.multiply(NEG_N1), C1.multiply(NEG_N1), C2.multiply(NEG_N1));
        Triad<Apfloat> vB44 = new Triad<>(C2, C0, C1);
        Triad<Apfloat> vB45 = new Triad<>(C2, C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB46 = new Triad<>(C2, C0.multiply(NEG_N1), C1);
        Triad<Apfloat> vB47 = new Triad<>(C2, C0.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB48 = new Triad<>(C2.multiply(NEG_N1), C0, C1);
        Triad<Apfloat> vB49 = new Triad<>(C2.multiply(NEG_N1), C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB50 = new Triad<>(C2.multiply(NEG_N1), C0.multiply(NEG_N1), C1);
        Triad<Apfloat> vB51 = new Triad<>(C2.multiply(NEG_N1), C0.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB52 = new Triad<>(C1, C2, C0);
        Triad<Apfloat> vB53 = new Triad<>(C1, C2, C0.multiply(NEG_N1));
        Triad<Apfloat> vB54 = new Triad<>(C1, C2.multiply(NEG_N1), C0);
        Triad<Apfloat> vB55 = new Triad<>(C1, C2.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB56 = new Triad<>(C1.multiply(NEG_N1), C2, C0);
        Triad<Apfloat> vB57 = new Triad<>(C1.multiply(NEG_N1), C2, C0.multiply(NEG_N1));
        Triad<Apfloat> vB58 = new Triad<>(C1.multiply(NEG_N1), C2.multiply(NEG_N1), C0);
        Triad<Apfloat> vB59 = new Triad<>(C1.multiply(NEG_N1), C2.multiply(NEG_N1), C0.multiply(NEG_N1));

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = mult(normalize(vB0),  super.getRadius().toString());
        Triad<Apfloat> vC1  = mult(normalize(vB1),  super.getRadius().toString());
        Triad<Apfloat> vC2  = mult(normalize(vB2),  super.getRadius().toString());
        Triad<Apfloat> vC3  = mult(normalize(vB3),  super.getRadius().toString());
        Triad<Apfloat> vC4  = mult(normalize(vB4),  super.getRadius().toString());
        Triad<Apfloat> vC5  = mult(normalize(vB5),  super.getRadius().toString());
        Triad<Apfloat> vC6  = mult(normalize(vB6),  super.getRadius().toString());
        Triad<Apfloat> vC7  = mult(normalize(vB7),  super.getRadius().toString());
        Triad<Apfloat> vC8  = mult(normalize(vB8),  super.getRadius().toString());
        Triad<Apfloat> vC9  = mult(normalize(vB9),  super.getRadius().toString());
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

        // ==== DODECAGONAL FACES ====
        Decad<Tuple<Apfloat>> dec0 = new Decad<>(vC0, vC2, vC14, vC38, vC46, vC22, vC20, vC44, vC36, vC12);
        Triad<Apfloat> dec0_norm = normalDeca(vC0, vC2, vC14, vC38, vC46, vC22, vC20, vC44, vC36, vC12, true);
        Decad<Tuple<Apfloat>> dec1 = new Decad<>(vC1, vC3, vC19, vC43, vC51, vC27, vC25, vC49, vC41, vC17);
        Triad<Apfloat> dec1_norm = normalDeca(vC1, vC3, vC19, vC43, vC51, vC27, vC25, vC49, vC41, vC17, true);

        Decad<Tuple<Apfloat>> dec2 = new Decad<>(vC2, vC0, vC16, vC40, vC48, vC24, vC26, vC50, vC42, vC18);
        Triad<Apfloat> dec2_norm = normalDeca(vC2, vC0, vC16, vC40, vC48, vC24, vC26, vC50, vC42, vC18, true);
        Decad<Tuple<Apfloat>> dec3 = new Decad<>(vC3, vC1, vC13, vC37, vC45, vC21, vC23, vC47, vC39, vC15);
        Triad<Apfloat> dec3_norm = normalDeca(vC3, vC1, vC13, vC37, vC45, vC21, vC23, vC47, vC39, vC15, true);
        Decad<Tuple<Apfloat>> dec4 = new Decad<>(vC4, vC5, vC21, vC45, vC53, vC29, vC28, vC52, vC44, vC20);
        Triad<Apfloat> dec4_norm = normalDeca(vC4, vC5, vC21, vC45, vC53, vC29, vC28, vC52, vC44, vC20, true);
        Decad<Tuple<Apfloat>> dec5 = new Decad<>(vC5, vC4, vC22, vC46, vC54, vC30, vC31, vC55, vC47, vC23);
        Triad<Apfloat> dec5_norm = normalDeca(vC5, vC4, vC22, vC46, vC54, vC30, vC31, vC55, vC47, vC23, true);
        Decad<Tuple<Apfloat>> dec6 = new Decad<>(vC6, vC7, vC27, vC51, vC59, vC35, vC34, vC58, vC50, vC26);
        Triad<Apfloat> dec6_norm = normalDeca(vC6, vC7, vC27, vC51, vC59, vC35, vC34, vC58, vC50, vC26, true);
        Decad<Tuple<Apfloat>> dec7 = new Decad<>(vC7, vC6, vC24, vC48, vC56, vC32, vC33, vC57, vC49, vC25);
        Triad<Apfloat> dec7_norm = normalDeca(vC7, vC6, vC24, vC48, vC56, vC32, vC33, vC57, vC49, vC25, true);
        Decad<Tuple<Apfloat>> dec8 = new Decad<>(vC8, vC10, vC32, vC56, vC40, vC16, vC12, vC36, vC52, vC28);
        Triad<Apfloat> dec8_norm = normalDeca(vC8, vC10, vC32, vC56, vC40, vC16, vC12, vC36, vC52, vC28, true);
        Decad<Tuple<Apfloat>> dec9 = new Decad<>(vC9, vC11, vC35, vC59, vC43, vC19, vC15, vC39, vC55, vC31);
        Triad<Apfloat> dec9_norm = normalDeca(vC9, vC11, vC35, vC59, vC43, vC19, vC15, vC39, vC55, vC31, true);
        Decad<Tuple<Apfloat>> dec10 = new Decad<>(vC10, vC8, vC29, vC53, vC37, vC13, vC17, vC41, vC57, vC33);
        Triad<Apfloat> dec10_norm = normalDeca(vC10, vC8, vC29, vC53, vC37, vC13, vC17, vC41, vC57, vC33, true);
        Decad<Tuple<Apfloat>> dec11 = new Decad<>(vC11, vC9, vC30, vC54, vC38, vC14, vC18, vC42, vC58, vC34);
        Triad<Apfloat> dec11_norm = normalDeca(vC11, vC9, vC30, vC54, vC38, vC14, vC18, vC42, vC58, vC34, true);
        faces_dec = new Dodecad<>(
                dec0, dec1, dec2,
                dec3, dec4, dec5,
                dec6, dec7, dec8,
                dec9, dec10, dec11
        );
        face_norms_dec = new Dodecad<>(
                dec0_norm, dec1_norm, dec2_norm,
                dec3_norm, dec4_norm, dec5_norm,
                dec6_norm, dec7_norm, dec8_norm,
                dec9_norm, dec10_norm, dec11_norm
        );

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC0, vC12, vC16);
        Triad<Apfloat> tri0_norm = normalTriple(vC0, vC12, vC16, true);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC1, vC17, vC13);
        Triad<Apfloat> tri1_norm = normalTriple(vC1, vC17, vC13, true);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC2, vC18, vC14);
        Triad<Apfloat> tri2_norm = normalTriple(vC2, vC18, vC14, true);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC3, vC15, vC19);
        Triad<Apfloat> tri3_norm = normalTriple(vC3, vC15, vC19, true);
        Triad<Tuple<Apfloat>> tri4 = new Triad<>(vC4, vC20, vC22);
        Triad<Apfloat> tri4_norm = normalTriple(vC4, vC20, vC22, true);
        Triad<Tuple<Apfloat>> tri5 = new Triad<>(vC5, vC23, vC21);
        Triad<Apfloat> tri5_norm = normalTriple(vC5, vC23, vC21, true);
        Triad<Tuple<Apfloat>> tri6 = new Triad<>(vC6, vC26, vC24);
        Triad<Apfloat> tri6_norm = normalTriple(vC6, vC26, vC24, true);
        Triad<Tuple<Apfloat>> tri7 = new Triad<>(vC7, vC25, vC27);
        Triad<Apfloat> tri7_norm = normalTriple(vC7, vC25, vC27, true);
        Triad<Tuple<Apfloat>> tri8 = new Triad<>(vC8, vC28, vC29);
        Triad<Apfloat> tri8_norm = normalTriple(vC8, vC28, vC29, true);
        Triad<Tuple<Apfloat>> tri9 = new Triad<>(vC9, vC31, vC30);
        Triad<Apfloat> tri9_norm = normalTriple(vC9, vC31, vC30, true);
        Triad<Tuple<Apfloat>> tri10 = new Triad<>(vC10, vC33, vC32);
        Triad<Apfloat> tri10_norm = normalTriple(vC10, vC33, vC32, true);
        Triad<Tuple<Apfloat>> tri11 = new Triad<>(vC11, vC34, vC35);
        Triad<Apfloat> tri11_norm = normalTriple(vC11, vC34, vC35, true);
        Triad<Tuple<Apfloat>> tri12 = new Triad<>(vC36, vC44, vC52);
        Triad<Apfloat> tri12_norm = normalTriple(vC36, vC44, vC52, true);
        Triad<Tuple<Apfloat>> tri13 = new Triad<>(vC37, vC53, vC45);
        Triad<Apfloat> tri13_norm = normalTriple(vC37, vC53, vC45, true);
        Triad<Tuple<Apfloat>> tri14 = new Triad<>(vC38, vC54, vC46);
        Triad<Apfloat> tri14_norm = normalTriple(vC38, vC54, vC46, true);
        Triad<Tuple<Apfloat>> tri15 = new Triad<>(vC39, vC47, vC55);
        Triad<Apfloat> tri15_norm = normalTriple(vC39, vC47, vC55, true);
        Triad<Tuple<Apfloat>> tri16 = new Triad<>(vC40, vC56, vC48);
        Triad<Apfloat> tri16_norm = normalTriple(vC40, vC56, vC48, true);
        Triad<Tuple<Apfloat>> tri17 = new Triad<>(vC41, vC49, vC57);
        Triad<Apfloat> tri17_norm = normalTriple(vC41, vC49, vC57, true);
        Triad<Tuple<Apfloat>> tri18 = new Triad<>(vC42, vC50, vC58);
        Triad<Apfloat> tri18_norm = normalTriple(vC42, vC50, vC58, true);
        Triad<Tuple<Apfloat>> tri19 = new Triad<>(vC43, vC59, vC51);
        Triad<Apfloat> tri19_norm = normalTriple(vC43, vC59, vC51, true);
        faces_tri = new Icosad<>(
                tri0, tri1, tri2, tri3,
                tri4, tri5, tri6, tri7,
                tri8, tri9, tri10, tri11,
                tri12, tri13, tri14, tri15,
                tri16, tri17, tri18, tri19
        );
        face_norms_tri = new Icosad<>(
                tri0_norm, tri1_norm, tri2_norm, tri3_norm,
                tri4_norm, tri5_norm, tri6_norm, tri7_norm,
                tri8_norm, tri9_norm, tri10_norm, tri11_norm,
                tri12_norm, tri13_norm, tri14_norm, tri15_norm,
                tri16_norm, tri17_norm, tri18_norm, tri19_norm
        );
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < 20; i++) {
            // === Check decagonal faces ===
            if (i < face_norms_dec.fetchSize()) {
                Decad<Tuple<Apfloat>> face = (Decad<Tuple<Apfloat>>) faces_dec.fetch(i);
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
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_dec.fetch(i);
                Apfloat d = dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }


            // === Check triangular faces ===
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
        return true;
    }
}
