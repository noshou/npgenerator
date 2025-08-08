package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/BiscribedTruncatedIcosahedron.txt>...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class IcosahedronTruncatedBiscribed extends Shape {

    // 20 hexagonal faces
    private final Icosad<Tuple<Tuple<Apfloat>>> faces_hex;
    private final Icosad<Tuple<Apfloat>> face_norms_hex;

    // 12 pentagonal faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_pnt;
    private final Dodecad<Tuple<Apfloat>> face_norms_pnt;

    // Apfloat constants
    private final Apfloat NEG_N1    = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3  = new Apfloat("3", super.precision);
    private final Apfloat SQRT3 = ApfloatMath.sqrt(N3);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(N5);
    private final Apfloat SQRT15 = ApfloatMath.sqrt(new Apfloat("15", super.precision));
    private final Apfloat SQRT_5_plus_2mulSQRT5 = ApfloatMath.sqrt(N5.add(SQRT5.multiply(N2)));
    private final Apfloat SQRT_5_mul_2mulSQRT5 = ApfloatMath.sqrt(N2.multiply(N5.add(SQRT5)));
    private final Apfloat SQRT_2_mul_5minSQRT5 = ApfloatMath.sqrt(N2.multiply(N5.subtract(SQRT5)));

    // C0 = (sqrt(2 * (5 + sqrt(5))) / 2) - sqrt(3)
    private final Apfloat C0 = SQRT_5_mul_2mulSQRT5.divide(N2).subtract(SQRT3);

    // C1 = (5 * sqrt(3) - sqrt(15) - sqrt(2 * (5 + sqrt(5)))) / 4
    private final Apfloat C1 = SQRT3.multiply(N5).subtract(SQRT15).subtract(SQRT_5_mul_2mulSQRT5).divide(N4);

    // C2 = (sqrt(15) - sqrt(5 + 2 * sqrt(5))) / 2
    private final Apfloat C2 = SQRT15.subtract(SQRT_5_plus_2mulSQRT5).divide(N2);

    // C3 = (sqrt(3) - sqrt(15) + sqrt(2 * (5 + sqrt(5)))) / 4
    private final Apfloat C3 = SQRT3.subtract(SQRT15).add(SQRT_5_mul_2mulSQRT5).divide(N4);

    // C4 = (sqrt(5 + 2 * sqrt(5)) - sqrt(3)) / 2
    private final Apfloat C4 = SQRT_5_plus_2mulSQRT5.subtract(SQRT3).divide(N2);

    // C5 = (sqrt(3) + sqrt(15) - sqrt(2 * (5 - sqrt(5)))) / 4
    private final Apfloat C5 = SQRT3.add(SQRT15).subtract(SQRT_2_mul_5minSQRT5).divide(N4);

    // C6 = (3 * sqrt(3) - sqrt(15) + sqrt(2 * (5 - sqrt(5)))) / 4
    private final Apfloat C6 = SQRT3.multiply(N3).subtract(SQRT15).add(SQRT_2_mul_5minSQRT5).divide(N4);

    // C7 = (sqrt(15) - sqrt(3)) / 2
    private final Apfloat C7 = SQRT15.subtract(SQRT3).divide(N2);


    public IcosahedronTruncatedBiscribed(
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
        Triad<Apfloat> vB0 = new Triad<>(C0, N0, C7);
        Triad<Apfloat> vB1 = new Triad<>(C0, N0, C7.multiply(NEG_N1));
        Triad<Apfloat> vB2 = new Triad<>(C0.multiply(NEG_N1), N0, C7);
        Triad<Apfloat> vB3 = new Triad<>(C0.multiply(NEG_N1), N0, C7.multiply(NEG_N1));
        Triad<Apfloat> vB4 = new Triad<>(C7, C0, N0);
        Triad<Apfloat> vB5 = new Triad<>(C7, C0.multiply(NEG_N1), N0);
        Triad<Apfloat> vB6 = new Triad<>(C7.multiply(NEG_N1), C0, N0);
        Triad<Apfloat> vB7 = new Triad<>(C7.multiply(NEG_N1), C0.multiply(NEG_N1), N0);
        Triad<Apfloat> vB8 = new Triad<>(N0, C7, C0);
        Triad<Apfloat> vB9 = new Triad<>(N0, C7, C0.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(N0, C7.multiply(NEG_N1), C0);
        Triad<Apfloat> vB11 = new Triad<>(N0, C7.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(C3, C2, C6);
        Triad<Apfloat> vB13 = new Triad<>(C3, C2, C6.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(C3, C2.multiply(NEG_N1), C6);
        Triad<Apfloat> vB15 = new Triad<>(C3, C2.multiply(NEG_N1), C6.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>(C3.multiply(NEG_N1), C2, C6);
        Triad<Apfloat> vB17 = new Triad<>(C3.multiply(NEG_N1), C2, C6.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>(C3.multiply(NEG_N1), C2.multiply(NEG_N1), C6);
        Triad<Apfloat> vB19 = new Triad<>(C3.multiply(NEG_N1), C2.multiply(NEG_N1), C6.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>(C6, C3, C2);
        Triad<Apfloat> vB21 = new Triad<>(C6, C3, C2.multiply(NEG_N1));
        Triad<Apfloat> vB22 = new Triad<>(C6, C3.multiply(NEG_N1), C2);
        Triad<Apfloat> vB23 = new Triad<>(C6, C3.multiply(NEG_N1), C2.multiply(NEG_N1));
        Triad<Apfloat> vB24 = new Triad<>(C6.multiply(NEG_N1), C3, C2);
        Triad<Apfloat> vB25 = new Triad<>(C6.multiply(NEG_N1), C3, C2.multiply(NEG_N1));
        Triad<Apfloat> vB26 = new Triad<>(C6.multiply(NEG_N1), C3.multiply(NEG_N1), C2);
        Triad<Apfloat> vB27 = new Triad<>(C6.multiply(NEG_N1), C3.multiply(NEG_N1), C2.multiply(NEG_N1));
        Triad<Apfloat> vB28 = new Triad<>(C2, C6, C3);
        Triad<Apfloat> vB29 = new Triad<>(C2, C6, C3.multiply(NEG_N1));
        Triad<Apfloat> vB30 = new Triad<>(C2, C6.multiply(NEG_N1), C3);
        Triad<Apfloat> vB31 = new Triad<>(C2, C6.multiply(NEG_N1), C3.multiply(NEG_N1));
        Triad<Apfloat> vB32 = new Triad<>(C2.multiply(NEG_N1), C6, C3);
        Triad<Apfloat> vB33 = new Triad<>(C2.multiply(NEG_N1), C6, C3.multiply(NEG_N1));
        Triad<Apfloat> vB34 = new Triad<>(C2.multiply(NEG_N1), C6.multiply(NEG_N1), C3);
        Triad<Apfloat> vB35 = new Triad<>(C2.multiply(NEG_N1), C6.multiply(NEG_N1), C3.multiply(NEG_N1));
        Triad<Apfloat> vB36 = new Triad<>(C1, C4, C5);
        Triad<Apfloat> vB37 = new Triad<>(C1, C4, C5.multiply(NEG_N1));
        Triad<Apfloat> vB38 = new Triad<>(C1, C4.multiply(NEG_N1), C5);
        Triad<Apfloat> vB39 = new Triad<>(C1, C4.multiply(NEG_N1), C5.multiply(NEG_N1));
        Triad<Apfloat> vB40 = new Triad<>(C1.multiply(NEG_N1), C4, C5);
        Triad<Apfloat> vB41 = new Triad<>(C1.multiply(NEG_N1), C4, C5.multiply(NEG_N1));
        Triad<Apfloat> vB42 = new Triad<>(C1.multiply(NEG_N1), C4.multiply(NEG_N1), C5);
        Triad<Apfloat> vB43 = new Triad<>(C1.multiply(NEG_N1), C4.multiply(NEG_N1), C5.multiply(NEG_N1));
        Triad<Apfloat> vB44 = new Triad<>(C5, C1, C4);
        Triad<Apfloat> vB45 = new Triad<>(C5, C1, C4.multiply(NEG_N1));
        Triad<Apfloat> vB46 = new Triad<>(C5, C1.multiply(NEG_N1), C4);
        Triad<Apfloat> vB47 = new Triad<>(C5, C1.multiply(NEG_N1), C4.multiply(NEG_N1));
        Triad<Apfloat> vB48 = new Triad<>(C5.multiply(NEG_N1), C1, C4);
        Triad<Apfloat> vB49 = new Triad<>(C5.multiply(NEG_N1), C1, C4.multiply(NEG_N1));
        Triad<Apfloat> vB50 = new Triad<>(C5.multiply(NEG_N1), C1.multiply(NEG_N1), C4);
        Triad<Apfloat> vB51 = new Triad<>(C5.multiply(NEG_N1), C1.multiply(NEG_N1), C4.multiply(NEG_N1));
        Triad<Apfloat> vB52 = new Triad<>(C4, C5, C1);
        Triad<Apfloat> vB53 = new Triad<>(C4, C5, C1.multiply(NEG_N1));
        Triad<Apfloat> vB54 = new Triad<>(C4, C5.multiply(NEG_N1), C1);
        Triad<Apfloat> vB55 = new Triad<>(C4, C5.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB56 = new Triad<>(C4.multiply(NEG_N1), C5, C1);
        Triad<Apfloat> vB57 = new Triad<>(C4.multiply(NEG_N1), C5, C1.multiply(NEG_N1));
        Triad<Apfloat> vB58 = new Triad<>(C4.multiply(NEG_N1), C5.multiply(NEG_N1), C1);
        Triad<Apfloat> vB59 = new Triad<>(C4.multiply(NEG_N1), C5.multiply(NEG_N1), C1.multiply(NEG_N1));

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = VectorMath.mult(VectorMath.normalize(vB0),  super.getRadius().toString());
        Triad<Apfloat> vC1  = VectorMath.mult(VectorMath.normalize(vB1),  super.getRadius().toString());
        Triad<Apfloat> vC2  = VectorMath.mult(VectorMath.normalize(vB2),  super.getRadius().toString());
        Triad<Apfloat> vC3  = VectorMath.mult(VectorMath.normalize(vB3),  super.getRadius().toString());
        Triad<Apfloat> vC4  = VectorMath.mult(VectorMath.normalize(vB4),  super.getRadius().toString());
        Triad<Apfloat> vC5  = VectorMath.mult(VectorMath.normalize(vB5),  super.getRadius().toString());
        Triad<Apfloat> vC6  = VectorMath.mult(VectorMath.normalize(vB6),  super.getRadius().toString());
        Triad<Apfloat> vC7  = VectorMath.mult(VectorMath.normalize(vB7),  super.getRadius().toString());
        Triad<Apfloat> vC8  = VectorMath.mult(VectorMath.normalize(vB8),  super.getRadius().toString());
        Triad<Apfloat> vC9  = VectorMath.mult(VectorMath.normalize(vB9),  super.getRadius().toString());
        Triad<Apfloat> vC10 = VectorMath.mult(VectorMath.normalize(vB10), super.getRadius().toString());
        Triad<Apfloat> vC11 = VectorMath.mult(VectorMath.normalize(vB11), super.getRadius().toString());
        Triad<Apfloat> vC12 = VectorMath.mult(VectorMath.normalize(vB12), super.getRadius().toString());
        Triad<Apfloat> vC13 = VectorMath.mult(VectorMath.normalize(vB13), super.getRadius().toString());
        Triad<Apfloat> vC14 = VectorMath.mult(VectorMath.normalize(vB14), super.getRadius().toString());
        Triad<Apfloat> vC15 = VectorMath.mult(VectorMath.normalize(vB15), super.getRadius().toString());
        Triad<Apfloat> vC16 = VectorMath.mult(VectorMath.normalize(vB16), super.getRadius().toString());
        Triad<Apfloat> vC17 = VectorMath.mult(VectorMath.normalize(vB17), super.getRadius().toString());
        Triad<Apfloat> vC18 = VectorMath.mult(VectorMath.normalize(vB18), super.getRadius().toString());
        Triad<Apfloat> vC19 = VectorMath.mult(VectorMath.normalize(vB19), super.getRadius().toString());
        Triad<Apfloat> vC20 = VectorMath.mult(VectorMath.normalize(vB20), super.getRadius().toString());
        Triad<Apfloat> vC21 = VectorMath.mult(VectorMath.normalize(vB21), super.getRadius().toString());
        Triad<Apfloat> vC22 = VectorMath.mult(VectorMath.normalize(vB22), super.getRadius().toString());
        Triad<Apfloat> vC23 = VectorMath.mult(VectorMath.normalize(vB23), super.getRadius().toString());
        Triad<Apfloat> vC24 = VectorMath.mult(VectorMath.normalize(vB24), super.getRadius().toString());
        Triad<Apfloat> vC25 = VectorMath.mult(VectorMath.normalize(vB25), super.getRadius().toString());
        Triad<Apfloat> vC26 = VectorMath.mult(VectorMath.normalize(vB26), super.getRadius().toString());
        Triad<Apfloat> vC27 = VectorMath.mult(VectorMath.normalize(vB27), super.getRadius().toString());
        Triad<Apfloat> vC28 = VectorMath.mult(VectorMath.normalize(vB28), super.getRadius().toString());
        Triad<Apfloat> vC29 = VectorMath.mult(VectorMath.normalize(vB29), super.getRadius().toString());
        Triad<Apfloat> vC30 = VectorMath.mult(VectorMath.normalize(vB30), super.getRadius().toString());
        Triad<Apfloat> vC31 = VectorMath.mult(VectorMath.normalize(vB31), super.getRadius().toString());
        Triad<Apfloat> vC32 = VectorMath.mult(VectorMath.normalize(vB32), super.getRadius().toString());
        Triad<Apfloat> vC33 = VectorMath.mult(VectorMath.normalize(vB33), super.getRadius().toString());
        Triad<Apfloat> vC34 = VectorMath.mult(VectorMath.normalize(vB34), super.getRadius().toString());
        Triad<Apfloat> vC35 = VectorMath.mult(VectorMath.normalize(vB35), super.getRadius().toString());
        Triad<Apfloat> vC36 = VectorMath.mult(VectorMath.normalize(vB36), super.getRadius().toString());
        Triad<Apfloat> vC37 = VectorMath.mult(VectorMath.normalize(vB37), super.getRadius().toString());
        Triad<Apfloat> vC38 = VectorMath.mult(VectorMath.normalize(vB38), super.getRadius().toString());
        Triad<Apfloat> vC39 = VectorMath.mult(VectorMath.normalize(vB39), super.getRadius().toString());
        Triad<Apfloat> vC40 = VectorMath.mult(VectorMath.normalize(vB40), super.getRadius().toString());
        Triad<Apfloat> vC41 = VectorMath.mult(VectorMath.normalize(vB41), super.getRadius().toString());
        Triad<Apfloat> vC42 = VectorMath.mult(VectorMath.normalize(vB42), super.getRadius().toString());
        Triad<Apfloat> vC43 = VectorMath.mult(VectorMath.normalize(vB43), super.getRadius().toString());
        Triad<Apfloat> vC44 = VectorMath.mult(VectorMath.normalize(vB44), super.getRadius().toString());
        Triad<Apfloat> vC45 = VectorMath.mult(VectorMath.normalize(vB45), super.getRadius().toString());
        Triad<Apfloat> vC46 = VectorMath.mult(VectorMath.normalize(vB46), super.getRadius().toString());
        Triad<Apfloat> vC47 = VectorMath.mult(VectorMath.normalize(vB47), super.getRadius().toString());
        Triad<Apfloat> vC48 = VectorMath.mult(VectorMath.normalize(vB48), super.getRadius().toString());
        Triad<Apfloat> vC49 = VectorMath.mult(VectorMath.normalize(vB49), super.getRadius().toString());
        Triad<Apfloat> vC50 = VectorMath.mult(VectorMath.normalize(vB50), super.getRadius().toString());
        Triad<Apfloat> vC51 = VectorMath.mult(VectorMath.normalize(vB51), super.getRadius().toString());
        Triad<Apfloat> vC52 = VectorMath.mult(VectorMath.normalize(vB52), super.getRadius().toString());
        Triad<Apfloat> vC53 = VectorMath.mult(VectorMath.normalize(vB53), super.getRadius().toString());
        Triad<Apfloat> vC54 = VectorMath.mult(VectorMath.normalize(vB54), super.getRadius().toString());
        Triad<Apfloat> vC55 = VectorMath.mult(VectorMath.normalize(vB55), super.getRadius().toString());
        Triad<Apfloat> vC56 = VectorMath.mult(VectorMath.normalize(vB56), super.getRadius().toString());
        Triad<Apfloat> vC57 = VectorMath.mult(VectorMath.normalize(vB57), super.getRadius().toString());
        Triad<Apfloat> vC58 = VectorMath.mult(VectorMath.normalize(vB58), super.getRadius().toString());
        Triad<Apfloat> vC59 = VectorMath.mult(VectorMath.normalize(vB59), super.getRadius().toString());

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC0, vC2, vC18, vC42, vC38, vC14);
        Triad<Apfloat> hex0_norm = VectorMath.normalHex(vC0, vC2, vC18, vC42, vC38, vC14, true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC1, vC3, vC17, vC41, vC37, vC13);
        Triad<Apfloat> hex1_norm = VectorMath.normalHex(vC1, vC3, vC17, vC41, vC37, vC13, true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC2, vC0, vC12, vC36, vC40, vC16);
        Triad<Apfloat> hex2_norm = VectorMath.normalHex(vC2, vC0, vC12, vC36, vC40, vC16, true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC3, vC1, vC15, vC39, vC43, vC19);
        Triad<Apfloat> hex3_norm = VectorMath.normalHex(vC3, vC1, vC15, vC39, vC43, vC19, true);
        Hexad<Tuple<Apfloat>> hex4 = new Hexad<>(vC4, vC5, vC23, vC47, vC45, vC21);
        Triad<Apfloat> hex4_norm = VectorMath.normalHex(vC4, vC5, vC23, vC47, vC45, vC21, true);
        Hexad<Tuple<Apfloat>> hex5 = new Hexad<>(vC5, vC4, vC20, vC44, vC46, vC22);
        Triad<Apfloat> hex5_norm = VectorMath.normalHex(vC5, vC4, vC20, vC44, vC46, vC22, true);
        Hexad<Tuple<Apfloat>> hex6 = new Hexad<>(vC6, vC7, vC26, vC50, vC48, vC24);
        Triad<Apfloat> hex6_norm = VectorMath.normalHex(vC6, vC7, vC26, vC50, vC48, vC24, true);
        Hexad<Tuple<Apfloat>> hex7 = new Hexad<>(vC7, vC6, vC25, vC49, vC51, vC27);
        Triad<Apfloat> hex7_norm = VectorMath.normalHex(vC7, vC6, vC25, vC49, vC51, vC27, true);
        Hexad<Tuple<Apfloat>> hex8 = new Hexad<>(vC8, vC9, vC33, vC57, vC56, vC32);
        Triad<Apfloat> hex8_norm = VectorMath.normalHex(vC8, vC9, vC33, vC57, vC56, vC32, true);
        Hexad<Tuple<Apfloat>> hex9 = new Hexad<>(vC9, vC8, vC28, vC52, vC53, vC29);
        Triad<Apfloat> hex9_norm = VectorMath.normalHex(vC9, vC8, vC28, vC52, vC53, vC29, true);
        Hexad<Tuple<Apfloat>> hex10 = new Hexad<>(vC10, vC11, vC31, vC55, vC54, vC30);
        Triad<Apfloat> hex10_norm = VectorMath.normalHex(vC10, vC11, vC31, vC55, vC54, vC30, true);
        Hexad<Tuple<Apfloat>> hex11 = new Hexad<>(vC11, vC10, vC34, vC58, vC59, vC35);
        Triad<Apfloat> hex11_norm = VectorMath.normalHex(vC11, vC10, vC34, vC58, vC59, vC35, true);
        Hexad<Tuple<Apfloat>> hex12 = new Hexad<>(vC12, vC44, vC20, vC52, vC28, vC36);
        Triad<Apfloat> hex12_norm = VectorMath.normalHex(vC12, vC44, vC20, vC52, vC28, vC36, true);
        Hexad<Tuple<Apfloat>> hex13 = new Hexad<>(vC13, vC37, vC29, vC53, vC21, vC45);
        Triad<Apfloat> hex13_norm = VectorMath.normalHex(vC13, vC37, vC29, vC53, vC21, vC45, true);
        Hexad<Tuple<Apfloat>> hex14 = new Hexad<>(vC14, vC38, vC30, vC54, vC22, vC46);
        Triad<Apfloat> hex14_norm = VectorMath.normalHex(vC14, vC38, vC30, vC54, vC22, vC46, true);
        Hexad<Tuple<Apfloat>> hex15 = new Hexad<>(vC15, vC47, vC23, vC55, vC31, vC39);
        Triad<Apfloat> hex15_norm = VectorMath.normalHex(vC15, vC47, vC23, vC55, vC31, vC39, true);
        Hexad<Tuple<Apfloat>> hex16 = new Hexad<>(vC16, vC40, vC32, vC56, vC24, vC48);
        Triad<Apfloat> hex16_norm = VectorMath.normalHex(vC16, vC40, vC32, vC56, vC24, vC48, true);
        Hexad<Tuple<Apfloat>> hex17 = new Hexad<>(vC17, vC49, vC25, vC57, vC33, vC41);
        Triad<Apfloat> hex17_norm = VectorMath.normalHex(vC17, vC49, vC25, vC57, vC33, vC41, true);
        Hexad<Tuple<Apfloat>> hex18 = new Hexad<>(vC18, vC50, vC26, vC58, vC34, vC42);
        Triad<Apfloat> hex18_norm = VectorMath.normalHex(vC18, vC50, vC26, vC58, vC34, vC42, true);
        Hexad<Tuple<Apfloat>> hex19 = new Hexad<>(vC19, vC43, vC35, vC59, vC27, vC51);
        Triad<Apfloat> hex19_norm = VectorMath.normalHex(vC19, vC43, vC35, vC59, vC27, vC51, true);
        faces_hex = new Icosad<>(
                hex0, hex1, hex2, hex3,
                hex4, hex5, hex6, hex7,
                hex8, hex9, hex10, hex11,
                hex12, hex13, hex14, hex15,
                hex16, hex17, hex18, hex19
        );
        face_norms_hex = new Icosad<>(
                hex0_norm, hex1_norm, hex2_norm, hex3_norm,
                hex4_norm, hex5_norm, hex6_norm, hex7_norm,
                hex8_norm, hex9_norm, hex10_norm, hex11_norm,
                hex12_norm, hex13_norm, hex14_norm, hex15_norm,
                hex16_norm, hex17_norm, hex18_norm, hex19_norm

        );

        // ==== PENTAGONAL FACES ====
        Pentad<Tuple<Apfloat>> pnt0 = new Pentad<>(vC0, vC14, vC46, vC44, vC12);
        Triad<Apfloat> pnt0_norm = VectorMath.normalPent(vC0, vC14, vC46, vC44, vC12, true);
        Pentad<Tuple<Apfloat>> pnt1 = new Pentad<>(vC1, vC13, vC45, vC47, vC15);
        Triad<Apfloat> pnt1_norm = VectorMath.normalPent(vC1, vC13, vC45, vC47, vC15, true);
        Pentad<Tuple<Apfloat>> pnt2 = new Pentad<>(vC2, vC16, vC48, vC50, vC18);
        Triad<Apfloat> pnt2_norm = VectorMath.normalPent(vC2, vC16, vC48, vC50, vC18, true);
        Pentad<Tuple<Apfloat>> pnt3 = new Pentad<>(vC3, vC19, vC51, vC49, vC17);
        Triad<Apfloat> pnt3_norm = VectorMath.normalPent(vC3, vC19, vC51, vC49, vC17, true);
        Pentad<Tuple<Apfloat>> pnt4 = new Pentad<>(vC4, vC21, vC53, vC52, vC20);
        Triad<Apfloat> pnt4_norm = VectorMath.normalPent(vC4, vC21, vC53, vC52, vC20, true);
        Pentad<Tuple<Apfloat>> pnt5 = new Pentad<>(vC5, vC22, vC54, vC55, vC23);
        Triad<Apfloat> pnt5_norm = VectorMath.normalPent(vC5, vC22, vC54, vC55, vC23, true);
        Pentad<Tuple<Apfloat>> pnt6 = new Pentad<>(vC6, vC24, vC56, vC57, vC25);
        Triad<Apfloat> pnt6_norm = VectorMath.normalPent(vC6, vC24, vC56, vC57, vC25, true);
        Pentad<Tuple<Apfloat>> pnt7 = new Pentad<>(vC7, vC27, vC59, vC58, vC26);
        Triad<Apfloat> pnt7_norm = VectorMath.normalPent(vC7, vC27, vC59, vC58, vC26, true);
        Pentad<Tuple<Apfloat>> pnt8 = new Pentad<>(vC8, vC32, vC40, vC36, vC28);
        Triad<Apfloat> pnt8_norm = VectorMath.normalPent(vC8, vC32, vC40, vC36, vC28, true);
        Pentad<Tuple<Apfloat>> pnt9 = new Pentad<>(vC9, vC29, vC37, vC41, vC33);
        Triad<Apfloat> pnt9_norm = VectorMath.normalPent(vC9, vC29, vC37, vC41, vC33, true);
        Pentad<Tuple<Apfloat>> pnt10 = new Pentad<>(vC10, vC30, vC38, vC42, vC34);
        Triad<Apfloat> pnt10_norm = VectorMath.normalPent(vC10, vC30, vC38, vC42, vC34, true);
        Pentad<Tuple<Apfloat>> pnt11 = new Pentad<>(vC11, vC35, vC43, vC39, vC31);
        Triad<Apfloat> pnt11_norm = VectorMath.normalPent(vC11, vC35, vC43, vC39, vC31, true);
        faces_pnt = new Dodecad<>(
                pnt0, pnt1, pnt2, pnt3,
                pnt4, pnt5, pnt6, pnt7,
                pnt8, pnt9, pnt10, pnt11
        );
        face_norms_pnt = new Dodecad<>(
                pnt0_norm, pnt1_norm, pnt2_norm, pnt3_norm,
                pnt4_norm, pnt5_norm, pnt6_norm, pnt7_norm,
                pnt8_norm, pnt9_norm, pnt10_norm, pnt11_norm

        );
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < 20; i++) {
            if (i < 12){
                // pentagonal faces
                Pentad<Tuple<Apfloat>> face = (Pentad<Tuple<Apfloat>>) faces_pnt.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_pnt.fetch(i);

                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Apfloat d = VectorMath.dot_prod(norm, m);

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
            Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Apfloat d = VectorMath.dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
