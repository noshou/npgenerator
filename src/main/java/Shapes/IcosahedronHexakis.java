package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/DisdyakisTriacontahedron.txt>...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class IcosahedronHexakis extends Shape {

    // 120 triangular faces
    private final ArrayList<Triad<Tuple<Apfloat>>> faces = new ArrayList<>();
    private final ArrayList<Triad<Apfloat>> face_norms = new ArrayList<>();

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N9 = new Apfloat("9", super.precision);
    private final Apfloat N10 = new Apfloat("10", super.precision);
    private final Apfloat N11 = new Apfloat("11", super.precision);
    private final Apfloat N15 = new Apfloat("15", super.precision);
    private final Apfloat N22 = new Apfloat("22", super.precision);
    private final Apfloat N44 = new Apfloat("44", super.precision);
    private final Apfloat N105 = new Apfloat("105", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(N5);

    // C0 = 3 * (15 + sqrt(5)) / 44
    private final Apfloat C0 = N3.multiply(N15.add(SQRT5)).divide(N44);

    // C1 = (5 - sqrt(5)) / 2
    private final Apfloat C1 = N5.subtract(SQRT5).divide(N2);

    // C2 = 3 * (5 + 4 * sqrt(5)) / 22
    private final Apfloat C2 = N3.multiply(N5.add(N4.multiply(SQRT5))).divide(N22);

    // C3 = 3 * (5 + sqrt(5)) / 10
    private final Apfloat C3 = N3.multiply(N5.add(SQRT5)).divide(N10);

    // C4 = sqrt(5)
    private final Apfloat C4 = SQRT5;

    // C5 = (75 + 27 * sqrt(5)) / 44
    private final Apfloat C5 = (N105.multiply(SQRT5)).divide(N44);

    // C6 = (15 + 9 * sqrt(5)) / 10
    private final Apfloat C6 = N15.add(N9.multiply(SQRT5)).divide(N10);

    // C7 = (5 + sqrt(5)) / 2
    private final Apfloat C7 = N5.add(SQRT5).divide(N2);

    // C8 = 3 * (5 + 4 * sqrt(5)) / 11
    private final Apfloat C8 = N3.multiply(N9.multiply(SQRT5)).divide(N11);


    public IcosahedronHexakis(
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
        Triad<Apfloat> vC60 = VectorMath.mult(VectorMath.normalize(vB60), super.getRadius().toString());
        Triad<Apfloat> vC61 = VectorMath.mult(VectorMath.normalize(vB61), super.getRadius().toString());

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC18, vC0, vC8);
        faces.add(tri0);
        Triad<Apfloat> tri0_norm = VectorMath.normalTriple(vC18, vC0, vC8, true);
        face_norms.add(tri0_norm);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC18, vC8, vC32);
        faces.add(tri1);
        Triad<Apfloat> tri1_norm = VectorMath.normalTriple(vC18, vC8, vC32, true);
        face_norms.add(tri1_norm);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC18, vC32, vC56);
        faces.add(tri2);
        Triad<Apfloat> tri2_norm = VectorMath.normalTriple(vC18, vC32, vC56, true);
        face_norms.add(tri2_norm);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC18, vC56, vC40);
        faces.add(tri3);
        Triad<Apfloat> tri3_norm = VectorMath.normalTriple(vC18, vC56, vC40, true);
        face_norms.add(tri3_norm);
        Triad<Tuple<Apfloat>> tri4 = new Triad<>(vC18, vC40, vC10);
        faces.add(tri4);
        Triad<Apfloat> tri4_norm = VectorMath.normalTriple(vC18, vC40, vC10, true);
        face_norms.add(tri4_norm);
        Triad<Tuple<Apfloat>> tri5 = new Triad<>(vC18, vC10, vC38);
        faces.add(tri5);
        Triad<Apfloat> tri5_norm = VectorMath.normalTriple(vC18, vC10, vC38, true);
        face_norms.add(tri5_norm);
        Triad<Tuple<Apfloat>> tri6 = new Triad<>(vC18, vC38, vC54);
        faces.add(tri6);
        Triad<Apfloat> tri6_norm = VectorMath.normalTriple(vC18, vC38, vC54, true);
        face_norms.add(tri6_norm);
        Triad<Tuple<Apfloat>> tri7 = new Triad<>(vC18, vC54, vC30);
        faces.add(tri7);
        Triad<Apfloat> tri7_norm = VectorMath.normalTriple(vC18, vC54, vC30, true);
        face_norms.add(tri7_norm);
        Triad<Tuple<Apfloat>> tri8 = new Triad<>(vC18, vC30, vC6);
        faces.add(tri8);
        Triad<Apfloat> tri8_norm = VectorMath.normalTriple(vC18, vC30, vC6, true);
        face_norms.add(tri8_norm);
        Triad<Tuple<Apfloat>> tri9 = new Triad<>(vC18, vC6, vC0);
        faces.add(tri9);
        Triad<Apfloat> tri9_norm = VectorMath.normalTriple(vC18, vC6, vC0, true);
        face_norms.add(tri9_norm);
        Triad<Tuple<Apfloat>> tri10 = new Triad<>(vC19, vC1, vC7);
        faces.add(tri10);
        Triad<Apfloat> tri10_norm = VectorMath.normalTriple(vC19, vC1, vC7, true);
        face_norms.add(tri10_norm);
        Triad<Tuple<Apfloat>> tri11 = new Triad<>(vC19, vC7, vC31);
        faces.add(tri11);
        Triad<Apfloat> tri11_norm = VectorMath.normalTriple(vC19, vC7, vC31, true);
        face_norms.add(tri11_norm);
        Triad<Tuple<Apfloat>> tri12 = new Triad<>(vC19, vC31, vC55);
        faces.add(tri12);
        Triad<Apfloat> tri12_norm = VectorMath.normalTriple(vC19, vC31, vC55, true);
        face_norms.add(tri12_norm);
        Triad<Tuple<Apfloat>> tri13 = new Triad<>(vC19, vC55, vC39);
        faces.add(tri13);
        Triad<Apfloat> tri13_norm = VectorMath.normalTriple(vC19, vC55, vC39, true);
        face_norms.add(tri13_norm);
        Triad<Tuple<Apfloat>> tri14 = new Triad<>(vC19, vC39, vC11);
        faces.add(tri14);
        Triad<Apfloat> tri14_norm = VectorMath.normalTriple(vC19, vC39, vC11, true);
        face_norms.add(tri14_norm);
        Triad<Tuple<Apfloat>> tri15 = new Triad<>(vC19, vC11, vC41);
        faces.add(tri15);
        Triad<Apfloat> tri15_norm = VectorMath.normalTriple(vC19, vC11, vC41, true);
        face_norms.add(tri15_norm);
        Triad<Tuple<Apfloat>> tri16 = new Triad<>(vC19, vC41, vC57);
        faces.add(tri16);
        Triad<Apfloat> tri16_norm = VectorMath.normalTriple(vC19, vC41, vC57, true);
        face_norms.add(tri16_norm);
        Triad<Tuple<Apfloat>> tri17 = new Triad<>(vC19, vC57, vC33);
        faces.add(tri17);
        Triad<Apfloat> tri17_norm = VectorMath.normalTriple(vC19, vC57, vC33, true);
        face_norms.add(tri17_norm);
        Triad<Tuple<Apfloat>> tri18 = new Triad<>(vC19, vC33, vC9);
        faces.add(tri18);
        Triad<Apfloat> tri18_norm = VectorMath.normalTriple(vC19, vC33, vC9, true);
        face_norms.add(tri18_norm);
        Triad<Tuple<Apfloat>> tri19 = new Triad<>(vC19, vC9, vC1);
        faces.add(tri19);
        Triad<Apfloat> tri19_norm = VectorMath.normalTriple(vC19, vC9, vC1, true);
        face_norms.add(tri19_norm);
        Triad<Tuple<Apfloat>> tri20 = new Triad<>(vC20, vC0, vC6);
        faces.add(tri20);
        Triad<Apfloat> tri20_norm = VectorMath.normalTriple(vC20, vC0, vC6, true);
        face_norms.add(tri20_norm);
        Triad<Tuple<Apfloat>> tri21 = new Triad<>(vC20, vC6, vC34);
        faces.add(tri21);
        Triad<Apfloat> tri21_norm = VectorMath.normalTriple(vC20, vC6, vC34, true);
        face_norms.add(tri21_norm);
        Triad<Tuple<Apfloat>> tri22 = new Triad<>(vC20, vC34, vC58);
        faces.add(tri22);
        Triad<Apfloat> tri22_norm = VectorMath.normalTriple(vC20, vC34, vC58, true);
        face_norms.add(tri22_norm);
        Triad<Tuple<Apfloat>> tri23 = new Triad<>(vC20, vC58, vC42);
        faces.add(tri23);
        Triad<Apfloat> tri23_norm = VectorMath.normalTriple(vC20, vC58, vC42, true);
        face_norms.add(tri23_norm);
        Triad<Tuple<Apfloat>> tri24 = new Triad<>(vC20, vC42, vC12);
        faces.add(tri24);
        Triad<Apfloat> tri24_norm = VectorMath.normalTriple(vC20, vC42, vC12, true);
        face_norms.add(tri24_norm);
        Triad<Tuple<Apfloat>> tri25 = new Triad<>(vC20, vC12, vC44);
        faces.add(tri25);
        Triad<Apfloat> tri25_norm = VectorMath.normalTriple(vC20, vC12, vC44, true);
        face_norms.add(tri25_norm);
        Triad<Tuple<Apfloat>> tri26 = new Triad<>(vC20, vC44, vC60);
        faces.add(tri26);
        Triad<Apfloat> tri26_norm = VectorMath.normalTriple(vC20, vC44, vC60, true);
        face_norms.add(tri26_norm);
        Triad<Tuple<Apfloat>> tri27 = new Triad<>(vC20, vC60, vC36);
        faces.add(tri27);
        Triad<Apfloat> tri27_norm = VectorMath.normalTriple(vC20, vC60, vC36, true);
        face_norms.add(tri27_norm);
        Triad<Tuple<Apfloat>> tri28 = new Triad<>(vC20, vC36, vC8);
        faces.add(tri28);
        Triad<Apfloat> tri28_norm = VectorMath.normalTriple(vC20, vC36, vC8, true);
        face_norms.add(tri28_norm);
        Triad<Tuple<Apfloat>> tri29 = new Triad<>(vC20, vC8, vC0);
        faces.add(tri29);
        Triad<Apfloat> tri29_norm = VectorMath.normalTriple(vC20, vC8, vC0, true);
        face_norms.add(tri29_norm);
        Triad<Tuple<Apfloat>> tri30 = new Triad<>(vC21, vC1, vC9);
        faces.add(tri30);
        Triad<Apfloat> tri30_norm = VectorMath.normalTriple(vC21, vC1, vC9, true);
        face_norms.add(tri30_norm);
        Triad<Tuple<Apfloat>> tri31 = new Triad<>(vC21, vC9, vC37);
        faces.add(tri31);
        Triad<Apfloat> tri31_norm = VectorMath.normalTriple(vC21, vC9, vC37, true);
        face_norms.add(tri31_norm);
        Triad<Tuple<Apfloat>> tri32 = new Triad<>(vC21, vC37, vC61);
        faces.add(tri32);
        Triad<Apfloat> tri32_norm = VectorMath.normalTriple(vC21, vC37, vC61, true);
        face_norms.add(tri32_norm);
        Triad<Tuple<Apfloat>> tri33 = new Triad<>(vC21, vC61, vC45);
        faces.add(tri33);
        Triad<Apfloat> tri33_norm = VectorMath.normalTriple(vC21, vC61, vC45, true);
        face_norms.add(tri33_norm);
        Triad<Tuple<Apfloat>> tri34 = new Triad<>(vC21, vC45, vC13);
        faces.add(tri34);
        Triad<Apfloat> tri34_norm = VectorMath.normalTriple(vC21, vC45, vC13, true);
        face_norms.add(tri34_norm);
        Triad<Tuple<Apfloat>> tri35 = new Triad<>(vC21, vC13, vC43);
        faces.add(tri35);
        Triad<Apfloat> tri35_norm = VectorMath.normalTriple(vC21, vC13, vC43, true);
        face_norms.add(tri35_norm);
        Triad<Tuple<Apfloat>> tri36 = new Triad<>(vC21, vC43, vC59);
        faces.add(tri36);
        Triad<Apfloat> tri36_norm = VectorMath.normalTriple(vC21, vC43, vC59, true);
        face_norms.add(tri36_norm);
        Triad<Tuple<Apfloat>> tri37 = new Triad<>(vC21, vC59, vC35);
        faces.add(tri37);
        Triad<Apfloat> tri37_norm = VectorMath.normalTriple(vC21, vC59, vC35, true);
        face_norms.add(tri37_norm);
        Triad<Tuple<Apfloat>> tri38 = new Triad<>(vC21, vC35, vC7);
        faces.add(tri38);
        Triad<Apfloat> tri38_norm = VectorMath.normalTriple(vC21, vC35, vC7, true);
        face_norms.add(tri38_norm);
        Triad<Tuple<Apfloat>> tri39 = new Triad<>(vC21, vC7, vC1);
        faces.add(tri39);
        Triad<Apfloat> tri39_norm = VectorMath.normalTriple(vC21, vC7, vC1, true);
        face_norms.add(tri39_norm);
        Triad<Tuple<Apfloat>> tri40 = new Triad<>(vC22, vC2, vC11);
        faces.add(tri40);
        Triad<Apfloat> tri40_norm = VectorMath.normalTriple(vC22, vC2, vC11, true);
        face_norms.add(tri40_norm);
        Triad<Tuple<Apfloat>> tri41 = new Triad<>(vC22, vC11, vC39);
        faces.add(tri41);
        Triad<Apfloat> tri41_norm = VectorMath.normalTriple(vC22, vC11, vC39, true);
        face_norms.add(tri41_norm);
        Triad<Tuple<Apfloat>> tri42 = new Triad<>(vC22, vC39, vC55);
        faces.add(tri42);
        Triad<Apfloat> tri42_norm = VectorMath.normalTriple(vC22, vC39, vC55, true);
        face_norms.add(tri42_norm);
        Triad<Tuple<Apfloat>> tri43 = new Triad<>(vC22, vC55, vC47);
        faces.add(tri43);
        Triad<Apfloat> tri43_norm = VectorMath.normalTriple(vC22, vC55, vC47, true);
        face_norms.add(tri43_norm);
        Triad<Tuple<Apfloat>> tri44 = new Triad<>(vC22, vC47, vC14);
        faces.add(tri44);
        Triad<Apfloat> tri44_norm = VectorMath.normalTriple(vC22, vC47, vC14, true);
        face_norms.add(tri44_norm);
        Triad<Tuple<Apfloat>> tri45 = new Triad<>(vC22, vC14, vC46);
        faces.add(tri45);
        Triad<Apfloat> tri45_norm = VectorMath.normalTriple(vC22, vC14, vC46, true);
        face_norms.add(tri45_norm);
        Triad<Tuple<Apfloat>> tri46 = new Triad<>(vC22, vC46, vC54);
        faces.add(tri46);
        Triad<Apfloat> tri46_norm = VectorMath.normalTriple(vC22, vC46, vC54, true);
        face_norms.add(tri46_norm);
        Triad<Tuple<Apfloat>> tri47 = new Triad<>(vC22, vC54, vC38);
        faces.add(tri47);
        Triad<Apfloat> tri47_norm = VectorMath.normalTriple(vC22, vC54, vC38, true);
        face_norms.add(tri47_norm);
        Triad<Tuple<Apfloat>> tri48 = new Triad<>(vC22, vC38, vC10);
        faces.add(tri48);
        Triad<Apfloat> tri48_norm = VectorMath.normalTriple(vC22, vC38, vC10, true);
        face_norms.add(tri48_norm);
        Triad<Tuple<Apfloat>> tri49 = new Triad<>(vC22, vC10, vC2);
        faces.add(tri49);
        Triad<Apfloat> tri49_norm = VectorMath.normalTriple(vC22, vC10, vC2, true);
        face_norms.add(tri49_norm);
        Triad<Tuple<Apfloat>> tri50 = new Triad<>(vC23, vC2, vC10);
        faces.add(tri50);
        Triad<Apfloat> tri50_norm = VectorMath.normalTriple(vC23, vC2, vC10, true);
        face_norms.add(tri50_norm);
        Triad<Tuple<Apfloat>> tri51 = new Triad<>(vC23, vC10, vC40);
        faces.add(tri51);
        Triad<Apfloat> tri51_norm = VectorMath.normalTriple(vC23, vC10, vC40, true);
        face_norms.add(tri51_norm);
        Triad<Tuple<Apfloat>> tri52 = new Triad<>(vC23, vC40, vC56);
        faces.add(tri52);
        Triad<Apfloat> tri52_norm = VectorMath.normalTriple(vC23, vC40, vC56, true);
        face_norms.add(tri52_norm);
        Triad<Tuple<Apfloat>> tri53 = new Triad<>(vC23, vC56, vC48);
        faces.add(tri53);
        Triad<Apfloat> tri53_norm = VectorMath.normalTriple(vC23, vC56, vC48, true);
        face_norms.add(tri53_norm);
        Triad<Tuple<Apfloat>> tri54 = new Triad<>(vC23, vC48, vC15);
        faces.add(tri54);
        Triad<Apfloat> tri54_norm = VectorMath.normalTriple(vC23, vC48, vC15, true);
        face_norms.add(tri54_norm);
        Triad<Tuple<Apfloat>> tri55 = new Triad<>(vC23, vC15, vC49);
        faces.add(tri55);
        Triad<Apfloat> tri55_norm = VectorMath.normalTriple(vC23, vC15, vC49, true);
        face_norms.add(tri55_norm);
        Triad<Tuple<Apfloat>> tri56 = new Triad<>(vC23, vC49, vC57);
        faces.add(tri56);
        Triad<Apfloat> tri56_norm = VectorMath.normalTriple(vC23, vC49, vC57, true);
        face_norms.add(tri56_norm);
        Triad<Tuple<Apfloat>> tri57 = new Triad<>(vC23, vC57, vC41);
        faces.add(tri57);
        Triad<Apfloat> tri57_norm = VectorMath.normalTriple(vC23, vC57, vC41, true);
        face_norms.add(tri57_norm);
        Triad<Tuple<Apfloat>> tri58 = new Triad<>(vC23, vC41, vC11);
        faces.add(tri58);
        Triad<Apfloat> tri58_norm = VectorMath.normalTriple(vC23, vC41, vC11, true);
        face_norms.add(tri58_norm);
        Triad<Tuple<Apfloat>> tri59 = new Triad<>(vC23, vC11, vC2);
        faces.add(tri59);
        Triad<Apfloat> tri59_norm = VectorMath.normalTriple(vC23, vC11, vC2, true);
        face_norms.add(tri59_norm);
        Triad<Tuple<Apfloat>> tri60 = new Triad<>(vC24, vC3, vC12);
        faces.add(tri60);
        Triad<Apfloat> tri60_norm = VectorMath.normalTriple(vC24, vC3, vC12, true);
        face_norms.add(tri60_norm);
        Triad<Tuple<Apfloat>> tri61 = new Triad<>(vC24, vC12, vC42);
        faces.add(tri61);
        Triad<Apfloat> tri61_norm = VectorMath.normalTriple(vC24, vC12, vC42, true);
        face_norms.add(tri61_norm);
        Triad<Tuple<Apfloat>> tri62 = new Triad<>(vC24, vC42, vC58);
        faces.add(tri62);
        Triad<Apfloat> tri62_norm = VectorMath.normalTriple(vC24, vC42, vC58, true);
        face_norms.add(tri62_norm);
        Triad<Tuple<Apfloat>> tri63 = new Triad<>(vC24, vC58, vC50);
        faces.add(tri63);
        Triad<Apfloat> tri63_norm = VectorMath.normalTriple(vC24, vC58, vC50, true);
        face_norms.add(tri63_norm);
        Triad<Tuple<Apfloat>> tri64 = new Triad<>(vC24, vC50, vC16);
        faces.add(tri64);
        Triad<Apfloat> tri64_norm = VectorMath.normalTriple(vC24, vC50, vC16, true);
        face_norms.add(tri64_norm);
        Triad<Tuple<Apfloat>> tri65 = new Triad<>(vC24, vC16, vC51);
        faces.add(tri65);
        Triad<Apfloat> tri65_norm = VectorMath.normalTriple(vC24, vC16, vC51, true);
        face_norms.add(tri65_norm);
        Triad<Tuple<Apfloat>> tri66 = new Triad<>(vC24, vC51, vC59);
        faces.add(tri66);
        Triad<Apfloat> tri66_norm = VectorMath.normalTriple(vC24, vC51, vC59, true);
        face_norms.add(tri66_norm);
        Triad<Tuple<Apfloat>> tri67 = new Triad<>(vC24, vC59, vC43);
        faces.add(tri67);
        Triad<Apfloat> tri67_norm = VectorMath.normalTriple(vC24, vC59, vC43, true);
        face_norms.add(tri67_norm);
        Triad<Tuple<Apfloat>> tri68 = new Triad<>(vC24, vC43, vC13);
        faces.add(tri68);
        Triad<Apfloat> tri68_norm = VectorMath.normalTriple(vC24, vC43, vC13, true);
        face_norms.add(tri68_norm);
        Triad<Tuple<Apfloat>> tri69 = new Triad<>(vC24, vC13, vC3);
        faces.add(tri69);
        Triad<Apfloat> tri69_norm = VectorMath.normalTriple(vC24, vC13, vC3, true);
        face_norms.add(tri69_norm);
        Triad<Tuple<Apfloat>> tri70 = new Triad<>(vC25, vC3, vC13);
        faces.add(tri70);
        Triad<Apfloat> tri70_norm = VectorMath.normalTriple(vC25, vC3, vC13, true);
        face_norms.add(tri70_norm);
        Triad<Tuple<Apfloat>> tri71 = new Triad<>(vC25, vC13, vC45);
        faces.add(tri71);
        Triad<Apfloat> tri71_norm = VectorMath.normalTriple(vC25, vC13, vC45, true);
        face_norms.add(tri71_norm);
        Triad<Tuple<Apfloat>> tri72 = new Triad<>(vC25, vC45, vC61);
        faces.add(tri72);
        Triad<Apfloat> tri72_norm = VectorMath.normalTriple(vC25, vC45, vC61, true);
        face_norms.add(tri72_norm);
        Triad<Tuple<Apfloat>> tri73 = new Triad<>(vC25, vC61, vC53);
        faces.add(tri73);
        Triad<Apfloat> tri73_norm = VectorMath.normalTriple(vC25, vC61, vC53, true);
        face_norms.add(tri73_norm);
        Triad<Tuple<Apfloat>> tri74 = new Triad<>(vC25, vC53, vC17);
        faces.add(tri74);
        Triad<Apfloat> tri74_norm = VectorMath.normalTriple(vC25, vC53, vC17, true);
        face_norms.add(tri74_norm);
        Triad<Tuple<Apfloat>> tri75 = new Triad<>(vC25, vC17, vC52);
        faces.add(tri75);
        Triad<Apfloat> tri75_norm = VectorMath.normalTriple(vC25, vC17, vC52, true);
        face_norms.add(tri75_norm);
        Triad<Tuple<Apfloat>> tri76 = new Triad<>(vC25, vC52, vC60);
        faces.add(tri76);
        Triad<Apfloat> tri76_norm = VectorMath.normalTriple(vC25, vC52, vC60, true);
        face_norms.add(tri76_norm);
        Triad<Tuple<Apfloat>> tri77 = new Triad<>(vC25, vC60, vC44);
        faces.add(tri77);
        Triad<Apfloat> tri77_norm = VectorMath.normalTriple(vC25, vC60, vC44, true);
        face_norms.add(tri77_norm);
        Triad<Tuple<Apfloat>> tri78 = new Triad<>(vC25, vC44, vC12);
        faces.add(tri78);
        Triad<Apfloat> tri78_norm = VectorMath.normalTriple(vC25, vC44, vC12, true);
        face_norms.add(tri78_norm);
        Triad<Tuple<Apfloat>> tri79 = new Triad<>(vC25, vC12, vC3);
        faces.add(tri79);
        Triad<Apfloat> tri79_norm = VectorMath.normalTriple(vC25, vC12, vC3, true);
        face_norms.add(tri79_norm);
        Triad<Tuple<Apfloat>> tri80 = new Triad<>(vC26, vC4, vC16);
        faces.add(tri80);
        Triad<Apfloat> tri80_norm = VectorMath.normalTriple(vC26, vC4, vC16, true);
        face_norms.add(tri80_norm);
        Triad<Tuple<Apfloat>> tri81 = new Triad<>(vC26, vC16, vC50);
        faces.add(tri81);
        Triad<Apfloat> tri81_norm = VectorMath.normalTriple(vC26, vC16, vC50, true);
        face_norms.add(tri81_norm);
        Triad<Tuple<Apfloat>> tri82 = new Triad<>(vC26, vC50, vC58);
        faces.add(tri82);
        Triad<Apfloat> tri82_norm = VectorMath.normalTriple(vC26, vC50, vC58, true);
        face_norms.add(tri82_norm);
        Triad<Tuple<Apfloat>> tri83 = new Triad<>(vC26, vC58, vC34);
        faces.add(tri83);
        Triad<Apfloat> tri83_norm = VectorMath.normalTriple(vC26, vC58, vC34, true);
        face_norms.add(tri83_norm);
        Triad<Tuple<Apfloat>> tri84 = new Triad<>(vC26, vC34, vC6);
        faces.add(tri84);
        Triad<Apfloat> tri84_norm = VectorMath.normalTriple(vC26, vC34, vC6, true);
        face_norms.add(tri84_norm);
        Triad<Tuple<Apfloat>> tri85 = new Triad<>(vC26, vC6, vC30);
        faces.add(tri85);
        Triad<Apfloat> tri85_norm = VectorMath.normalTriple(vC26, vC6, vC30, true);
        face_norms.add(tri85_norm);
        Triad<Tuple<Apfloat>> tri86 = new Triad<>(vC26, vC30, vC54);
        faces.add(tri86);
        Triad<Apfloat> tri86_norm = VectorMath.normalTriple(vC26, vC30, vC54, true);
        face_norms.add(tri86_norm);
        Triad<Tuple<Apfloat>> tri87 = new Triad<>(vC26, vC54, vC46);
        faces.add(tri87);
        Triad<Apfloat> tri87_norm = VectorMath.normalTriple(vC26, vC54, vC46, true);
        face_norms.add(tri87_norm);
        Triad<Tuple<Apfloat>> tri88 = new Triad<>(vC26, vC46, vC14);
        faces.add(tri88);
        Triad<Apfloat> tri88_norm = VectorMath.normalTriple(vC26, vC46, vC14, true);
        face_norms.add(tri88_norm);
        Triad<Tuple<Apfloat>> tri89 = new Triad<>(vC26, vC14, vC4);
        faces.add(tri89);
        Triad<Apfloat> tri89_norm = VectorMath.normalTriple(vC26, vC14, vC4, true);
        face_norms.add(tri89_norm);
        Triad<Tuple<Apfloat>> tri90 = new Triad<>(vC27, vC4, vC14);
        faces.add(tri90);
        Triad<Apfloat> tri90_norm = VectorMath.normalTriple(vC27, vC4, vC14, true);
        face_norms.add(tri90_norm);
        Triad<Tuple<Apfloat>> tri91 = new Triad<>(vC27, vC14, vC47);
        faces.add(tri91);
        Triad<Apfloat> tri91_norm = VectorMath.normalTriple(vC27, vC14, vC47, true);
        face_norms.add(tri91_norm);
        Triad<Tuple<Apfloat>> tri92 = new Triad<>(vC27, vC47, vC55);
        faces.add(tri92);
        Triad<Apfloat> tri92_norm = VectorMath.normalTriple(vC27, vC47, vC55, true);
        face_norms.add(tri92_norm);
        Triad<Tuple<Apfloat>> tri93 = new Triad<>(vC27, vC55, vC31);
        faces.add(tri93);
        Triad<Apfloat> tri93_norm = VectorMath.normalTriple(vC27, vC55, vC31, true);
        face_norms.add(tri93_norm);
        Triad<Tuple<Apfloat>> tri94 = new Triad<>(vC27, vC31, vC7);
        faces.add(tri94);
        Triad<Apfloat> tri94_norm = VectorMath.normalTriple(vC27, vC31, vC7, true);
        face_norms.add(tri94_norm);
        Triad<Tuple<Apfloat>> tri95 = new Triad<>(vC27, vC7, vC35);
        faces.add(tri95);
        Triad<Apfloat> tri95_norm = VectorMath.normalTriple(vC27, vC7, vC35, true);
        face_norms.add(tri95_norm);
        Triad<Tuple<Apfloat>> tri96 = new Triad<>(vC27, vC35, vC59);
        faces.add(tri96);
        Triad<Apfloat> tri96_norm = VectorMath.normalTriple(vC27, vC35, vC59, true);
        face_norms.add(tri96_norm);
        Triad<Tuple<Apfloat>> tri97 = new Triad<>(vC27, vC59, vC51);
        faces.add(tri97);
        Triad<Apfloat> tri97_norm = VectorMath.normalTriple(vC27, vC59, vC51, true);
        face_norms.add(tri97_norm);
        Triad<Tuple<Apfloat>> tri98 = new Triad<>(vC27, vC51, vC16);
        faces.add(tri98);
        Triad<Apfloat> tri98_norm = VectorMath.normalTriple(vC27, vC51, vC16, true);
        face_norms.add(tri98_norm);
        Triad<Tuple<Apfloat>> tri99 = new Triad<>(vC27, vC16, vC4);
        faces.add(tri99);
        Triad<Apfloat> tri99_norm = VectorMath.normalTriple(vC27, vC16, vC4, true);
        face_norms.add(tri99_norm);
        Triad<Tuple<Apfloat>> tri100 = new Triad<>(vC28, vC5, vC15);
        faces.add(tri100);
        Triad<Apfloat> tri100_norm = VectorMath.normalTriple(vC28, vC5, vC15, true);
        face_norms.add(tri100_norm);
        Triad<Tuple<Apfloat>> tri101 = new Triad<>(vC28, vC15, vC48);
        faces.add(tri101);
        Triad<Apfloat> tri101_norm = VectorMath.normalTriple(vC28, vC15, vC48, true);
        face_norms.add(tri101_norm);
        Triad<Tuple<Apfloat>> tri102 = new Triad<>(vC28, vC48, vC56);
        faces.add(tri102);
        Triad<Apfloat> tri102_norm = VectorMath.normalTriple(vC28, vC48, vC56, true);
        face_norms.add(tri102_norm);
        Triad<Tuple<Apfloat>> tri103 = new Triad<>(vC28, vC56, vC32);
        faces.add(tri103);
        Triad<Apfloat> tri103_norm = VectorMath.normalTriple(vC28, vC56, vC32, true);
        face_norms.add(tri103_norm);
        Triad<Tuple<Apfloat>> tri104 = new Triad<>(vC28, vC32, vC8);
        faces.add(tri104);
        Triad<Apfloat> tri104_norm = VectorMath.normalTriple(vC28, vC32, vC8, true);
        face_norms.add(tri104_norm);
        Triad<Tuple<Apfloat>> tri105 = new Triad<>(vC28, vC8, vC36);
        faces.add(tri105);
        Triad<Apfloat> tri105_norm = VectorMath.normalTriple(vC28, vC8, vC36, true);
        face_norms.add(tri105_norm);
        Triad<Tuple<Apfloat>> tri106 = new Triad<>(vC28, vC36, vC60);
        faces.add(tri106);
        Triad<Apfloat> tri106_norm = VectorMath.normalTriple(vC28, vC36, vC60, true);
        face_norms.add(tri106_norm);
        Triad<Tuple<Apfloat>> tri107 = new Triad<>(vC28, vC60, vC52);
        faces.add(tri107);
        Triad<Apfloat> tri107_norm = VectorMath.normalTriple(vC28, vC60, vC52, true);
        face_norms.add(tri107_norm);
        Triad<Tuple<Apfloat>> tri108 = new Triad<>(vC28, vC52, vC17);
        faces.add(tri108);
        Triad<Apfloat> tri108_norm = VectorMath.normalTriple(vC28, vC52, vC17, true);
        face_norms.add(tri108_norm);
        Triad<Tuple<Apfloat>> tri109 = new Triad<>(vC28, vC17, vC5);
        faces.add(tri109);
        Triad<Apfloat> tri109_norm = VectorMath.normalTriple(vC28, vC17, vC5, true);
        face_norms.add(tri109_norm);
        Triad<Tuple<Apfloat>> tri110 = new Triad<>(vC29, vC5, vC17);
        faces.add(tri110);
        Triad<Apfloat> tri110_norm = VectorMath.normalTriple(vC29, vC5, vC17, true);
        face_norms.add(tri110_norm);
        Triad<Tuple<Apfloat>> tri111 = new Triad<>(vC29, vC17, vC53);
        faces.add(tri111);
        Triad<Apfloat> tri111_norm = VectorMath.normalTriple(vC29, vC17, vC53, true);
        face_norms.add(tri111_norm);
        Triad<Tuple<Apfloat>> tri112 = new Triad<>(vC29, vC53, vC61);
        faces.add(tri112);
        Triad<Apfloat> tri112_norm = VectorMath.normalTriple(vC29, vC53, vC61, true);
        face_norms.add(tri112_norm);
        Triad<Tuple<Apfloat>> tri113 = new Triad<>(vC29, vC61, vC37);
        faces.add(tri113);
        Triad<Apfloat> tri113_norm = VectorMath.normalTriple(vC29, vC61, vC37, true);
        face_norms.add(tri113_norm);
        Triad<Tuple<Apfloat>> tri114 = new Triad<>(vC29, vC37, vC9);
        faces.add(tri114);
        Triad<Apfloat> tri114_norm = VectorMath.normalTriple(vC29, vC37, vC9, true);
        face_norms.add(tri114_norm);
        Triad<Tuple<Apfloat>> tri115 = new Triad<>(vC29, vC9, vC33);
        faces.add(tri115);
        Triad<Apfloat> tri115_norm = VectorMath.normalTriple(vC29, vC9, vC33, true);
        face_norms.add(tri115_norm);
        Triad<Tuple<Apfloat>> tri116 = new Triad<>(vC29, vC33, vC57);
        faces.add(tri116);
        Triad<Apfloat> tri116_norm = VectorMath.normalTriple(vC29, vC33, vC57, true);
        face_norms.add(tri116_norm);
        Triad<Tuple<Apfloat>> tri117 = new Triad<>(vC29, vC57, vC49);
        faces.add(tri117);
        Triad<Apfloat> tri117_norm = VectorMath.normalTriple(vC29, vC57, vC49, true);
        face_norms.add(tri117_norm);
        Triad<Tuple<Apfloat>> tri118 = new Triad<>(vC29, vC49, vC15);
        faces.add(tri118);
        Triad<Apfloat> tri118_norm = VectorMath.normalTriple(vC29, vC49, vC15, true);
        face_norms.add(tri118_norm);
        Triad<Tuple<Apfloat>> tri119 = new Triad<>(vC29, vC15, vC5);
        faces.add(tri119);
        Triad<Apfloat> tri119_norm = VectorMath.normalTriple(vC29, vC15, vC5, true);
        face_norms.add(tri119_norm);
    }

    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        // check if each point lies within bounds of each face
        for (int i = 0; i < faces.size(); i++) {
            // get vertices of vertA
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces.get(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x - vertA_z)
            Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> face_norm = (Triad<Apfloat>) face_norms.get(i);
            Apfloat d = VectorMath.dot_prod(face_norm, m);

            // if d < 0,  point lies behind face (within bounds)
            // if d == 0, point lies on face (within bounds)
            // if d > 0,  point lies in front of face (out of bounds)
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
