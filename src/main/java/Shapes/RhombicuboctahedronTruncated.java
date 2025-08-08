package Shapes;

import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;
import Lattice.*;
import Atom.*;
import java.util.*;

/**
 * {@link:<a href="https://dmccooey.com/polyhedra/CanonicalTruncatedRhombicuboctahedron.html">...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class RhombicuboctahedronTruncated extends Shape {

    // 24 Quadrilateral Faces
    private final Polyad<Tuple<Tuple<Apfloat>>> faces_sqr;
    private final Polyad<Tuple<Apfloat>> face_norms_sqr;

    // 8 Hexagonal Faces
    private final Octad<Tuple<Tuple<Apfloat>>> faces_hex;
    private final Octad<Tuple<Apfloat>> face_norms_hex;

    // 18 octagonal Faces
    private final Octakaidecad<Tuple<Tuple<Apfloat>>> faces_oct;
    private final Octakaidecad<Tuple<Apfloat>> face_norms_oct;

    // Apfloat Constants
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat SQRT2 = ApfloatMath.sqrt(N2);
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N8 = new Apfloat("8", super.precision);
    private final Apfloat N12 = new Apfloat("12", super.precision);
    private final Apfloat N18 = new Apfloat("18", super.precision);
    private final Apfloat N24 = new Apfloat("24", super.precision);
    private final Apfloat N28 = new Apfloat("28", super.precision);
    private final Apfloat N29 = new Apfloat("29", super.precision);
    private final Apfloat N30 = new Apfloat("30", super.precision);
    private final Apfloat N33 = new Apfloat("33", super.precision);
    private final Apfloat N46 = new Apfloat("46", super.precision);
    private final Apfloat N50 = new Apfloat("50", super.precision);
    private final Apfloat N103 = new Apfloat("103", super.precision); // 103
    private final Apfloat N112 = new Apfloat("112", super.precision); // 112
    private final Apfloat N117 = new Apfloat("117", super.precision); // 117
    private final Apfloat N134 = new Apfloat("134", super.precision); // 134
    private final Apfloat N213 = new Apfloat("213", super.precision);
    private final Apfloat N2327 = new Apfloat("2327", super.precision); // 2327
    private final Apfloat N1644 = new Apfloat("1644", super.precision);   // 1644


    public RhombicuboctahedronTruncated(
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

        // ==== CONSTANTS ====
        // s0 = sqrt(3 * (50 - 33 * sqrt(2) - sqrt(2 * (2327 - 1644 * sqrt(2))))) / 12
        Apfloat s0_00 = N2327.subtract(N1644.multiply(SQRT2));
        Apfloat s0_01 = ApfloatMath.sqrt(N2.multiply(s0_00));
        Apfloat s0_02 = N50.subtract(N33.multiply(SQRT2));
        Apfloat s0_03 = s0_02.subtract(s0_01);
        Apfloat s0 = ApfloatMath.sqrt(N3.multiply(s0_03)).divide(N12);

        // s1 = sqrt(3 * (8 - 3 * sqrt(2) + sqrt(2 * (29 - 18 * sqrt(2))))) / 12
        Apfloat s1_00 = N29.subtract(N18.multiply(SQRT2));
        Apfloat s1_01 = ApfloatMath.sqrt(N2.multiply(s1_00));
        Apfloat s1_02 = N8.subtract(N3.multiply(SQRT2));
        Apfloat s1_03 = s1_02.add(s1_01);
        Apfloat s1 = ApfloatMath.sqrt(N3.multiply(s1_03)).divide(N12);

        // s2 = sqrt(3 * (18 + sqrt(2) - sqrt(2 * (103 - 24 * sqrt(2))))) / 12
        Apfloat s2_00 = N103.subtract(N24.multiply(SQRT2));
        Apfloat s2_01 = ApfloatMath.sqrt(N2.multiply(s2_00));
        Apfloat s2_02 = N18.add(SQRT2);
        Apfloat s2_03 = s2_02.subtract(s2_01);
        Apfloat s2 = ApfloatMath.sqrt(N3.multiply(s2_03)).divide(N12);

        // s3 = sqrt(28 + 5 * sqrt(2) + sqrt(2 * (117 - 46 * sqrt(2)))) / 12
        Apfloat s3_00 = N117.subtract(N46.multiply(SQRT2));
        Apfloat s3_01 = ApfloatMath.sqrt(N2.multiply(s3_00));
        Apfloat s3_02 = N28.add(N5.multiply(SQRT2));
        Apfloat s3_03 = s3_02.add(s3_01);
        Apfloat s3 = ApfloatMath.sqrt(s3_03).divide(N12);

        // s4 = sqrt(112 - sqrt(2) - sqrt(2 * (213 + 134 * sqrt(2)))) / 12
        Apfloat s4_00 = N213.add(N134.multiply(SQRT2));
        Apfloat s4_01 = ApfloatMath.sqrt(N2.multiply(s4_00));
        Apfloat s4_02 = N112.subtract(SQRT2);
        Apfloat s4_03 = s4_02.subtract(s4_01);
        Apfloat s4 = ApfloatMath.sqrt(s4_03).divide(N12);

        // s5 = sqrt(3 * (30 - sqrt(2) + sqrt(2 * (103 - 24 * sqrt(2))))) / 12
        Apfloat s5_00 = N103.subtract(N24.multiply(SQRT2));
        Apfloat s5_01 = ApfloatMath.sqrt(N2.multiply(s5_00));
        Apfloat s5_02 = N30.subtract(SQRT2);
        Apfloat s5_03 = s5_02.add(s5_01);
        Apfloat s5 = ApfloatMath.sqrt(N3.multiply(s5_03)).divide(N12);

        // ==== CANONICAL BASIS VERTICES ====
        Triad<Apfloat> vB0  = new Triad<>(s2, s0, s5);                                                      // ( s2,  s0,  s5)
        Triad<Apfloat> vB1  = new Triad<>(s2, s0, s5.multiply(NEG_N1));                                     // ( s2,  s0, -s5)
        Triad<Apfloat> vB2  = new Triad<>(s2, s0.multiply(NEG_N1), s5);                                     // ( s2, -s0,  s5)
        Triad<Apfloat> vB3  = new Triad<>(s2, s0.multiply(NEG_N1), s5.multiply(NEG_N1));                    // ( s2, -s0, -s5)
        Triad<Apfloat> vB4  = new Triad<>(s2.multiply(NEG_N1), s0, s5);                                     // (-s2,  s0,  s5)
        Triad<Apfloat> vB5  = new Triad<>(s2.multiply(NEG_N1), s0, s5.multiply(NEG_N1));                    // (-s2,  s0, -s5)
        Triad<Apfloat> vB6  = new Triad<>(s2.multiply(NEG_N1), s0.multiply(NEG_N1), s5);                    // (-s2, -s0,  s5)
        Triad<Apfloat> vB7  = new Triad<>(s2.multiply(NEG_N1), s0.multiply(NEG_N1), s5.multiply(NEG_N1));   // (-s2, -s0, -s5)
        Triad<Apfloat> vB8  = new Triad<>(s5, s2, s0);                                                      // ( s5,  s2,  s0)
        Triad<Apfloat> vB9  = new Triad<>(s5, s2, s0.multiply(NEG_N1));                                     // ( s5,  s2, -s0)
        Triad<Apfloat> vB10 = new Triad<>(s5, s2.multiply(NEG_N1), s0);                                     // ( s5, -s2,  s0)
        Triad<Apfloat> vB11 = new Triad<>(s5, s2.multiply(NEG_N1), s0.multiply(NEG_N1));                    // ( s5, -s2, -s0)
        Triad<Apfloat> vB12 = new Triad<>(s5.multiply(NEG_N1), s2, s0);                                     // (-s5,  s2,  s0)
        Triad<Apfloat> vB13 = new Triad<>(s5.multiply(NEG_N1), s2, s0.multiply(NEG_N1));                    // (-s5,  s2, -s0)
        Triad<Apfloat> vB14 = new Triad<>(s5.multiply(NEG_N1), s2.multiply(NEG_N1), s0);                    // (-s5, -s2,  s0)
        Triad<Apfloat> vB15 = new Triad<>(s5.multiply(NEG_N1), s2.multiply(NEG_N1), s0.multiply(NEG_N1));   // (-s5, -s2, -s0)
        Triad<Apfloat> vB16 = new Triad<>(s0, s5, s2);                                                      // ( s0,  s5,  s2)
        Triad<Apfloat> vB17 = new Triad<>(s0, s5, s2.multiply(NEG_N1));                                     // ( s0,  s5, -s2)
        Triad<Apfloat> vB18 = new Triad<>(s0, s5.multiply(NEG_N1), s2);                                     // ( s0, -s5,  s2)
        Triad<Apfloat> vB19 = new Triad<>(s0, s5.multiply(NEG_N1), s2.multiply(NEG_N1));                    // ( s0, -s5, -s2)
        Triad<Apfloat> vB20 = new Triad<>(s0.multiply(NEG_N1), s5, s2);                                     // (-s0,  s5,  s2)
        Triad<Apfloat> vB21 = new Triad<>(s0.multiply(NEG_N1), s5, s2.multiply(NEG_N1));                    // (-s0,  s5, -s2)
        Triad<Apfloat> vB22 = new Triad<>(s0.multiply(NEG_N1), s5.multiply(NEG_N1), s2);                    // (-s0, -s5,  s2)
        Triad<Apfloat> vB23 = new Triad<>(s0.multiply(NEG_N1), s5.multiply(NEG_N1), s2.multiply(NEG_N1));   // (-s0, -s5, -s2)
        Triad<Apfloat> vB24 = new Triad<>(s0, s2, s5);                                                      // ( s0,  s2,  s5)
        Triad<Apfloat> vB25 = new Triad<>(s0, s2, s5.multiply(NEG_N1));                                     // ( s0,  s2, -s5)
        Triad<Apfloat> vB26 = new Triad<>(s0, s2.multiply(NEG_N1), s5);                                     // ( s0, -s2,  s5)
        Triad<Apfloat> vB27 = new Triad<>(s0, s2.multiply(NEG_N1), s5.multiply(NEG_N1));                    // ( s0, -s2, -s5)
        Triad<Apfloat> vB28 = new Triad<>(s0.multiply(NEG_N1), s2, s5);                                     // (-s0,  s2,  s5)
        Triad<Apfloat> vB29 = new Triad<>(s0.multiply(NEG_N1), s2, s5.multiply(NEG_N1));                    // (-s0,  s2, -s5)
        Triad<Apfloat> vB30 = new Triad<>(s0.multiply(NEG_N1), s2.multiply(NEG_N1), s5);                    // (-s0, -s2,  s5)
        Triad<Apfloat> vB31 = new Triad<>(s0.multiply(NEG_N1), s2.multiply(NEG_N1), s5.multiply(NEG_N1));   // (-s0, -s2, -s5)
        Triad<Apfloat> vB32 = new Triad<>(s5, s0, s2);                                                      // ( s5,  s0,  s2)
        Triad<Apfloat> vB33 = new Triad<>(s5, s0, s2.multiply(NEG_N1));                                     // ( s5,  s0, -s2)
        Triad<Apfloat> vB34 = new Triad<>(s5, s0.multiply(NEG_N1), s2);                                     // ( s5, -s0,  s2)
        Triad<Apfloat> vB35 = new Triad<>(s5, s0.multiply(NEG_N1), s2.multiply(NEG_N1));                    // ( s5, -s0, -s2)
        Triad<Apfloat> vB36 = new Triad<>(s5.multiply(NEG_N1), s0, s2);                                     // (-s5,  s0,  s2)
        Triad<Apfloat> vB37 = new Triad<>(s5.multiply(NEG_N1), s0, s2.multiply(NEG_N1));                    // (-s5,  s0, -s2)
        Triad<Apfloat> vB38 = new Triad<>(s5.multiply(NEG_N1), s0.multiply(NEG_N1), s2);                    // (-s5, -s0,  s2)
        Triad<Apfloat> vB39 = new Triad<>(s5.multiply(NEG_N1), s0.multiply(NEG_N1), s2.multiply(NEG_N1));   // (-s5, -s0, -s2)
        Triad<Apfloat> vB40 = new Triad<>(s2, s5, s0);                                                      // ( s2,  s5,  s0)
        Triad<Apfloat> vB41 = new Triad<>(s2, s5, s0.multiply(NEG_N1));                                     // ( s2,  s5, -s0)
        Triad<Apfloat> vB42 = new Triad<>(s2, s5.multiply(NEG_N1), s0);                                     // ( s2, -s5,  s0)
        Triad<Apfloat> vB43 = new Triad<>(s2, s5.multiply(NEG_N1), s0.multiply(NEG_N1));                    // ( s2, -s5, -s0)
        Triad<Apfloat> vB44 = new Triad<>(s2.multiply(NEG_N1), s5, s0);                                     // (-s2,  s5,  s0)
        Triad<Apfloat> vB45 = new Triad<>(s2.multiply(NEG_N1), s5, s0.multiply(NEG_N1));                    // (-s2,  s5, -s0)
        Triad<Apfloat> vB46 = new Triad<>(s2.multiply(NEG_N1), s5.multiply(NEG_N1), s0);                    // (-s2, -s5,  s0)
        Triad<Apfloat> vB47 = new Triad<>(s2.multiply(NEG_N1), s5.multiply(NEG_N1), s0.multiply(NEG_N1));   // (-s2, -s5, -s0)
        Triad<Apfloat> vB48 = new Triad<>(s3, s1, s4);                                                      // ( s3,  s1,  s4)
        Triad<Apfloat> vB49 = new Triad<>(s3, s1, s4.multiply(NEG_N1));                                     // ( s3,  s1, -s4)
        Triad<Apfloat> vB50 = new Triad<>(s3, s1.multiply(NEG_N1), s4);                                     // ( s3, -s1,  s4)
        Triad<Apfloat> vB51 = new Triad<>(s3, s1.multiply(NEG_N1), s4.multiply(NEG_N1));                    // ( s3, -s1, -s4)
        Triad<Apfloat> vB52 = new Triad<>(s3.multiply(NEG_N1), s1, s4);                                     // (-s3,  s1,  s4)
        Triad<Apfloat> vB53 = new Triad<>(s3.multiply(NEG_N1), s1, s4.multiply(NEG_N1));                    // (-s3,  s1, -s4)
        Triad<Apfloat> vB54 = new Triad<>(s3.multiply(NEG_N1), s1.multiply(NEG_N1), s4);                    // (-s3, -s1,  s4)
        Triad<Apfloat> vB55 = new Triad<>(s3.multiply(NEG_N1), s1.multiply(NEG_N1), s4.multiply(NEG_N1));   // (-s3, -s1, -s4)
        Triad<Apfloat> vB56 = new Triad<>(s4, s3, s1);                                                      // ( s4,  s3,  s1)
        Triad<Apfloat> vB57 = new Triad<>(s4, s3, s1.multiply(NEG_N1));                                     // ( s4,  s3, -s1)
        Triad<Apfloat> vB58 = new Triad<>(s4, s3.multiply(NEG_N1), s1);                                     // ( s4, -s3,  s1)
        Triad<Apfloat> vB59 = new Triad<>(s4, s3.multiply(NEG_N1), s1.multiply(NEG_N1));                    // ( s4, -s3, -s1)
        Triad<Apfloat> vB60 = new Triad<>(s4.multiply(NEG_N1), s3, s1);                                     // (-s4,  s3,  s1)
        Triad<Apfloat> vB61 = new Triad<>(s4.multiply(NEG_N1), s3, s1.multiply(NEG_N1));                    // (-s4,  s3, -s1)
        Triad<Apfloat> vB62 = new Triad<>(s4.multiply(NEG_N1), s3.multiply(NEG_N1), s1);                    // (-s4, -s3,  s1)
        Triad<Apfloat> vB63 = new Triad<>(s4.multiply(NEG_N1), s3.multiply(NEG_N1), s1.multiply(NEG_N1));   // (-s4, -s3, -s1)
        Triad<Apfloat> vB64 = new Triad<>(s1, s4, s3);                                                      // ( s1,  s4,  s3)
        Triad<Apfloat> vB65 = new Triad<>(s1, s4, s3.multiply(NEG_N1));                                     // ( s1,  s4, -s3)
        Triad<Apfloat> vB66 = new Triad<>(s1, s4.multiply(NEG_N1), s3);                                     // ( s1, -s4,  s3)
        Triad<Apfloat> vB67 = new Triad<>(s1, s4.multiply(NEG_N1), s3.multiply(NEG_N1));                    // ( s1, -s4, -s3)
        Triad<Apfloat> vB68 = new Triad<>(s1.multiply(NEG_N1), s4, s3);                                     // (-s1,  s4,  s3)
        Triad<Apfloat> vB69 = new Triad<>(s1.multiply(NEG_N1), s4, s3.multiply(NEG_N1));                    // (-s1,  s4, -s3)
        Triad<Apfloat> vB70 = new Triad<>(s1.multiply(NEG_N1), s4.multiply(NEG_N1), s3);                    // (-s1, -s4,  s3)
        Triad<Apfloat> vB71 = new Triad<>(s1.multiply(NEG_N1), s4.multiply(NEG_N1), s3.multiply(NEG_N1));   // (-s1, -s4, -s3)
        Triad<Apfloat> vB72 = new Triad<>(s1, s3, s4);                                                      // ( s1,  s3,  s4)
        Triad<Apfloat> vB73 = new Triad<>(s1, s3, s4.multiply(NEG_N1));                                     // ( s1,  s3, -s4)
        Triad<Apfloat> vB74 = new Triad<>(s1, s3.multiply(NEG_N1), s4);                                     // ( s1, -s3,  s4)
        Triad<Apfloat> vB75 = new Triad<>(s1, s3.multiply(NEG_N1), s4.multiply(NEG_N1));                    // ( s1, -s3, -s4)
        Triad<Apfloat> vB76 = new Triad<>(s1.multiply(NEG_N1), s3, s4);                                     // (-s1,  s3,  s4)
        Triad<Apfloat> vB77 = new Triad<>(s1.multiply(NEG_N1), s3, s4.multiply(NEG_N1));                    // (-s1,  s3, -s4)
        Triad<Apfloat> vB78 = new Triad<>(s1.multiply(NEG_N1), s3.multiply(NEG_N1), s4);                    // (-s1, -s3,  s4)
        Triad<Apfloat> vB79 = new Triad<>(s1.multiply(NEG_N1), s3.multiply(NEG_N1), s4.multiply(NEG_N1));   // (-s1, -s3, -s4)
        Triad<Apfloat> vB80 = new Triad<>(s4, s1, s3);                                                      // ( s4,  s1,  s3)
        Triad<Apfloat> vB81 = new Triad<>(s4, s1, s3.multiply(NEG_N1));                                     // ( s4,  s1, -s3)
        Triad<Apfloat> vB82 = new Triad<>(s4, s1.multiply(NEG_N1), s3);                                     // ( s4, -s1,  s3)
        Triad<Apfloat> vB83 = new Triad<>(s4, s1.multiply(NEG_N1), s3.multiply(NEG_N1));                    // ( s4, -s1, -s3)
        Triad<Apfloat> vB84 = new Triad<>(s4.multiply(NEG_N1), s1, s3);                                     // (-s4,  s1,  s3)
        Triad<Apfloat> vB85 = new Triad<>(s4.multiply(NEG_N1), s1, s3.multiply(NEG_N1));                    // (-s4,  s1, -s3)
        Triad<Apfloat> vB86 = new Triad<>(s4.multiply(NEG_N1), s1.multiply(NEG_N1), s3);                    // (-s4, -s1,  s3)
        Triad<Apfloat> vB87 = new Triad<>(s4.multiply(NEG_N1), s1.multiply(NEG_N1), s3.multiply(NEG_N1));   // (-s4, -s1, -s3)
        Triad<Apfloat> vB88 = new Triad<>(s3, s4, s1);                                                      // ( s3,  s4,  s1)
        Triad<Apfloat> vB89 = new Triad<>(s3, s4, s1.multiply(NEG_N1));                                     // ( s3,  s4, -s1)
        Triad<Apfloat> vB90 = new Triad<>(s3, s4.multiply(NEG_N1), s1);                                     // ( s3, -s4,  s1)
        Triad<Apfloat> vB91 = new Triad<>(s3, s4.multiply(NEG_N1), s1.multiply(NEG_N1));                    // ( s3, -s4, -s1)
        Triad<Apfloat> vB92 = new Triad<>(s3.multiply(NEG_N1), s4, s1);                                     // (-s3,  s4,  s1)
        Triad<Apfloat> vB93 = new Triad<>(s3.multiply(NEG_N1), s4, s1.multiply(NEG_N1));                    // (-s3,  s4, -s1)
        Triad<Apfloat> vB94 = new Triad<>(s3.multiply(NEG_N1), s4.multiply(NEG_N1), s1);                    // (-s3, -s4,  s1)
        Triad<Apfloat> vB95 = new Triad<>(s3.multiply(NEG_N1), s4.multiply(NEG_N1), s1.multiply(NEG_N1));   // (-s3, -s4, -s1)

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = VectorMath.mult(VectorMath.normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1  = VectorMath.mult(VectorMath.normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2  = VectorMath.mult(VectorMath.normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3  = VectorMath.mult(VectorMath.normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4  = VectorMath.mult(VectorMath.normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5  = VectorMath.mult(VectorMath.normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6  = VectorMath.mult(VectorMath.normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7  = VectorMath.mult(VectorMath.normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8  = VectorMath.mult(VectorMath.normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9  = VectorMath.mult(VectorMath.normalize(vB9), super.getRadius().toString());
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
        Triad<Apfloat> vC62 = VectorMath.mult(VectorMath.normalize(vB62), super.getRadius().toString());
        Triad<Apfloat> vC63 = VectorMath.mult(VectorMath.normalize(vB63), super.getRadius().toString());
        Triad<Apfloat> vC64 = VectorMath.mult(VectorMath.normalize(vB64), super.getRadius().toString());
        Triad<Apfloat> vC65 = VectorMath.mult(VectorMath.normalize(vB65), super.getRadius().toString());
        Triad<Apfloat> vC66 = VectorMath.mult(VectorMath.normalize(vB66), super.getRadius().toString());
        Triad<Apfloat> vC67 = VectorMath.mult(VectorMath.normalize(vB67), super.getRadius().toString());
        Triad<Apfloat> vC68 = VectorMath.mult(VectorMath.normalize(vB68), super.getRadius().toString());
        Triad<Apfloat> vC69 = VectorMath.mult(VectorMath.normalize(vB69), super.getRadius().toString());
        Triad<Apfloat> vC70 = VectorMath.mult(VectorMath.normalize(vB70), super.getRadius().toString());
        Triad<Apfloat> vC71 = VectorMath.mult(VectorMath.normalize(vB71), super.getRadius().toString());
        Triad<Apfloat> vC72 = VectorMath.mult(VectorMath.normalize(vB72), super.getRadius().toString());
        Triad<Apfloat> vC73 = VectorMath.mult(VectorMath.normalize(vB73), super.getRadius().toString());
        Triad<Apfloat> vC74 = VectorMath.mult(VectorMath.normalize(vB74), super.getRadius().toString());
        Triad<Apfloat> vC75 = VectorMath.mult(VectorMath.normalize(vB75), super.getRadius().toString());
        Triad<Apfloat> vC76 = VectorMath.mult(VectorMath.normalize(vB76), super.getRadius().toString());
        Triad<Apfloat> vC77 = VectorMath.mult(VectorMath.normalize(vB77), super.getRadius().toString());
        Triad<Apfloat> vC78 = VectorMath.mult(VectorMath.normalize(vB78), super.getRadius().toString());
        Triad<Apfloat> vC79 = VectorMath.mult(VectorMath.normalize(vB79), super.getRadius().toString());
        Triad<Apfloat> vC80 = VectorMath.mult(VectorMath.normalize(vB80), super.getRadius().toString());
        Triad<Apfloat> vC81 = VectorMath.mult(VectorMath.normalize(vB81), super.getRadius().toString());
        Triad<Apfloat> vC82 = VectorMath.mult(VectorMath.normalize(vB82), super.getRadius().toString());
        Triad<Apfloat> vC83 = VectorMath.mult(VectorMath.normalize(vB83), super.getRadius().toString());
        Triad<Apfloat> vC84 = VectorMath.mult(VectorMath.normalize(vB84), super.getRadius().toString());
        Triad<Apfloat> vC85 = VectorMath.mult(VectorMath.normalize(vB85), super.getRadius().toString());
        Triad<Apfloat> vC86 = VectorMath.mult(VectorMath.normalize(vB86), super.getRadius().toString());
        Triad<Apfloat> vC87 = VectorMath.mult(VectorMath.normalize(vB87), super.getRadius().toString());
        Triad<Apfloat> vC88 = VectorMath.mult(VectorMath.normalize(vB88), super.getRadius().toString());
        Triad<Apfloat> vC89 = VectorMath.mult(VectorMath.normalize(vB89), super.getRadius().toString());
        Triad<Apfloat> vC90 = VectorMath.mult(VectorMath.normalize(vB90), super.getRadius().toString());
        Triad<Apfloat> vC91 = VectorMath.mult(VectorMath.normalize(vB91), super.getRadius().toString());
        Triad<Apfloat> vC92 = VectorMath.mult(VectorMath.normalize(vB92), super.getRadius().toString());
        Triad<Apfloat> vC93 = VectorMath.mult(VectorMath.normalize(vB93), super.getRadius().toString());
        Triad<Apfloat> vC94 = VectorMath.mult(VectorMath.normalize(vB94), super.getRadius().toString());
        Triad<Apfloat> vC95 = VectorMath.mult(VectorMath.normalize(vB95), super.getRadius().toString());

        // ==== OCTAGONAL FACES ====
        Octad<Tuple<Apfloat>> oct0 = new Octad<>(vC0, vC24, vC28, vC4, vC6, vC30, vC26, vC2);
        Triad<Apfloat> oct0_norm = VectorMath.normalOct(vC0, vC24, vC28, vC4, vC6, vC30, vC26, vC2, true);
        Octad<Tuple<Apfloat>> oct1 = new Octad<>(vC1, vC3, vC27, vC31, vC7, vC5, vC29, vC25);
        Triad<Apfloat> oct1_norm = VectorMath.normalOct(vC1, vC3, vC27, vC31, vC7, vC5, vC29, vC25, true);
        Octad<Tuple<Apfloat>> oct2 = new Octad<>(vC8, vC32, vC34, vC10, vC11, vC35, vC33, vC9);
        Triad<Apfloat> oct2_norm = VectorMath.normalOct(vC8, vC32, vC34, vC10, vC11, vC35, vC33, vC9, true);
        Octad<Tuple<Apfloat>> oct3 = new Octad<>(vC12, vC13, vC37, vC39, vC15, vC14, vC38, vC36);
        Triad<Apfloat> oct3_norm = VectorMath.normalOct(vC12, vC13, vC37, vC39, vC15, vC14, vC38, vC36, true);
        Octad<Tuple<Apfloat>> oct4 = new Octad<>(vC16, vC40, vC41, vC17, vC21, vC45, vC44, vC20);
        Triad<Apfloat> oct4_norm = VectorMath.normalOct(vC16, vC40, vC41, vC17, vC21, vC45, vC44, vC20, true);
        Octad<Tuple<Apfloat>> oct5 = new Octad<>(vC18, vC22, vC46, vC47, vC23, vC19, vC43, vC42);
        Triad<Apfloat> oct5_norm = VectorMath.normalOct(vC18, vC22, vC46, vC47, vC23, vC19, vC43, vC42, true);
        Octad<Tuple<Apfloat>> oct6 = new Octad<>(vC0, vC2, vC50, vC82, vC34, vC32, vC80, vC48);
        Triad<Apfloat> oct6_norm = VectorMath.normalOct(vC0, vC2, vC50, vC82, vC34, vC32, vC80, vC48, true);
        Octad<Tuple<Apfloat>> oct7 = new Octad<>(vC1, vC49, vC81, vC33, vC35, vC83, vC51, vC3);
        Triad<Apfloat> oct7_norm = VectorMath.normalOct(vC1, vC49, vC81, vC33, vC35, vC83, vC51, vC3, true);
        Octad<Tuple<Apfloat>> oct8 = new Octad<>(vC4, vC52, vC84, vC36, vC38, vC86, vC54, vC6);
        Triad<Apfloat> oct8_norm = VectorMath.normalOct(vC4, vC52, vC84, vC36, vC38, vC86, vC54, vC6, true);
        Octad<Tuple<Apfloat>> oct9 = new Octad<>(vC5, vC7, vC55, vC87, vC39, vC37, vC85, vC53);
        Triad<Apfloat> oct9_norm = VectorMath.normalOct(vC5, vC7, vC55, vC87, vC39, vC37, vC85, vC53, true);
        Octad<Tuple<Apfloat>> oct10 = new Octad<>(vC8, vC9, vC57, vC89, vC41, vC40, vC88, vC56);
        Triad<Apfloat> oct10_norm = VectorMath.normalOct(vC8, vC9, vC57, vC89, vC41, vC40, vC88, vC56, true);
        Octad<Tuple<Apfloat>> oct11 = new Octad<>(vC10, vC58, vC90, vC42, vC43, vC91, vC59, vC11);
        Triad<Apfloat> oct11_norm = VectorMath.normalOct(vC10, vC58, vC90, vC42, vC43, vC91, vC59, vC11, true);
        Octad<Tuple<Apfloat>> oct12 = new Octad<>(vC12, vC60, vC92, vC44, vC45, vC93, vC61, vC13);
        Triad<Apfloat> oct12_norm = VectorMath.normalOct(vC12, vC60, vC92, vC44, vC45, vC93, vC61, vC13, true);
        Octad<Tuple<Apfloat>> oct13 = new Octad<>(vC14, vC15, vC63, vC95, vC47, vC46, vC94, vC62);
        Triad<Apfloat> oct13_norm = VectorMath.normalOct(vC14, vC15, vC63, vC95, vC47, vC46, vC94, vC62, true);
        Octad<Tuple<Apfloat>> oct14 = new Octad<>(vC16, vC20, vC68, vC76, vC28, vC24, vC72, vC64);
        Triad<Apfloat> oct14_norm = VectorMath.normalOct(vC16, vC20, vC68, vC76, vC28, vC24, vC72, vC64, true);
        Octad<Tuple<Apfloat>> oct15 = new Octad<>(vC17, vC65, vC73, vC25, vC29, vC77, vC69, vC21);
        Triad<Apfloat> oct15_norm = VectorMath.normalOct(vC17, vC65, vC73, vC25, vC29, vC77, vC69, vC21, true);
        Octad<Tuple<Apfloat>> oct16 = new Octad<>(vC18, vC66, vC74, vC26, vC30, vC78, vC70, vC22);
        Triad<Apfloat> oct16_norm = VectorMath.normalOct(vC18, vC66, vC74, vC26, vC30, vC78, vC70, vC22, true);
        Octad<Tuple<Apfloat>> oct17 = new Octad<>(vC19, vC23, vC71, vC79, vC31, vC27, vC75, vC67);
        Triad<Apfloat> oct17_norm = VectorMath.normalOct(vC19, vC23, vC71, vC79, vC31, vC27, vC75, vC67, true);
        faces_oct = new Octakaidecad<>(
                oct0, oct1, oct2, oct3,
                oct4, oct5, oct6, oct7,
                oct8, oct9, oct10, oct11,
                oct12, oct13, oct14, oct15,
                oct16, oct17
        );
        face_norms_oct = new Octakaidecad<>(
                oct0_norm, oct1_norm, oct2_norm, oct3_norm,
                oct4_norm, oct5_norm, oct6_norm, oct7_norm,
                oct8_norm, oct9_norm, oct10_norm, oct11_norm,
                oct12_norm, oct13_norm, oct14_norm, oct15_norm,
                oct16_norm, oct17_norm
        );

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC48, vC80, vC56, vC88, vC64, vC72);
        Triad<Apfloat> hex0_norm = VectorMath.normalHex(vC48, vC80, vC56, vC88, vC64, vC72, true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC49, vC73, vC65, vC89, vC57, vC81);
        Triad<Apfloat> hex1_norm = VectorMath.normalHex(vC49, vC73, vC65, vC89, vC57, vC81, true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC50, vC74, vC66, vC90, vC58, vC82);
        Triad<Apfloat> hex2_norm = VectorMath.normalHex(vC50, vC74, vC66, vC90, vC58, vC82, true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC51, vC83, vC59, vC91, vC67, vC75);
        Triad<Apfloat> hex3_norm = VectorMath.normalHex(vC51, vC83, vC59, vC91, vC67, vC75, true);
        Hexad<Tuple<Apfloat>> hex4 = new Hexad<>(vC52, vC76, vC68, vC92, vC60, vC84);
        Triad<Apfloat> hex4_norm = VectorMath.normalHex(vC52, vC76, vC68, vC92, vC60, vC84, true);
        Hexad<Tuple<Apfloat>> hex5 = new Hexad<>(vC53, vC85, vC61, vC93, vC69, vC77);
        Triad<Apfloat> hex5_norm = VectorMath.normalHex(vC53, vC85, vC61, vC93, vC69, vC77, true);
        Hexad<Tuple<Apfloat>> hex6 = new Hexad<>(vC54, vC86, vC62, vC94, vC70, vC78);
        Triad<Apfloat> hex6_norm = VectorMath.normalHex(vC54, vC86, vC62, vC94, vC70, vC78, true);
        Hexad<Tuple<Apfloat>> hex7 = new Hexad<>(vC55, vC79, vC71, vC95, vC63, vC87);
        Triad<Apfloat> hex7_norm = VectorMath.normalHex(vC55, vC79, vC71, vC95, vC63, vC87, true);
        faces_hex = new Octad<>(hex0,hex1,hex2,hex3,hex4,hex5,hex6,hex7);
        face_norms_hex = new Octad<>(hex0_norm,hex1_norm,hex2_norm,hex3_norm,hex4_norm,hex5_norm,hex6_norm,hex7_norm);

        // ==== QUADRILATERAL FACES ====

        List<Tetrad<Tuple<Apfloat>>> sqr = new ArrayList<>();
        List<Triad<Apfloat>> nrm = new ArrayList<>();

        Tetrad<Tuple<Apfloat>> sqr0 = new Tetrad<>(vC0, vC48, vC72, vC24);
        Triad<Apfloat> sqr0_norm = VectorMath.normalQuad(vC0, vC48, vC72, vC24, true);
        sqr.add(sqr0);
        nrm.add(sqr0_norm);

        Tetrad<Tuple<Apfloat>> sqr1 = new Tetrad<>(vC1, vC25, vC73, vC49);
        Triad<Apfloat> sqr1_norm = VectorMath.normalQuad(vC1, vC25, vC73, vC49, true);
        sqr.add(sqr1);
        nrm.add(sqr1_norm);

        Tetrad<Tuple<Apfloat>> sqr2 = new Tetrad<>(vC2, vC26, vC74, vC50);
        Triad<Apfloat> sqr2_norm = VectorMath.normalQuad(vC2, vC26, vC74, vC50, true);
        sqr.add(sqr2);
        nrm.add(sqr2_norm);

        Tetrad<Tuple<Apfloat>> sqr3 = new Tetrad<>(vC3, vC51, vC75, vC27);
        Triad<Apfloat> sqr3_norm = VectorMath.normalQuad(vC3, vC51, vC75, vC27, true);
        sqr.add(sqr3);
        nrm.add(sqr3_norm);

        Tetrad<Tuple<Apfloat>> sqr4 = new Tetrad<>(vC4, vC28, vC76, vC52);
        Triad<Apfloat> sqr4_norm = VectorMath.normalQuad(vC4, vC28, vC76, vC52, true);
        sqr.add(sqr4);
        nrm.add(sqr4_norm);

        Tetrad<Tuple<Apfloat>> sqr5 = new Tetrad<>(vC5, vC53, vC77, vC29);
        Triad<Apfloat> sqr5_norm = VectorMath.normalQuad(vC5, vC53, vC77, vC29, true);
        sqr.add(sqr5);
        nrm.add(sqr5_norm);

        Tetrad<Tuple<Apfloat>> sqr6 = new Tetrad<>(vC6, vC54, vC78, vC30);
        Triad<Apfloat> sqr6_norm = VectorMath.normalQuad(vC6, vC54, vC78, vC30, true);
        sqr.add(sqr6);
        nrm.add(sqr6_norm);

        Tetrad<Tuple<Apfloat>> sqr7 = new Tetrad<>(vC7, vC31, vC79, vC55);
        Triad<Apfloat> sqr7_norm = VectorMath.normalQuad(vC7, vC31, vC79, vC55, true);
        sqr.add(sqr7);
        nrm.add(sqr7_norm);

        Tetrad<Tuple<Apfloat>> sqr8 = new Tetrad<>(vC8, vC56, vC80, vC32);
        Triad<Apfloat> sqr8_norm = VectorMath.normalQuad(vC8, vC56, vC80, vC32, true);
        sqr.add(sqr8);
        nrm.add(sqr8_norm);

        Tetrad<Tuple<Apfloat>> sqr9 = new Tetrad<>(vC9, vC33, vC81, vC57);
        Triad<Apfloat> sqr9_norm = VectorMath.normalQuad(vC9, vC33, vC81, vC57, true);
        sqr.add(sqr9);
        nrm.add(sqr9_norm);

        Tetrad<Tuple<Apfloat>> sqr10 = new Tetrad<>(vC10, vC34, vC82, vC58);
        Triad<Apfloat> sqr10_norm = VectorMath.normalQuad(vC10, vC34, vC82, vC58, true);
        sqr.add(sqr10);
        nrm.add(sqr10_norm);

        Tetrad<Tuple<Apfloat>> sqr11 = new Tetrad<>(vC11, vC59, vC83, vC35);
        Triad<Apfloat> sqr11_norm = VectorMath.normalQuad(vC11, vC59, vC83, vC35, true);
        sqr.add(sqr11);
        nrm.add(sqr11_norm);

        Tetrad<Tuple<Apfloat>> sqr12 = new Tetrad<>(vC12, vC36, vC84, vC60);
        Triad<Apfloat> sqr12_norm = VectorMath.normalQuad(vC12, vC36, vC84, vC60, true);
        sqr.add(sqr12);
        nrm.add(sqr12_norm);

        Tetrad<Tuple<Apfloat>> sqr13 = new Tetrad<>(vC13, vC61, vC85, vC37);
        Triad<Apfloat> sqr13_norm = VectorMath.normalQuad(vC13, vC61, vC85, vC37, true);
        sqr.add(sqr13);
        nrm.add(sqr13_norm);

        Tetrad<Tuple<Apfloat>> sqr14 = new Tetrad<>(vC14, vC62, vC86, vC38);
        Triad<Apfloat> sqr14_norm = VectorMath.normalQuad(vC14, vC62, vC86, vC38, true);
        sqr.add(sqr14);
        nrm.add(sqr14_norm);

        Tetrad<Tuple<Apfloat>> sqr15 = new Tetrad<>(vC15, vC39, vC87, vC63);
        Triad<Apfloat> sqr15_norm = VectorMath.normalQuad(vC15, vC39, vC87, vC63, true);
        sqr.add(sqr15);
        nrm.add(sqr15_norm);

        Tetrad<Tuple<Apfloat>> sqr16 = new Tetrad<>(vC16, vC64, vC88, vC40);
        Triad<Apfloat> sqr16_norm = VectorMath.normalQuad(vC16, vC64, vC88, vC40, true);
        sqr.add(sqr16);
        nrm.add(sqr16_norm);

        Tetrad<Tuple<Apfloat>> sqr17 = new Tetrad<>(vC17, vC41, vC89, vC65);
        Triad<Apfloat> sqr17_norm = VectorMath.normalQuad(vC17, vC41, vC89, vC65, true);
        sqr.add(sqr17);
        nrm.add(sqr17_norm);

        Tetrad<Tuple<Apfloat>> sqr18 = new Tetrad<>(vC18, vC42, vC90, vC66);
        Triad<Apfloat> sqr18_norm = VectorMath.normalQuad(vC18, vC42, vC90, vC66, true);
        sqr.add(sqr18);
        nrm.add(sqr18_norm);

        Tetrad<Tuple<Apfloat>> sqr19 = new Tetrad<>(vC19, vC67, vC91, vC43);
        Triad<Apfloat> sqr19_norm = VectorMath.normalQuad(vC19, vC67, vC91, vC43, true);
        sqr.add(sqr19);
        nrm.add(sqr19_norm);

        Tetrad<Tuple<Apfloat>> sqr20 = new Tetrad<>(vC20, vC44, vC92, vC68);
        Triad<Apfloat> sqr20_norm = VectorMath.normalQuad(vC20, vC44, vC92, vC68, true);
        sqr.add(sqr20);
        nrm.add(sqr20_norm);

        Tetrad<Tuple<Apfloat>> sqr21 = new Tetrad<>(vC21, vC69, vC93, vC45);
        Triad<Apfloat> sqr21_norm = VectorMath.normalQuad(vC21, vC69, vC93, vC45, true);
        sqr.add(sqr21);
        nrm.add(sqr21_norm);

        Tetrad<Tuple<Apfloat>> sqr22 = new Tetrad<>(vC22, vC70, vC94, vC46);
        Triad<Apfloat> sqr22_norm = VectorMath.normalQuad(vC22, vC70, vC94, vC46, true);
        sqr.add(sqr22);
        nrm.add(sqr22_norm);

        Tetrad<Tuple<Apfloat>> sqr23 = new Tetrad<>(vC23, vC47, vC95, vC71);
        Triad<Apfloat> sqr23_norm = VectorMath.normalQuad(vC23, vC47, vC95, vC71, true);
        sqr.add(sqr23);
        nrm.add(sqr23_norm);

        Tetrad<Tuple<Apfloat>>[] face = sqr.toArray(Tetrad[]::new);
        Triad<Apfloat>[] nrms = nrm.toArray(Triad[]::new);
        faces_sqr = new Polyad<>(face);
        face_norms_sqr = new Polyad<>(nrms);
    }


    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_sqr.fetchSize(); i++){
            // === Check hexagonal faces ===
            if (i < face_norms_hex.fetchSize()) {
                Hexad<Tuple<Apfloat>> face = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_hex.fetch(i);
                Apfloat d = VectorMath.dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check octagonal faces ===
            if (i < faces_oct.fetchSize()) {
                Octad<Tuple<Apfloat>> face = (Octad<Tuple<Apfloat>>) faces_oct.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_oct.fetch(i);
                Apfloat d = VectorMath.dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check square faces ===
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
            Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_sqr.fetch(i);
            Apfloat d = VectorMath.dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }

        return true;
    }
}
