package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.NotNull;
import io.github.noshou.npg.shapes.archimedean.*;

/**
 * Represents a <b>Biscribed <i>levo</i> Pentagonal Icositetrahedron</b>
 * <p> The Biscribed {@link IcositetrahedronPentagonalDextro}.
 * <p> It is the dual of the {@link CuboctahedronSnubLevoBiscribed}. </p>
 * <p> It is the right-handed (<i>dextro</i>) enantiomorph of the {@link IcositetrahedronPentagonal}
 * @see <a href="https://dmccooey.com/polyhedra/BiscribedRpentagonalIcositetrahedron.html">
 *      Biscribed <i>dextro</i> Pentagonal Icositetrahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class IcositetrahedronPentagonalDextroBiscribed extends IcositetrahedronPentagonalDextro {

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N6 = new Apfloat("6", super.precision);
    private final Apfloat N9 = new Apfloat("9", super.precision);
    private final Apfloat N15 = new Apfloat("15", super.precision);
    private final Apfloat N22 = new Apfloat("22", super.precision);
    private final Apfloat N52 = new Apfloat("52", super.precision);
    private final Apfloat N89 = new Apfloat("89", super.precision);
    private final Apfloat SQRT3 = ApfloatMath.sqrt(N3);

    // 4 * sqrt(3) / 6
    private final Apfloat FRAC_N4_mul_SQRT3_over_N6 = N4.multiply(SQRT3).divide(N6);

    // sqrt(3 * (52 * sqrt(3) - 89)) / 6
    private final Apfloat FRAC_SQRT_N3_mul_N52_mul_SQRT3_min_89_over_N6 = ApfloatMath.sqrt(
            N3.multiply(N52.multiply(SQRT3).subtract(N89))
    ).divide(N6);

    // sqrt(6 * (15 * sqrt(3) - 22)) / 6
    private final Apfloat FRAC_SQRT_N6_mul_N15_mul_SQRT3_min_N22_over_N6 = ApfloatMath.sqrt(
            N6.multiply(N15.multiply(SQRT3).subtract(N22))
    ).divide(N6);

    // 9 / 6
    private final Apfloat FRAC_N9_over_N6 = N9.divide(N6);

    // (5*sqrt(3)/6)
    private final Apfloat FRAC_N5_mul_SQRT3_over_N6 = N5.multiply(SQRT3).divide(N6);

    // C0
    private final Apfloat C0 =  FRAC_N4_mul_SQRT3_over_N6.multiply(N2)
                                .subtract(N15.divide(N6))
                                .add(FRAC_SQRT_N3_mul_N52_mul_SQRT3_min_89_over_N6);
    /**
     * @return (8 * sqrt(3) - 15 + sqrt(3 * (52 * sqrt(3) - 89))) / 6
     */
    @Override
    public Apfloat C0() {return C0;}
    /**
     * @return -((8 * sqrt(3) - 15 + sqrt(3 * (52 * sqrt(3) - 89))) / 6)
     */
    @Override
    public Apfloat NEG_C0() {
        return C0().multiply(NEG_N1);
    }

    // C1
    private final Apfloat C1 = SQRT3.divide(N3);
    /**
     * @return sqrt(3) / 3
     */
    @Override
    public Apfloat C1() {return C1;}
    /**
     * @return -(sqrt(3) / 3)
     */
    @Override
    public Apfloat NEG_C1() {
        return C1().multiply(NEG_N1);
    }

    // C2
    private final Apfloat C2 =  FRAC_N9_over_N6
                                .subtract(FRAC_N4_mul_SQRT3_over_N6)
                                .add(FRAC_SQRT_N3_mul_N52_mul_SQRT3_min_89_over_N6);
    /**
     * @return (9 - 4 * sqrt(3) + sqrt(3 * (52 * sqrt(3) - 89))) / 6
     */
    @Override
    public Apfloat C2() {return C2;}
    /**
     * @return -((9 - 4 * sqrt(3) + sqrt(3 * (52 * sqrt(3) - 89))) / 6)
     */
    @Override
    public Apfloat NEG_C2() {
        return C2().multiply(NEG_N1);
    }

    // C3
    private final Apfloat C3 =  FRAC_N5_mul_SQRT3_over_N6
                                .subtract(FRAC_N9_over_N6)
                                .add(FRAC_SQRT_N6_mul_N15_mul_SQRT3_min_N22_over_N6);
    /**
     * @return (5 * sqrt(3) - 9 + sqrt(6 * (15 * sqrt(3) - 22))) / 6
     */
    @Override
    public Apfloat C3() {return C3;}
    /**
     * @return -((5 * sqrt(3) - 9 + sqrt(6 * (15 * sqrt(3) - 22))) / 6)
     */
    @Override
    public Apfloat NEG_C3() {
        return C3().multiply(NEG_N1);
    }

    public IcositetrahedronPentagonalDextroBiscribed (
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
        vB0 = new Triad<>(N0(), N0(), N1);
        vB1 = new Triad<>(N0(), N0(), NEG_N1);
        vB2 = new Triad<>(N1, N0(), N0());
        vB3 = new Triad<>(NEG_N1, N0(), N0());
        vB4 = new Triad<>(N0(), N1, N0());
        vB5 = new Triad<>(N0(), NEG_N1, N0());
        vB6 = new Triad<>(C2(), NEG_C0(), C3());
        vB7 = new Triad<>(C2(), C0(), NEG_C3());
        vB8 = new Triad<>(NEG_C2(), C0(), C3());
        vB9 = new Triad<>(NEG_C2(), NEG_C0(), NEG_C3());
        vB10 = new Triad<>(C3(), NEG_C2(), C0());
        vB11 = new Triad<>(C3(), C2(), NEG_C0());
        vB12 = new Triad<>(NEG_C3(), C2(), C0());
        vB13 = new Triad<>(NEG_C3(), NEG_C2(), NEG_C0());
        vB14 = new Triad<>(C0(), NEG_C3(), C2());
        vB15 = new Triad<>(C0(), C3(), NEG_C2());
        vB16 = new Triad<>(NEG_C0(), C3(), C2());
        vB17 = new Triad<>(NEG_C0(), NEG_C3(), NEG_C2());
        vB18 = new Triad<>(C0(), C2(), C3());
        vB19 = new Triad<>(C0(), NEG_C2(), NEG_C3());
        vB20 = new Triad<>(NEG_C0(), NEG_C2(), C3());
        vB21 = new Triad<>(NEG_C0(), C2(), NEG_C3());
        vB22 = new Triad<>(C3(), C0(), C2());
        vB23 = new Triad<>(C3(), NEG_C0(), NEG_C2());
        vB24 = new Triad<>(NEG_C3(), NEG_C0(), C2());
        vB25 = new Triad<>(NEG_C3(), C0(), NEG_C2());
        vB26 = new Triad<>(C2(), C3(), C0());
        vB27 = new Triad<>(C2(), NEG_C3(), NEG_C0());
        vB28 = new Triad<>(NEG_C2(), NEG_C3(), C0());
        vB29 = new Triad<>(NEG_C2(), C3(), NEG_C0());
        vB30 = new Triad<>(C1(), C1(), C1());
        vB31 = new Triad<>(C1(), C1(), NEG_C1());
        vB32 = new Triad<>(C1(), NEG_C1(), C1());
        vB33 = new Triad<>(C1(), NEG_C1(), NEG_C1());
        vB34 = new Triad<>(NEG_C1(), C1(), C1());
        vB35 = new Triad<>(NEG_C1(), C1(), NEG_C1());
        vB36 = new Triad<>(NEG_C1(), NEG_C1(), C1());
        vB37 = new Triad<>(NEG_C1(), NEG_C1(), NEG_C1());
        setVerts();
    }
}
