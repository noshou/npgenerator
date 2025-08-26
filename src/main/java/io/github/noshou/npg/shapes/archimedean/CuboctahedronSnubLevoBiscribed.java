package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.Polyad;
import io.github.noshou.tuple.Triad;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import io.github.noshou.npg.shapes.catalan.*;

/**
 * Represents a <b>Biscribed <i>levo</i>-Snub Cuboctahedron</b>
 * <p> The biscribed {@link CuboctahedronSnubLevo}.
 * <p> It is the dual of the {@link IcositetrahedronPentagonalDextroBiscribed}. </p>
 * <p> It is the left-handed (<i>levo</i>) enantiomorph of the {@link CuboctahedronSnub}
 * @see <a href="https://dmccooey.com/polyhedra/BiscribedLsnubCube.html">
 *      Biscribed <i>levo</i>-Snub Cuboctahedron (David McCooey)</a>
 */
public class CuboctahedronSnubLevoBiscribed extends CuboctahedronSnubLevo {

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N8 = new Apfloat("8", super.precision);
    private final Apfloat SQRT3 = ApfloatMath.sqrt(N3);
    private final Apfloat N1_add_N8_mul_SQRT3 = N1.add(N8.multiply(SQRT3));


    // C0
    private final Apfloat C0 = (N1.add(N2.multiply(SQRT3.subtract(ApfloatMath.sqrt(N1_add_N8_mul_SQRT3))))).divide(N2);

    /**
     * @return  = (1 + 2 * sqrt(3) - sqrt(1 + 8 * sqrt(3))) / 2
     */
    @Override
    public Apfloat C0() {return C0;}
    /**
     * @return  -(1 + 2 * sqrt(3) - sqrt(1 + 8 * sqrt(3))) / 2
     */
    @Override
    public Apfloat NEG_C0() { return C0().multiply(NEG_N1);}

    // C1
    private final Apfloat C1 = (ApfloatMath.sqrt(N1_add_N8_mul_SQRT3).subtract(N3)).divide(N2);
    /**
     * @return (sqrt(1 + 8 * sqrt(3)) - 3) / 2
     */
    @Override
    public Apfloat C1() {return C1;}

    /**
     * @return  -((sqrt(1 + 8 * sqrt(3)) - 3) / 2)
     */
    @Override
    public Apfloat NEG_C1() {
        return C1().multiply(NEG_N1);
    }

    /**
     * no C2; returns null
     * @return {@code null}
     */
    @Override
    public Apfloat C2(){return null;}

    /**
     * no C2; returns null
     * @return {@code null}
     */
    @Override
    public Apfloat NEG_C2(){return null;}

    public CuboctahedronSnubLevoBiscribed(
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
        vB0  = new Triad<>(C1(),C0(),N1);
        vB1  = new Triad<>(C1(),NEG_C0(), NEG_N1);
        vB2  = new Triad<>(NEG_C1(), NEG_C0(), N1);
        vB3  = new Triad<>(NEG_C1(), C0(),NEG_N1);
        vB4  = new Triad<>(N1,C1(),C0());
        vB5  = new Triad<>(N1,NEG_C1(), NEG_C0());
        vB6  = new Triad<>(NEG_N1,NEG_C1(), C0());
        vB7  = new Triad<>(NEG_N1,C1(),NEG_C0());
        vB8  = new Triad<>(C0(),N1,C1());
        vB9  = new Triad<>(C0(),NEG_N1,NEG_C1());
        vB10 = new Triad<>(NEG_C0(), NEG_N1,C1());
        vB11 = new Triad<>(NEG_C0(), N1,NEG_C1());
        vB12 = new Triad<>(C0(),NEG_C1(), N1);
        vB13 = new Triad<>(C0(),C1(),NEG_N1);
        vB14 = new Triad<>(NEG_C0(), C1(),N1);
        vB15 = new Triad<>(NEG_C0(), NEG_C1(), NEG_N1);
        vB16 = new Triad<>(N1,NEG_C0(), C1());
        vB17 = new Triad<>(N1,C0(),NEG_C1());
        vB18 = new Triad<>(NEG_N1,C0(),C1());
        vB19 = new Triad<>(NEG_N1,NEG_C0(), NEG_C1());
        vB20 = new Triad<>(C1(),NEG_N1,C0());
        vB21 = new Triad<>(C1(),N1,NEG_C0());
        vB22 = new Triad<>(NEG_C1(), N1,C0());
        vB23 = new Triad<>(NEG_C1(), NEG_N1,NEG_C0());

        setVerts();
    }
}
