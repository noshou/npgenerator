package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.Polyad;
import io.github.noshou.tuple.Tetrad;
import io.github.noshou.tuple.Triad;
import io.github.noshou.tuple.Tuple;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;

/**
 * Represents a <b>Canonical <i>levo</i>-Snub Cuboctahedron</b>
 * <p> The canonical {@link CuboctahedronSnubLevo}.
 * <p> It is the dual of the {@link IcositetrahedronPentagonalDextroCanonical}.
 * <p> It is the left-handed (<i>levo</i>) enantiomorph of the {@link CuboctahedronSnub}.
 * @see <a href="https://dmccooey.com/polyhedra/LsnubCube.html">
 *      Canonical <i>levo</i>-Snub Cuboctahedron (David McCooey)</a>
 */
public class CuboctahedronSnubLevoCanonical extends CuboctahedronSnubLevo {

    public CuboctahedronSnubLevoCanonical(
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
        vB0 = new Triad<>(C1(), C0(), C2());
        vB1 = new Triad<>(C1(), NEG_C0(), NEG_C2());
        vB2 = new Triad<>(NEG_C1(), NEG_C0(), C2());
        vB3 = new Triad<>(NEG_C1(), C0(), NEG_C2());
        vB4 = new Triad<>(C2(), C1(), C0());
        vB5 = new Triad<>(C2(), NEG_C1(), NEG_C0());
        vB6 = new Triad<>(NEG_C2(), NEG_C1(), C0());
        vB7 = new Triad<>(NEG_C2(), C1(), NEG_C0());
        vB8 = new Triad<>(C0(), C2(), C1());
        vB9 = new Triad<>(C0(), NEG_C2(), NEG_C1());
        vB10 = new Triad<>(NEG_C0(), NEG_C2(), C1());
        vB11 = new Triad<>(NEG_C0(), C2(), NEG_C1());
        vB12 = new Triad<>(C0(), NEG_C1(), C2());
        vB13 = new Triad<>(C0(), C1(), NEG_C2());
        vB14 = new Triad<>(NEG_C0(), C1(), C2());
        vB15 = new Triad<>(NEG_C0(), NEG_C1(), NEG_C2());
        vB16 = new Triad<>(C2(), NEG_C0(), C1());
        vB17 = new Triad<>(C2(), C0(), NEG_C1());
        vB18 = new Triad<>(NEG_C2(), C0(), C1());
        vB19 = new Triad<>(NEG_C2(), NEG_C0(), NEG_C1());
        vB20 = new Triad<>(C1(), NEG_C2(), C0());
        vB21 = new Triad<>(C1(), C2(), NEG_C0());
        vB22 = new Triad<>(NEG_C1(), C2(), C0());
        vB23 = new Triad<>(NEG_C1(), NEG_C2(), NEG_C0());

        setVerts();
    }
}
