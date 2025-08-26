package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.Polyad;
import io.github.noshou.tuple.Triad;
import org.jetbrains.annotations.NotNull;
import io.github.noshou.npg.shapes.archimedean.*;

/**
 * Represents a <b><i>levo<i/>Pentagonal Icositetrahedron</b>.
 * <p> The left-handed entantiomporph of the {@link IcositetrahedronPentagonal }.
 * <p> It is the dual of the {@link CuboctahedronSnubDextro}.</p>
 * @see <a href="https://dmccooey.com/polyhedra/RpentagonalIcositetrahedron.html">
 * <i>dextro</i>-Pentagonal Hexecontahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class IcositetrahedronPentagonalLevo extends IcositetrahedronPentagonal {
    public IcositetrahedronPentagonalLevo (
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
        vB0 = new Triad<>(N0(), N0(), NEG_C3());
        vB1 = new Triad<>(N0(), N0(), C3());
        vB2 = new Triad<>(NEG_C3(), N0(), N0());
        vB3  = new Triad<>(C3(), N0(), N0());
        vB4  = new Triad<>(N0(), NEG_C3(), N0());
        vB5  = new Triad<>(N0(), C3(), N0());
        vB6  = new Triad<>(NEG_C1(), C0(), NEG_C2());
        vB7  = new Triad<>(NEG_C1(), NEG_C0(), C2());
        vB8  = new Triad<>(C1(), NEG_C0(), NEG_C2());
        vB9  = new Triad<>(C1(), C0(), C2());
        vB10 = new Triad<>(NEG_C2(), C1(), NEG_C0());
        vB11 = new Triad<>(NEG_C2(), NEG_C1(), C0());
        vB12 = new Triad<>(C2(), NEG_C1(), NEG_C0());
        vB13 = new Triad<>(C2(), C1(), C0());
        vB14 = new Triad<>(NEG_C0(), C2(), NEG_C1());
        vB15 = new Triad<>(NEG_C0(), NEG_C2(), C1());
        vB16 = new Triad<>(C0(), NEG_C2(), NEG_C1());
        vB17 = new Triad<>(C0(), C2(), C1());
        vB18 = new Triad<>(NEG_C0(), NEG_C1(), NEG_C2());
        vB19 = new Triad<>(NEG_C0(), C1(), C2());
        vB20 = new Triad<>(C0(), C1(), NEG_C2());
        vB21 = new Triad<>(C0(), NEG_C1(), C2());
        vB22 = new Triad<>(NEG_C2(), NEG_C0(), NEG_C1());
        vB23 = new Triad<>(NEG_C2(), C0(), C1());
        vB24 = new Triad<>(C2(), C0(), NEG_C1());
        vB25 = new Triad<>(C2(), NEG_C0(), C1());
        vB26 = new Triad<>(NEG_C1(), NEG_C2(), NEG_C0());
        vB27 = new Triad<>(NEG_C1(), C2(), C0());
        vB28 = new Triad<>(C1(), C2(), NEG_C0());
        vB29 = new Triad<>(C1(), NEG_C2(), C0());
        vB30 = new Triad<>(NEG_C1(), NEG_C1(), NEG_C1());
        vB31 = new Triad<>(NEG_C1(), NEG_C1(), C1());
        vB32 = new Triad<>(NEG_C1(), C1(), NEG_C1());
        vB33 = new Triad<>(NEG_C1(), C1(), C1());
        vB34 = new Triad<>(C1(), NEG_C1(), NEG_C1());
        vB35 = new Triad<>(C1(), NEG_C1(), C1());
        vB36 = new Triad<>(C1(), C1(), NEG_C1());
        setVerts();
    }
}
