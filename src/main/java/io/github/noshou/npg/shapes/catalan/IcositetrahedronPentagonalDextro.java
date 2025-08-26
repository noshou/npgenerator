package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;
import io.github.noshou.npg.shapes.archimedean.*;

/**
 * Represents a <b><i>dextro</i>-Snub Cuboctahedron</b>
 * <p> The right-handed (<i>levo</i>) enantiomorph of the {@link IcositetrahedronPentagonal}.
 * <p> It is the dual of the {@link CuboctahedronSnubLevo}. </p>
 * <p> May be circumscribed ({@link IcositetrahedronPentagonalDextroCanonical}) or
 * circumscribed and inscribed ({@link IcositetrahedronPentagonalDextroBiscribed}).
 */
@SuppressWarnings("FieldCanBeLocal")
public abstract class IcositetrahedronPentagonalDextro extends IcositetrahedronPentagonal {

    public IcositetrahedronPentagonalDextro(
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
    }
}
