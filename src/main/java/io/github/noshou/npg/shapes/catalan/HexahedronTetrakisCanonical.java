package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.archimedean.*;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;

/**
 * Represents a <b>Canonical Tetrakis Hexahedron</b>.
 * <p> The canonical {@link HexahedronTetrakis}.
 * <p> It is the dual of the {@link OctahedronTruncatedCanonical}.
 * @see <a href="https://https://dmccooey.com/polyhedra/TetrakisHexahedron.html">
 *      Canonical Tetrakis Hexahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class HexahedronTetrakisCanonical extends HexahedronTetrakis {

    // Apfloat Constants
    private final Apfloat C0 = N3.multiply(SQRT2).divide(N4);
    private final Apfloat C1 = N9.multiply(SQRT2).divide(N8);

    public HexahedronTetrakisCanonical (
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

        setVertices(
                new Triad<>(N0, N0, C1),
                new Triad<>(C1, N0, N0),
                new Triad<>(C1.multiply(NEG_N1), N0, N0),
                new Triad<>(C1.multiply(NEG_N1), N0, N0),
                new Triad<>(N0, C1, N0),
                new Triad<>(N0, C1.multiply(NEG_N1), N0),
                new Triad<>(C0, C0, C0),
                new Triad<>(C0, C0, C0.multiply(NEG_N1)),
                new Triad<>(C0, C0.multiply(NEG_N1), C0),
                new Triad<>(C0, C0.multiply(NEG_N1), C0.multiply(NEG_N1)),
                new Triad<>(C0.multiply(NEG_N1), C0, C0),
                new Triad<>(C0.multiply(NEG_N1), C0, C0.multiply(NEG_N1)),
                new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), C0),
                new Triad<>(C0.multiply(NEG_N1), C0.multiply(NEG_N1), C0.multiply(NEG_N1))
        );
    }
}
