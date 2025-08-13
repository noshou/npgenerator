package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.catalan.*;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;


/**
 * Represents a <b>Canonical Truncated Octahedron</b>.
 * <p> The canonical circumscribed {@link OctahedronTruncated}.
 * <p> It is the dual of the {@link HexahedronTetrakisCanonical}.
 * @see <a href="https://dmccooey.com/polyhedra/TruncatedOctahedron.html">
 *      Canonical Truncated Octahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class OctahedronTruncatedCanonical extends OctahedronTruncated {

    // Apfloat Constants
    private final Apfloat C0 = SQRT2.divide(N2);
    private final Apfloat C1 = SQRT2;

    public OctahedronTruncatedCanonical(
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

        // ==== CANONICAL VERTICES ====
        setVertices(
                new Triad<>(C0, N0, C1),
                new Triad<>(C0, N0, C1.multiply(NEG_N1)),
                new Triad<>(C0.multiply(NEG_N1), N0, C1),
                new Triad<>(C0.multiply(NEG_N1), N0, C1.multiply(NEG_N1)),
                new Triad<>(C1, C0, N0),
                new Triad<>(C1, C0.multiply(NEG_N1), N0),
                new Triad<>(C1.multiply(NEG_N1), C0, N0),
                new Triad<>(C1.multiply(NEG_N1), C0.multiply(NEG_N1), N0),
                new Triad<>(N0, C1, C0),
                new Triad<>(N0, C1, C0.multiply(NEG_N1)),
                new Triad<>(N0, C1.multiply(NEG_N1), C0),
                new Triad<>(N0, C1.multiply(NEG_N1), C0.multiply(NEG_N1)),
                new Triad<>(N0, C0, C1),
                new Triad<>(N0, C0, C1.multiply(NEG_N1)),
                new Triad<>(N0, C0.multiply(NEG_N1), C1),
                new Triad<>(N0, C0.multiply(NEG_N1), C1.multiply(NEG_N1)),
                new Triad<>(C1, N0, C0),
                new Triad<>(C1, N0, C0.multiply(NEG_N1)),
                new Triad<>(C1.multiply(NEG_N1), N0, C0),
                new Triad<>(C1.multiply(NEG_N1), N0, C0.multiply(NEG_N1)),
                new Triad<>(C0, C1, N0),
                new Triad<>(C0, C1.multiply(NEG_N1), N0),
                new Triad<>(C0.multiply(NEG_N1), C1, N0),
                new Triad<>(C0.multiply(NEG_N1), C1.multiply(NEG_N1), N0)
        );
    }
}
