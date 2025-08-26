package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.archimedean.IcosidodecahedronTruncatedCanonical;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.*;

/**
 * Represents a <b>Canonical {@link TriacontahedronDisdyakis}</b>.
 * <p>This Catalan solid is the dual  of the {@link IcosidodecahedronTruncatedCanonical}.
 * It is characterized by 120 scalene triangular faces, 180 edges, and 62 vertices.
 * @see <a href="https://dmccooey.com/polyhedra/DisdyakisTriacontahedron.html">
 *      Disdyakis Triacontahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class TriacontahedronDisdyakisCanonical extends TriacontahedronDisdyakis {

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

    public TriacontahedronDisdyakisCanonical(
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
        vB0 = new Triad<>(N0, N0, C8);
        vB1 = new Triad<>(N0, N0, C8.multiply(NEG_N1));
        vB2  = new Triad<>(C8, N0, N0);
        vB3  = new Triad<>(C8.multiply(NEG_N1), N0, N0);
        vB4  = new Triad<>(N0, C8, N0);
        vB5  = new Triad<>(N0, C8.multiply(NEG_N1), N0);
        vB6  = new Triad<>(N0, C1, C7);
        vB7  = new Triad<>(N0, C1, C7.multiply(NEG_N1));
        vB8  = new Triad<>(N0, C1.multiply(NEG_N1), C7);
        vB9  = new Triad<>(N0, C1.multiply(NEG_N1), C7.multiply(NEG_N1));
        vB10 = new Triad<>(C7, N0, C1);
        vB11 = new Triad<>(C7, N0, C1.multiply(NEG_N1));
        vB12 = new Triad<>(C7.multiply(NEG_N1), N0, C1);
        vB13 = new Triad<>(C7.multiply(NEG_N1), N0, C1.multiply(NEG_N1));
        vB14 = new Triad<>(C1, C7, N0);
        vB15 = new Triad<>(C1, C7.multiply(NEG_N1), N0);
        vB16 = new Triad<>(C1.multiply(NEG_N1), C7, N0);
        vB17 = new Triad<>(C1.multiply(NEG_N1), C7.multiply(NEG_N1), N0);
        vB18 = new Triad<>(C3, N0, C6);
        vB19 = new Triad<>(C3, N0, C6.multiply(NEG_N1));
        vB20 = new Triad<>(C3.multiply(NEG_N1), N0, C6);
        vB21 = new Triad<>(C3.multiply(NEG_N1), N0, C6.multiply(NEG_N1));
        vB22 = new Triad<>(C6, C3, N0);
        vB23 = new Triad<>(C6, C3.multiply(NEG_N1), N0);
        vB24 = new Triad<>(C6.multiply(NEG_N1), C3, N0);
        vB25 = new Triad<>(C6.multiply(NEG_N1), C3.multiply(NEG_N1), N0);
        vB26 = new Triad<>(N0, C6, C3);
        vB27 = new Triad<>(N0, C6, C3.multiply(NEG_N1));
        vB28 = new Triad<>(N0, C6.multiply(NEG_N1), C3);
        vB29 = new Triad<>(N0, C6.multiply(NEG_N1), C3.multiply(NEG_N1));
        vB30 = new Triad<>(C0, C2, C5);
        vB31 = new Triad<>(C0, C2, C5.multiply(NEG_N1));
        vB32 = new Triad<>(C0, C2.multiply(NEG_N1), C5);
        vB33 = new Triad<>(C0, C2.multiply(NEG_N1), C5.multiply(NEG_N1));
        vB34 = new Triad<>(C0.multiply(NEG_N1), C2, C5);
        vB35 = new Triad<>(C0.multiply(NEG_N1), C2, C5.multiply(NEG_N1));
        vB36 = new Triad<>(C0.multiply(NEG_N1), C2.multiply(NEG_N1), C5);
        vB37 = new Triad<>(C0.multiply(NEG_N1), C2.multiply(NEG_N1), C5.multiply(NEG_N1));
        vB38 = new Triad<>(C5, C0, C2);
        vB39 = new Triad<>(C5, C0, C2.multiply(NEG_N1));
        vB40 = new Triad<>(C5, C0.multiply(NEG_N1), C2);
        vB41 = new Triad<>(C5, C0.multiply(NEG_N1), C2.multiply(NEG_N1));
        vB42 = new Triad<>(C5.multiply(NEG_N1), C0, C2);
        vB43 = new Triad<>(C5.multiply(NEG_N1), C0, C2.multiply(NEG_N1));
        vB44 = new Triad<>(C5.multiply(NEG_N1), C0.multiply(NEG_N1), C2);
        vB45 = new Triad<>(C5.multiply(NEG_N1), C0.multiply(NEG_N1), C2.multiply(NEG_N1));
        vB46 = new Triad<>(C2, C5, C0);
        vB47 = new Triad<>(C2, C5, C0.multiply(NEG_N1));
        vB48 = new Triad<>(C2, C5.multiply(NEG_N1), C0);
        vB49 = new Triad<>(C2, C5.multiply(NEG_N1), C0.multiply(NEG_N1));
        vB50 = new Triad<>(C2.multiply(NEG_N1), C5, C0);
        vB51 = new Triad<>(C2.multiply(NEG_N1), C5, C0.multiply(NEG_N1));
        vB52 = new Triad<>(C2.multiply(NEG_N1), C5.multiply(NEG_N1), C0);
        vB53 = new Triad<>(C2.multiply(NEG_N1), C5.multiply(NEG_N1), C0.multiply(NEG_N1));
        vB54 = new Triad<>(C4, C4, C4);
        vB55 = new Triad<>(C4, C4, C4.multiply(NEG_N1));
        vB56 = new Triad<>(C4, C4.multiply(NEG_N1), C4);
        vB57 = new Triad<>(C4, C4.multiply(NEG_N1), C4.multiply(NEG_N1));
        vB58 = new Triad<>(C4.multiply(NEG_N1), C4, C4);
        vB59 = new Triad<>(C4.multiply(NEG_N1), C4, C4.multiply(NEG_N1));
        vB60 = new Triad<>(C4.multiply(NEG_N1), C4.multiply(NEG_N1), C4);
        vB61 = new Triad<>(C4.multiply(NEG_N1), C4.multiply(NEG_N1), C4.multiply(NEG_N1));

        setVerts();

    }

}
