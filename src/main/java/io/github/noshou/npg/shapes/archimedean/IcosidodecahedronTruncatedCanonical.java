
 package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.catalan.TriacontahedronDisdyakisCanonical;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;


/**
 * Represents a <b> Canonical Truncated Icosidodecahedron</b>, an Archimedean solid.
 * <p>This polyhedron is formed by truncating (cutting off) the vertices of the
 * {@link Icosidodecahedron}. It consists of 62 faces: 30 squares, 20 hexagons,
 * and 12 decagons, with 180 edges and 120 vertices.
 * <p> It is the dual of the {@link TriacontahedronDisdyakisCanonical}.
 * @see <a href="https://dmccooey.com/polyhedra/TruncatedIcosidodecahedron.html">
 *      Truncated Icosidodecahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class IcosidodecahedronTruncatedCanonical extends IcosidodecahedronTruncated {

    // C0 = (3 + sqrt(5)) / 4
    private final Apfloat C0 = (N3.add(SQRT5)).divide(N4);

    // C1 = (1 + sqrt(5)) / 2
    private final Apfloat C1 = (N1.add(SQRT5)).divide(N2);

    // C2 = (5 + sqrt(5)) / 4
    private final Apfloat C2 = (N5.add(SQRT5)).divide(N4);

    // C3 = (2 + sqrt(5)) / 2
    private final Apfloat C3 = (N2.add(SQRT5)).divide(N2);

    // C4 =  3 * (1 + sqrt(5)) / 4
    private final Apfloat C4 = (N3.multiply(N1.add(SQRT5))).divide(N4);

    // C5 = (3 + sqrt(5)) / 2
    private final Apfloat C5 = (N3.add(SQRT5)).divide(N2);

    // C6 = (5 + 3 * sqrt(5)) / 4
    private final Apfloat C6 = (N5.add(N3.multiply(SQRT5))).divide(N4);

    // C7 =  (4 + sqrt(5)) / 2
    private final Apfloat C7 = (N4.add(SQRT5)).divide(N2);

    // C8 = (7 + 3 * sqrt(5)) / 4
    private final Apfloat C8 = (N7.add(N3.multiply(SQRT5))).divide(N4);

    // C9 = (3 + 2 * sqrt(5)) / 2
    private final Apfloat C9 = (N3.add(N2.multiply(SQRT5))).divide(N2);

    public IcosidodecahedronTruncatedCanonical(
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
        vB0 = new Triad<>(HALF, HALF, C9);
        vB1 = new Triad<>(HALF, HALF, C9.multiply(NEG_N1));
        vB2 = new Triad<>(HALF, HALF.multiply(NEG_N1), C9);
        vB3 = new Triad<>(HALF, HALF.multiply(NEG_N1), C9.multiply(NEG_N1));
        vB4 = new Triad<>(HALF.multiply(NEG_N1), HALF, C9);
        vB5 = new Triad<>(HALF.multiply(NEG_N1), HALF, C9.multiply(NEG_N1));
        vB6 = new Triad<>(HALF.multiply(NEG_N1), HALF.multiply(NEG_N1), C9);
        vB7 = new Triad<>(HALF.multiply(NEG_N1), HALF.multiply(NEG_N1), C9.multiply(NEG_N1));
        vB8 = new Triad<>(C9, HALF, HALF);
        vB9 = new Triad<>(C9, HALF, HALF.multiply(NEG_N1));
        vB10 = new Triad<>(C9, HALF.multiply(NEG_N1), HALF);
        vB11 = new Triad<>(C9, HALF.multiply(NEG_N1), HALF.multiply(NEG_N1));
        vB12 = new Triad<>(C9.multiply(NEG_N1), HALF, HALF);
        vB13 = new Triad<>(C9.multiply(NEG_N1), HALF, HALF.multiply(NEG_N1));
        vB14 = new Triad<>(C9.multiply(NEG_N1), HALF.multiply(NEG_N1), HALF);
        vB15 = new Triad<>(C9.multiply(NEG_N1), HALF.multiply(NEG_N1), HALF.multiply(NEG_N1));
        vB16 = new Triad<>(HALF, C9, HALF);
        vB17 = new Triad<>(HALF, C9, HALF.multiply(NEG_N1));
        vB18 = new Triad<>(HALF, C9.multiply(NEG_N1), HALF);
        vB19 = new Triad<>(HALF, C9.multiply(NEG_N1), HALF.multiply(NEG_N1));
        vB20 = new Triad<>(HALF.multiply(NEG_N1), C9, HALF);
        vB21 = new Triad<>(HALF.multiply(NEG_N1), C9, HALF.multiply(NEG_N1));
        vB22 = new Triad<>(HALF.multiply(NEG_N1), C9.multiply(NEG_N1), HALF);
        vB23 = new Triad<>(HALF.multiply(NEG_N1), C9.multiply(NEG_N1), HALF.multiply(NEG_N1));
        vB24 = new Triad<>(N1, C0, C8);
        vB25 = new Triad<>(N1, C0, C8.multiply(NEG_N1));
        vB26 = new Triad<>(N1, C0.multiply(NEG_N1), C8);
        vB27 = new Triad<>(N1, C0.multiply(NEG_N1), C8.multiply(NEG_N1));
        vB28 = new Triad<>(NEG_N1, C0, C8);
        vB29 = new Triad<>(NEG_N1, C0, C8.multiply(NEG_N1));
        vB30 = new Triad<>(NEG_N1, C0.multiply(NEG_N1), C8);
        vB31 = new Triad<>(NEG_N1, C0.multiply(NEG_N1), C8.multiply(NEG_N1));
        vB32 = new Triad<>(C8, N1, C0);
        vB33 = new Triad<>(C8, N1, C0.multiply(NEG_N1));
        vB34 = new Triad<>(C8, NEG_N1, C0);
        vB35 = new Triad<>(C8, NEG_N1, C0.multiply(NEG_N1));
        vB36 = new Triad<>(C8.multiply(NEG_N1), N1, C0);
        vB37 = new Triad<>(C8.multiply(NEG_N1), N1, C0.multiply(NEG_N1));
        vB38 = new Triad<>(C8.multiply(NEG_N1), NEG_N1, C0);
        vB39 = new Triad<>(C8.multiply(NEG_N1), NEG_N1, C0.multiply(NEG_N1));
        vB40 = new Triad<>(C0, C8, N1);
        vB41 = new Triad<>(C0, C8, NEG_N1);
        vB42 = new Triad<>(C0, C8.multiply(NEG_N1), N1);
        vB43 = new Triad<>(C0, C8.multiply(NEG_N1), NEG_N1);
        vB44 = new Triad<>(C0.multiply(NEG_N1), C8, N1);
        vB45 = new Triad<>(C0.multiply(NEG_N1), C8, NEG_N1);
        vB46 = new Triad<>(C0.multiply(NEG_N1), C8.multiply(NEG_N1), N1);
        vB47 = new Triad<>(C0.multiply(NEG_N1), C8.multiply(NEG_N1), NEG_N1);
        vB48 = new Triad<>(HALF, C3, C7);
        vB49 = new Triad<>(HALF, C3, C7.multiply(NEG_N1));
        vB50 = new Triad<>(HALF, C3.multiply(NEG_N1), C7);
        vB51 = new Triad<>(HALF, C3.multiply(NEG_N1), C7.multiply(NEG_N1));
        vB52 = new Triad<>(HALF.multiply(NEG_N1), C3, C7);
        vB53 = new Triad<>(HALF.multiply(NEG_N1), C3, C7.multiply(NEG_N1));
        vB54 = new Triad<>(HALF.multiply(NEG_N1), C3.multiply(NEG_N1), C7);
        vB55 = new Triad<>(HALF.multiply(NEG_N1), C3.multiply(NEG_N1), C7.multiply(NEG_N1));
        vB56 = new Triad<>(C7, HALF, C3);
        vB57 = new Triad<>(C7, HALF, C3.multiply(NEG_N1));
        vB58 = new Triad<>(C7, HALF.multiply(NEG_N1), C3);
        vB59 = new Triad<>(C7, HALF.multiply(NEG_N1), C3.multiply(NEG_N1));
        vB60 = new Triad<>(C7.multiply(NEG_N1), HALF, C3);
        vB61 = new Triad<>(C7.multiply(NEG_N1), HALF, C3.multiply(NEG_N1));
        vB62 = new Triad<>(C7.multiply(NEG_N1), HALF.multiply(NEG_N1), C3);
        vB63 = new Triad<>(C7.multiply(NEG_N1), HALF.multiply(NEG_N1), C3.multiply(NEG_N1));
        vB64 = new Triad<>(C3, C7, HALF);
        vB65 = new Triad<>(C3, C7, HALF.multiply(NEG_N1));
        vB66 = new Triad<>(C3, C7.multiply(NEG_N1), HALF);
        vB67 = new Triad<>(C3, C7.multiply(NEG_N1), HALF.multiply(NEG_N1));
        vB68 = new Triad<>(C3.multiply(NEG_N1), C7, HALF);
        vB69 = new Triad<>(C3.multiply(NEG_N1), C7, HALF.multiply(NEG_N1));
        vB70 = new Triad<>(C3.multiply(NEG_N1), C7.multiply(NEG_N1), HALF);
        vB71 = new Triad<>(C3.multiply(NEG_N1), C7.multiply(NEG_N1), HALF.multiply(NEG_N1));
        vB72 = new Triad<>(C2, C1, C6);
        vB73 = new Triad<>(C2, C1, C6.multiply(NEG_N1));
        vB74 = new Triad<>(C2, C1.multiply(NEG_N1), C6);
        vB75 = new Triad<>(C2, C1.multiply(NEG_N1), C6.multiply(NEG_N1));
        vB76 = new Triad<>(C2.multiply(NEG_N1), C1, C6);
        vB77 = new Triad<>(C2.multiply(NEG_N1), C1, C6.multiply(NEG_N1));
        vB78 = new Triad<>(C2.multiply(NEG_N1), C1.multiply(NEG_N1), C6);
        vB79 = new Triad<>(C2.multiply(NEG_N1), C1.multiply(NEG_N1), C6.multiply(NEG_N1));
        vB80 = new Triad<>(C6, C2, C1);
        vB81 = new Triad<>(C6, C2, C1.multiply(NEG_N1));
        vB82 = new Triad<>(C6, C2.multiply(NEG_N1), C1);
        vB83 = new Triad<>(C6, C2.multiply(NEG_N1), C1.multiply(NEG_N1));
        vB84 = new Triad<>(C6.multiply(NEG_N1), C2, C1);
        vB85 = new Triad<>(C6.multiply(NEG_N1), C2, C1.multiply(NEG_N1));
        vB86 = new Triad<>(C6.multiply(NEG_N1), C2.multiply(NEG_N1), C1);
        vB87 = new Triad<>(C6.multiply(NEG_N1), C2.multiply(NEG_N1), C1.multiply(NEG_N1));
        vB88 = new Triad<>(C1, C6, C2);
        vB89 = new Triad<>(C1, C6, C2.multiply(NEG_N1));
        vB90 = new Triad<>(C1, C6.multiply(NEG_N1), C2);
        vB91 = new Triad<>(C1, C6.multiply(NEG_N1), C2.multiply(NEG_N1));
        vB92 = new Triad<>(C1.multiply(NEG_N1), C6, C2);
        vB93 = new Triad<>(C1.multiply(NEG_N1), C6, C2.multiply(NEG_N1));
        vB94 = new Triad<>(C1.multiply(NEG_N1), C6.multiply(NEG_N1), C2);
        vB95 = new Triad<>(C1.multiply(NEG_N1), C6.multiply(NEG_N1), C2.multiply(NEG_N1));
        vB96 = new Triad<>(C0, C4, C5);
        vB97 = new Triad<>(C0, C4, C5.multiply(NEG_N1));
        vB98 = new Triad<>(C0, C4.multiply(NEG_N1), C5);
        vB99 = new Triad<>(C0, C4.multiply(NEG_N1), C5.multiply(NEG_N1));
        vB100 = new Triad<>(C0.multiply(NEG_N1), C4, C5);
        vB101 = new Triad<>(C0.multiply(NEG_N1), C4, C5.multiply(NEG_N1));
        vB102 = new Triad<>(C0.multiply(NEG_N1), C4.multiply(NEG_N1), C5);
        vB103 = new Triad<>(C0.multiply(NEG_N1), C4.multiply(NEG_N1), C5.multiply(NEG_N1));
        vB104 = new Triad<>(C5, C0, C4);
        vB105 = new Triad<>(C5, C0, C4.multiply(NEG_N1));
        vB106 = new Triad<>(C5, C0.multiply(NEG_N1), C4);
        vB107 = new Triad<>(C5, C0.multiply(NEG_N1), C4.multiply(NEG_N1));
        vB108 = new Triad<>(C5.multiply(NEG_N1), C0, C4);
        vB109 = new Triad<>(C5.multiply(NEG_N1), C0, C4.multiply(NEG_N1));
        vB110 = new Triad<>(C5.multiply(NEG_N1), C0.multiply(NEG_N1), C4);
        vB111 = new Triad<>(C5.multiply(NEG_N1), C0.multiply(NEG_N1), C4.multiply(NEG_N1));
        vB112 = new Triad<>(C4, C5, C0);
        vB113 = new Triad<>(C4, C5, C0.multiply(NEG_N1));
        vB114 = new Triad<>(C4, C5.multiply(NEG_N1), C0);
        vB115 = new Triad<>(C4, C5.multiply(NEG_N1), C0.multiply(NEG_N1));
        vB116 = new Triad<>(C4.multiply(NEG_N1), C5, C0);
        vB117 = new Triad<>(C4.multiply(NEG_N1), C5, C0.multiply(NEG_N1));
        vB118 = new Triad<>(C4.multiply(NEG_N1), C5.multiply(NEG_N1), C0);
        vB119 = new Triad<>(C4.multiply(NEG_N1), C5.multiply(NEG_N1), C0.multiply(NEG_N1));
    }
}