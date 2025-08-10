
package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.Polyad;
import io.github.noshou.tuple.Triad;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;

/**
* Represents a <b> Biscribed Truncated Icosidodecahedron</b>, an Archimedean solid.
* <p>This polyhedron is formed by biscribing and truncating (cutting off) the vertices of the
* {@link Icosidodecahedron}. It consists of 62 faces: 30 squares, 20 hexagons,
* and 12 decagons, with 180 edges and 120 vertices.
* <p> It is the dual of the {@link IcosahedronHexasisBiscribed}.
* @see <a href="https://dmccooey.com/polyhedra/BiscribedTruncatedIcosidodecahedron.html">
*      Truncated Icosidodecahedron (David McCooey)</a>
*/
@SuppressWarnings("FieldCanBeLocal")
public class IcosidodecahedronTruncatedBiscribed extends IcosidodecahedronTruncated {


    // C0  = (4 - 3 * sqrt(3) + 4 * sqrt(5) - sqrt(15) - sqrt(2*(5 + sqrt(5)))) / 4
    private final Apfloat C0 = N4
            .subtract(N3.multiply(SQRT3))
            .add(N4.multiply(SQRT5))
            .subtract(SQRT15)
            .subtract(SQRT_2_TIMES_5_PLUS_SQRT5)
            .divide(N4);

    // C1  = (sqrt(3) - 3 - sqrt(5) + sqrt(15)) / 2
    private final Apfloat C1 = SQRT3
            .subtract(N3)
            .subtract(SQRT5)
            .add(SQRT15)
            .divide(N2);

    // C2  = (3 - sqrt(3) + sqrt(5) - sqrt(5 + 2 * sqrt(5))) / 2
    private final Apfloat C2 = N3
            .subtract(SQRT3)
            .add(SQRT5)
            .subtract(SQRT_5_PLUS_2_SQRT5)
            .divide(N2);

    // C3  = (sqrt(2 * (5 + sqrt(5))) - 1 - sqrt(5)) / 2
    private final Apfloat C3 = SQRT_2_TIMES_5_PLUS_SQRT5
            .subtract(N1)
            .subtract(SQRT5)
            .divide(N2);

    // C4  = (2 - 3 * sqrt(3) + 2 * sqrt(5) - sqrt(15) + sqrt(2*(5 + sqrt(5)))) / 4
    private final Apfloat C4 = N2
            .subtract(N3.multiply(SQRT3))
            .add(N2.multiply(SQRT5))
            .subtract(SQRT15)
            .add(SQRT_2_TIMES_5_PLUS_SQRT5)
            .divide(N4);

    // C5  = (3 * sqrt(3) - 4 + sqrt(15) - sqrt(2 * (5 + sqrt(5)))) / 4
    private final Apfloat C5 = N3.multiply(SQRT3)
            .subtract(N4)
            .add(SQRT15)
            .subtract(SQRT_2_TIMES_5_PLUS_SQRT5)
            .divide(N4);

    // C6  = (2 + sqrt(3) - sqrt(5 + 2 * sqrt(5))) / 2
    private final Apfloat C6 = N2
            .add(SQRT3)
            .subtract(SQRT_5_PLUS_2_SQRT5)
            .divide(N2);

    // C7  = (3 * sqrt(3) - 6 - 2 * sqrt(5) + sqrt(15) + sqrt(2*(5 + sqrt(5)))) / 4
    private final Apfloat C7 = N3.multiply(SQRT3)
            .subtract(N6)
            .subtract(N2.multiply(SQRT5))
            .add(SQRT15)
            .add(SQRT_2_TIMES_5_PLUS_SQRT5)
            .divide(N4);

    // C8  = (6 - sqrt(3) + 2 * sqrt(5) - sqrt(15) - sqrt(2 * (5 - sqrt(5)))) / 4
    private final Apfloat C8 = N6
            .subtract(SQRT3)
            .add(N2.multiply(SQRT5))
            .subtract(SQRT15)
            .subtract(SQRT_2_TIMES_5_MINUS_SQRT5)
            .divide(N4);

    // C9  = (sqrt(5 + 2 * sqrt(5)) - sqrt(3)) / 2
    private final Apfloat C9 = SQRT_5_PLUS_2_SQRT5
            .subtract(SQRT3)
            .divide(N2);

    // C10 = (sqrt(3) - 1 - sqrt(5) + sqrt(5 + 2 * sqrt(5))) / 2
    private final Apfloat C10 = SQRT3
            .subtract(N1)
            .subtract(SQRT5)
            .add(SQRT_5_PLUS_2_SQRT5)
            .divide(N2);

    // C11 = (2 - sqrt(3) + 2 * sqrt(5) - sqrt(15) + sqrt(2 * (5 - sqrt(5)))) / 4
    private final Apfloat C11 = N2
            .subtract(SQRT3)
            .add(N2.multiply(SQRT5))
            .subtract(SQRT15)
            .add(SQRT_2_TIMES_5_MINUS_SQRT5)
            .divide(N4);

    // C12 = (sqrt(3) + sqrt(15) - sqrt(2 * (5 - sqrt(5)))) / 4
    private final Apfloat C12 = SQRT3
            .add(SQRT15)
            .subtract(SQRT_2_TIMES_5_MINUS_SQRT5)
            .divide(N4);

    // C13 = (sqrt(3) - 4 + sqrt(15) + sqrt(2 * (5 - sqrt(5)))) / 4
    private final Apfloat C13 = SQRT3
            .subtract(N4)
            .add(SQRT15)
            .add(SQRT_2_TIMES_5_MINUS_SQRT5)
            .divide(N4);

    public IcosidodecahedronTruncatedBiscribed(
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
        vB0 = new Triad<>(C3, C1, N1);
        vB1 = new Triad<>(C3, C1, NEG_N1);
        vB2 = new Triad<>(C3, C1.multiply(NEG_N1), N1);
        vB3  = new Triad<>(C3, C1.multiply(NEG_N1), NEG_N1);
        vB4  = new Triad<>(C3.multiply(NEG_N1), C1, N1);
        vB5  = new Triad<>(C3.multiply(NEG_N1), C1, NEG_N1);
        vB6  = new Triad<>(C3.multiply(NEG_N1), C1.multiply(NEG_N1), N1);
        vB7  = new Triad<>(C3.multiply(NEG_N1), C1.multiply(NEG_N1), NEG_N1);
        vB8  = new Triad<>(N1, C3, C1);
        vB9  = new Triad<>(N1, C3, C1.multiply(NEG_N1));
        vB10 = new Triad<>(N1, C3.multiply(NEG_N1), C1);
        vB11 = new Triad<>(N1, C3.multiply(NEG_N1), C1.multiply(NEG_N1));
        vB12 = new Triad<>(NEG_N1, C3, C1);
        vB13 = new Triad<>(NEG_N1, C3, C1.multiply(NEG_N1));
        vB14 = new Triad<>(NEG_N1, C3.multiply(NEG_N1), C1);
        vB15 = new Triad<>(NEG_N1, C3.multiply(NEG_N1), C1.multiply(NEG_N1));
        vB16 = new Triad<>(C1, N1, C3);
        vB17 = new Triad<>(C1, N1, C3.multiply(NEG_N1));
        vB18 = new Triad<>(C1, NEG_N1, C3);
        vB19 = new Triad<>(C1, NEG_N1, C3.multiply(NEG_N1));
        vB20 = new Triad<>(C1.multiply(NEG_N1), N1, C3);
        vB21 = new Triad<>(C1.multiply(NEG_N1), N1, C3.multiply(NEG_N1));
        vB22 = new Triad<>(C1.multiply(NEG_N1), NEG_N1, C3);
        vB23 = new Triad<>(C1.multiply(NEG_N1), NEG_N1, C3.multiply(NEG_N1));
        vB24 = new Triad<>(C4, C2, C13);
        vB25 = new Triad<>(C4, C2, C13.multiply(NEG_N1));
        vB26 = new Triad<>(C4, C2.multiply(NEG_N1), C13);
        vB27 = new Triad<>(C4, C2.multiply(NEG_N1), C13.multiply(NEG_N1));
        vB28 = new Triad<>(C4.multiply(NEG_N1), C2, C13);
        vB29 = new Triad<>(C4.multiply(NEG_N1), C2, C13.multiply(NEG_N1));
        vB30 = new Triad<>(C4.multiply(NEG_N1), C2.multiply(NEG_N1), C13);
        vB31 = new Triad<>(C4.multiply(NEG_N1), C2.multiply(NEG_N1), C13.multiply(NEG_N1));
        vB32 = new Triad<>(C13, C4, C2);
        vB33 = new Triad<>(C13, C4, C2.multiply(NEG_N1));
        vB34 = new Triad<>(C13, C4.multiply(NEG_N1), C2);
        vB35 = new Triad<>(C13, C4.multiply(NEG_N1), C2.multiply(NEG_N1));
        vB36 = new Triad<>(C13.multiply(NEG_N1), C4, C2);
        vB37 = new Triad<>(C13.multiply(NEG_N1), C4, C2.multiply(NEG_N1));
        vB38 = new Triad<>(C13.multiply(NEG_N1), C4.multiply(NEG_N1), C2);
        vB39 = new Triad<>(C13.multiply(NEG_N1), C4.multiply(NEG_N1), C2.multiply(NEG_N1));
        vB40 = new Triad<>(C2, C13, C4);
        vB41 = new Triad<>(C2, C13, C4.multiply(NEG_N1));
        vB42 = new Triad<>(C2, C13.multiply(NEG_N1), C4);
        vB43 = new Triad<>(C2, C13.multiply(NEG_N1), C4.multiply(NEG_N1));
        vB44 = new Triad<>(C2.multiply(NEG_N1), C13, C4);
        vB45 = new Triad<>(C2.multiply(NEG_N1), C13, C4.multiply(NEG_N1));
        vB46 = new Triad<>(C2.multiply(NEG_N1), C13.multiply(NEG_N1), C4);
        vB47 = new Triad<>(C2.multiply(NEG_N1), C13.multiply(NEG_N1), C4.multiply(NEG_N1));
        vB48 = new Triad<>(C0, C9, C12);
        vB49 = new Triad<>(C0, C9, C12.multiply(NEG_N1));
        vB50 = new Triad<>(C0, C9.multiply(NEG_N1), C12);
        vB51 = new Triad<>(C0, C9.multiply(NEG_N1), C12.multiply(NEG_N1));
        vB52 = new Triad<>(C0.multiply(NEG_N1), C9, C12);
        vB53 = new Triad<>(C0.multiply(NEG_N1), C9, C12.multiply(NEG_N1));
        vB54 = new Triad<>(C0.multiply(NEG_N1), C9.multiply(NEG_N1), C12);
        vB55 = new Triad<>(C0.multiply(NEG_N1), C9.multiply(NEG_N1), C12.multiply(NEG_N1));
        vB56 = new Triad<>(C12, C0, C9);
        vB57 = new Triad<>(C12, C0, C9.multiply(NEG_N1));
        vB58 = new Triad<>(C12, C0.multiply(NEG_N1), C9);
        vB59 = new Triad<>(C12, C0.multiply(NEG_N1), C9.multiply(NEG_N1));
        vB60 = new Triad<>(C12.multiply(NEG_N1), C0, C9);
        vB61 = new Triad<>(C12.multiply(NEG_N1), C0, C9.multiply(NEG_N1));
        vB62 = new Triad<>(C12.multiply(NEG_N1), C0.multiply(NEG_N1), C9);
        vB63 = new Triad<>(C12.multiply(NEG_N1), C0.multiply(NEG_N1), C9.multiply(NEG_N1));
        vB64 = new Triad<>(C9, C12, C0);
        vB65 = new Triad<>(C9, C12, C0.multiply(NEG_N1));
        vB66 = new Triad<>(C9, C12.multiply(NEG_N1), C0);
        vB67 = new Triad<>(C9, C12.multiply(NEG_N1), C0.multiply(NEG_N1));
        vB68 = new Triad<>(C9.multiply(NEG_N1), C12, C0);
        vB69 = new Triad<>(C9.multiply(NEG_N1), C12, C0.multiply(NEG_N1));
        vB70 = new Triad<>(C9.multiply(NEG_N1), C12.multiply(NEG_N1), C0);
        vB71 = new Triad<>(C9.multiply(NEG_N1), C12.multiply(NEG_N1), C0.multiply(NEG_N1));
        vB72 = new Triad<>(C7, C6, C11);
        vB73 = new Triad<>(C7, C6, C11.multiply(NEG_N1));
        vB74 = new Triad<>(C7, C6.multiply(NEG_N1), C11);
        vB75 = new Triad<>(C7, C6.multiply(NEG_N1), C11.multiply(NEG_N1));
        vB76 = new Triad<>(C7.multiply(NEG_N1), C6, C11);
        vB77 = new Triad<>(C7.multiply(NEG_N1), C6, C11.multiply(NEG_N1));
        vB78 = new Triad<>(C7.multiply(NEG_N1), C6.multiply(NEG_N1), C11);
        vB79 = new Triad<>(C7.multiply(NEG_N1), C6.multiply(NEG_N1), C11.multiply(NEG_N1));
        vB80 = new Triad<>(C11, C7, C6);
        vB81 = new Triad<>(C11, C7, C6.multiply(NEG_N1));
        vB82 = new Triad<>(C11, C7.multiply(NEG_N1), C6);
        vB83 = new Triad<>(C11, C7.multiply(NEG_N1), C6.multiply(NEG_N1));
        vB84 = new Triad<>(C11.multiply(NEG_N1), C7, C6);
        vB85 = new Triad<>(C11.multiply(NEG_N1), C7, C6.multiply(NEG_N1));
        vB86 = new Triad<>(C11.multiply(NEG_N1), C7.multiply(NEG_N1), C6);
        vB87 = new Triad<>(C11.multiply(NEG_N1), C7.multiply(NEG_N1), C6.multiply(NEG_N1));
        vB88 = new Triad<>(C6, C11, C7);
        vB89 = new Triad<>(C6, C11, C7.multiply(NEG_N1));
        vB90 = new Triad<>(C6, C11.multiply(NEG_N1), C7);
        vB91 = new Triad<>(C6, C11.multiply(NEG_N1), C7.multiply(NEG_N1));
        vB92 = new Triad<>(C6.multiply(NEG_N1), C11, C7);
        vB93 = new Triad<>(C6.multiply(NEG_N1), C11, C7.multiply(NEG_N1));
        vB94 = new Triad<>(C6.multiply(NEG_N1), C11.multiply(NEG_N1), C7);
        vB95 = new Triad<>(C6.multiply(NEG_N1), C11.multiply(NEG_N1), C7.multiply(NEG_N1));
        vB96 = new Triad<>(C5, C10, C8);
        vB97 = new Triad<>(C5, C10, C8.multiply(NEG_N1));
        vB98 = new Triad<>(C5, C10.multiply(NEG_N1), C8);
        vB99 = new Triad<>(C5, C10.multiply(NEG_N1), C8.multiply(NEG_N1));
        vB100 = new Triad<>(C5.multiply(NEG_N1), C10, C8);
        vB101 = new Triad<>(C5.multiply(NEG_N1), C10, C8.multiply(NEG_N1));
        vB102 = new Triad<>(C5.multiply(NEG_N1), C10.multiply(NEG_N1), C8);
        vB103 = new Triad<>(C5.multiply(NEG_N1), C10.multiply(NEG_N1), C8.multiply(NEG_N1));
        vB104 = new Triad<>(C8, C5, C10);
        vB105 = new Triad<>(C8, C5, C10.multiply(NEG_N1));
        vB106 = new Triad<>(C8, C5.multiply(NEG_N1), C10);
        vB107 = new Triad<>(C8, C5.multiply(NEG_N1), C10.multiply(NEG_N1));
        vB108 = new Triad<>(C8.multiply(NEG_N1), C5, C10);
        vB109 = new Triad<>(C8.multiply(NEG_N1), C5, C10.multiply(NEG_N1));
        vB110 = new Triad<>(C8.multiply(NEG_N1), C5.multiply(NEG_N1), C10);
        vB111 = new Triad<>(C8.multiply(NEG_N1), C5.multiply(NEG_N1), C10.multiply(NEG_N1));
        vB112 = new Triad<>(C10, C8, C5);
        vB113 = new Triad<>(C10, C8, C5.multiply(NEG_N1));
        vB114 = new Triad<>(C10, C8.multiply(NEG_N1), C5);
        vB115 = new Triad<>(C10, C8.multiply(NEG_N1), C5.multiply(NEG_N1));
        vB116 = new Triad<>(C10.multiply(NEG_N1), C8, C5);
        vB117 = new Triad<>(C10.multiply(NEG_N1), C8, C5.multiply(NEG_N1));
        vB118 = new Triad<>(C10.multiply(NEG_N1), C8.multiply(NEG_N1), C5);
        vB119 = new Triad<>(C10.multiply(NEG_N1), C8.multiply(NEG_N1), C5.multiply(NEG_N1));
    }
}