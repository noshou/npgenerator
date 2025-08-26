package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.archimedean.IcosidodecahedronTruncatedBiscribed;
import io.github.noshou.tuple.Polyad;
import io.github.noshou.tuple.Triad;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * Represents a <b>Biscribed {@link TriacontahedronDisdyakis}</b>.
 * <p>This Catalan solid is the dual  of the {@link IcosidodecahedronTruncatedBiscribed}.
 * It is characterized by 120 scalene triangular faces, 180 edges, and 62 vertices.
 * @see <a href="https://dmccooey.com/polyhedra/BiscribedDisdyakisTriacontahedron.html">
 *      Disdyakis Triacontahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class TriacontahedronDisdyakisBiscribed extends TriacontahedronDisdyakis {


    // C0 = (sqrt(5) - 1) / 4
    private final Apfloat C0 = (SQRT5.subtract(N1)).divide(N4);

    // C1 = (sqrt(15) - sqrt(3)) / 6
    private final Apfloat C1 = (SQRT15.subtract(SQRT3)).divide(N6);

    // C2 = sqrt(10 * (5 - sqrt(5))) / 10
    private final Apfloat C2 = ApfloatMath.sqrt(N10.multiply(N5.subtract(SQRT5))).divide(N10);

    // C3 = sqrt(3) / 3
    private final Apfloat C3 = SQRT3.divide(N3);

    // C4 =(1 + sqrt(5)) / 4
    private final Apfloat C4 = (N1.add(SQRT5)).divide(N4);

    // C5 = sqrt(10 * (5 + sqrt(5))) / 10
    private final Apfloat C5 = ApfloatMath.sqrt(N10.multiply(N5.add(SQRT5))).divide(N10);

    // C6 = (sqrt(3) + sqrt(15)) / 6
    private final Apfloat C6 = (SQRT3.add(SQRT15)).divide(N6);

    public TriacontahedronDisdyakisBiscribed(
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
        vB0  = new Triad<>(N0, N0, N1);
        vB1  = new Triad<>(N0, N0, NEG_N1);
        vB2  = new Triad<>(N1, N0, N0);
        vB3  = new Triad<>(NEG_N1, N0, N0);
        vB4  = new Triad<>(N0, N1, N0);
        vB5  = new Triad<>(N0, N1, N0);
        vB6  = new Triad<>(N0, C1, C6);
        vB7  = new Triad<>(N0, C1, C6.multiply(NEG_N1));
        vB8  = new Triad<>(N0, C1.multiply(NEG_N1), C6);
        vB9  = new Triad<>(N0, C1.multiply(NEG_N1), C6.multiply(NEG_N1));
        vB10 = new Triad<>(C6, N0, C1);
        vB11 = new Triad<>(C6, N0, C1.multiply(NEG_N1));
        vB12 = new Triad<>(C6.multiply(NEG_N1), N0, C1);
        vB13 = new Triad<>(C6.multiply(NEG_N1), N0, C1.multiply(NEG_N1));
        vB14 = new Triad<>(C1, C6, N0);
        vB15 = new Triad<>(C1, C6.multiply(NEG_N1), N0);
        vB16 = new Triad<>(C1.multiply(NEG_N1), C6, N0);
        vB17 = new Triad<>(C1.multiply(NEG_N1), C6.multiply(NEG_N1), N0);
        vB18 = new Triad<>(C2, N0, C5);
        vB19 = new Triad<>(C2, N0, C5.multiply(NEG_N1));
        vB20 = new Triad<>(C2.multiply(NEG_N1), N0, C5);
        vB21 = new Triad<>(C2.multiply(NEG_N1), N0, C5.multiply(NEG_N1));
        vB22 = new Triad<>(C5, C2, N0);
        vB23 = new Triad<>(C5, C2.multiply(NEG_N1), N0);
        vB24 = new Triad<>(C5.multiply(NEG_N1), C2, N0);
        vB25 = new Triad<>(C5.multiply(NEG_N1), C2.multiply(NEG_N1), N0);
        vB26 = new Triad<>(N0, C5, C2);
        vB27 = new Triad<>(N0, C5, C2.multiply(NEG_N1));
        vB28 = new Triad<>(N0, C5.multiply(NEG_N1), C2);
        vB29 = new Triad<>(N0, C5.multiply(NEG_N1), C2.multiply(NEG_N1));
        vB30 = new Triad<>(C0, HALF, C4);
        vB31 = new Triad<>(C0, HALF, C4.multiply(NEG_N1));
        vB32 = new Triad<>(C0, HALF.multiply(NEG_N1), C4);
        vB33 = new Triad<>(C0, HALF.multiply(NEG_N1), C4.multiply(NEG_N1));
        vB34 = new Triad<>(C0.multiply(NEG_N1), HALF, C4);
        vB35 = new Triad<>(C0.multiply(NEG_N1), HALF, C4.multiply(NEG_N1));
        vB36 = new Triad<>(C0.multiply(NEG_N1), HALF.multiply(NEG_N1), C4);
        vB37 = new Triad<>(C0.multiply(NEG_N1), HALF.multiply(NEG_N1), C4.multiply(NEG_N1));
        vB38 = new Triad<>(C4, C0, HALF);
        vB39 = new Triad<>(C4, C0, HALF.multiply(NEG_N1));
        vB40 = new Triad<>(C4, C0.multiply(NEG_N1), HALF);
        vB41 = new Triad<>(C4, C0.multiply(NEG_N1), HALF.multiply(NEG_N1));
        vB42 = new Triad<>(C4.multiply(NEG_N1), C0, HALF);
        vB43 = new Triad<>(C4.multiply(NEG_N1), C0, HALF.multiply(NEG_N1));
        vB44 = new Triad<>(C4.multiply(NEG_N1), C0.multiply(NEG_N1), HALF);
        vB45 = new Triad<>(C4.multiply(NEG_N1), C0.multiply(NEG_N1), HALF.multiply(NEG_N1));
        vB46 = new Triad<>(HALF, C4, C0);
        vB47 = new Triad<>(HALF, C4, C0.multiply(NEG_N1));
        vB48 = new Triad<>(HALF, C4.multiply(NEG_N1), C0);
        vB49 = new Triad<>(HALF, C4.multiply(NEG_N1), C0.multiply(NEG_N1));
        vB50 = new Triad<>(HALF.multiply(NEG_N1), C4, C0);
        vB51 = new Triad<>(HALF.multiply(NEG_N1), C4, C0.multiply(NEG_N1));
        vB52 = new Triad<>(HALF.multiply(NEG_N1), C4.multiply(NEG_N1), C0);
        vB53 = new Triad<>(HALF.multiply(NEG_N1), C4.multiply(NEG_N1), C0.multiply(NEG_N1));
        vB54 = new Triad<>(C3, C3, C3);
        vB55 = new Triad<>(C3, C3, C3.multiply(NEG_N1));
        vB56 = new Triad<>(C3, C3.multiply(NEG_N1), C3);
        vB57 = new Triad<>(C3, C3.multiply(NEG_N1), C3.multiply(NEG_N1));
        vB58 = new Triad<>(C3.multiply(NEG_N1), C3, C3);
        vB59 = new Triad<>(C3.multiply(NEG_N1), C3, C3.multiply(NEG_N1));
        vB60 = new Triad<>(C3.multiply(NEG_N1), C3.multiply(NEG_N1), C3);
        vB61 = new Triad<>(C3.multiply(NEG_N1), C3.multiply(NEG_N1), C3.multiply(NEG_N1));

        setVerts();
    }
}
