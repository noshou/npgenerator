package io.github.noshou.npg.shapes.archimedean;
import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.npg.shapes.catalan.TriacontahedronDisdyakis;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b>Truncated Icosidodecahedron</b>.
 * <p>This polyhedron is formed by truncating (cutting off) the vertices of an
 *  {@link Icosidodecahedron}, resulting in a convex solid with 120 vertices, 180 edges,
 * and 62 faces composed of regular polygons: 12 decagons, 30 squares, and 20 hexagons.
 * The truncated icosidodecahedron can appear in two common forms:
 * <ul>
 *   <li>{@link IcosidodecahedronTruncatedCanonical} is the standard, symmetric shape with equal edge lengths.</li>
 *   <li>{@link IcosidodecahedronTruncatedBiscribed} is inscribed between
 *       two concentric spheres—one touching all vertices and the other touching
 *       all face centers—resulting in slightly different vertex positions and
 *       subtle shape variations.</li>
 * </ul>
 * <p> It is the dual of the {@link TriacontahedronDisdyakis}.
 */
public abstract class IcosidodecahedronTruncated extends Shape {

    // basis vertices
    protected Triad<Apfloat> vB0;
    protected Triad<Apfloat> vB1;
    protected Triad<Apfloat> vB2;
    protected Triad<Apfloat> vB3;
    protected Triad<Apfloat> vB4;
    protected Triad<Apfloat> vB5;
    protected Triad<Apfloat> vB6;
    protected Triad<Apfloat> vB7;
    protected Triad<Apfloat> vB8;
    protected Triad<Apfloat> vB9;
    protected Triad<Apfloat> vB10;
    protected Triad<Apfloat> vB11;
    protected Triad<Apfloat> vB12;
    protected Triad<Apfloat> vB13;
    protected Triad<Apfloat> vB14;
    protected Triad<Apfloat> vB15;
    protected Triad<Apfloat> vB16;
    protected Triad<Apfloat> vB17;
    protected Triad<Apfloat> vB18;
    protected Triad<Apfloat> vB19;
    protected Triad<Apfloat> vB20;
    protected Triad<Apfloat> vB21;
    protected Triad<Apfloat> vB22;
    protected Triad<Apfloat> vB23;
    protected Triad<Apfloat> vB24;
    protected Triad<Apfloat> vB25;
    protected Triad<Apfloat> vB26;
    protected Triad<Apfloat> vB27;
    protected Triad<Apfloat> vB28;
    protected Triad<Apfloat> vB29;
    protected Triad<Apfloat> vB30;
    protected Triad<Apfloat> vB31;
    protected Triad<Apfloat> vB32;
    protected Triad<Apfloat> vB33;
    protected Triad<Apfloat> vB34;
    protected Triad<Apfloat> vB35;
    protected Triad<Apfloat> vB36;
    protected Triad<Apfloat> vB37;
    protected Triad<Apfloat> vB38;
    protected Triad<Apfloat> vB39;
    protected Triad<Apfloat> vB40;
    protected Triad<Apfloat> vB41;
    protected Triad<Apfloat> vB42;
    protected Triad<Apfloat> vB43;
    protected Triad<Apfloat> vB44;
    protected Triad<Apfloat> vB45;
    protected Triad<Apfloat> vB46;
    protected Triad<Apfloat> vB47;
    protected Triad<Apfloat> vB48;
    protected Triad<Apfloat> vB49;
    protected Triad<Apfloat> vB50;
    protected Triad<Apfloat> vB51;
    protected Triad<Apfloat> vB52;
    protected Triad<Apfloat> vB53;
    protected Triad<Apfloat> vB54;
    protected Triad<Apfloat> vB55;
    protected Triad<Apfloat> vB56;
    protected Triad<Apfloat> vB57;
    protected Triad<Apfloat> vB58;
    protected Triad<Apfloat> vB59;
    protected Triad<Apfloat> vB60;
    protected Triad<Apfloat> vB61;
    protected Triad<Apfloat> vB62;
    protected Triad<Apfloat> vB63;
    protected Triad<Apfloat> vB64;
    protected Triad<Apfloat> vB65;
    protected Triad<Apfloat> vB66;
    protected Triad<Apfloat> vB67;
    protected Triad<Apfloat> vB68;
    protected Triad<Apfloat> vB69;
    protected Triad<Apfloat> vB70;
    protected Triad<Apfloat> vB71;
    protected Triad<Apfloat> vB72;
    protected Triad<Apfloat> vB73;
    protected Triad<Apfloat> vB74;
    protected Triad<Apfloat> vB75;
    protected Triad<Apfloat> vB76;
    protected Triad<Apfloat> vB77;
    protected Triad<Apfloat> vB78;
    protected Triad<Apfloat> vB79;
    protected Triad<Apfloat> vB80;
    protected Triad<Apfloat> vB81;
    protected Triad<Apfloat> vB82;
    protected Triad<Apfloat> vB83;
    protected Triad<Apfloat> vB84;
    protected Triad<Apfloat> vB85;
    protected Triad<Apfloat> vB86;
    protected Triad<Apfloat> vB87;
    protected Triad<Apfloat> vB88;
    protected Triad<Apfloat> vB89;
    protected Triad<Apfloat> vB90;
    protected Triad<Apfloat> vB91;
    protected Triad<Apfloat> vB92;
    protected Triad<Apfloat> vB93;
    protected Triad<Apfloat> vB94;
    protected Triad<Apfloat> vB95;
    protected Triad<Apfloat> vB96;
    protected Triad<Apfloat> vB97;
    protected Triad<Apfloat> vB98;
    protected Triad<Apfloat> vB99;
    protected Triad<Apfloat> vB100;
    protected Triad<Apfloat> vB101;
    protected Triad<Apfloat> vB102;
    protected Triad<Apfloat> vB103;
    protected Triad<Apfloat> vB104;
    protected Triad<Apfloat> vB105;
    protected Triad<Apfloat> vB106;
    protected Triad<Apfloat> vB107;
    protected Triad<Apfloat> vB108;
    protected Triad<Apfloat> vB109;
    protected Triad<Apfloat> vB110;
    protected Triad<Apfloat> vB111;
    protected Triad<Apfloat> vB112;
    protected Triad<Apfloat> vB113;
    protected Triad<Apfloat> vB114;
    protected Triad<Apfloat> vB115;
    protected Triad<Apfloat> vB116;
    protected Triad<Apfloat> vB117;
    protected Triad<Apfloat> vB118;
    protected Triad<Apfloat> vB119;

    // 30 square faces
    private final ArrayList<Tetrad<Tuple<Apfloat>>> faces_sqr = new ArrayList<>();
    private final ArrayList<Triad<Apfloat>> face_norms_sqr = new ArrayList<>();

    // 20 hexagonal faces
    private final Icosad<Tuple<Tuple<Apfloat>>> faces_hex;
    private final Icosad<Tuple<Apfloat>> face_norms_hex;

    // 12 dodecagonal faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_dec;
    private final Dodecad<Tuple<Apfloat>> face_norms_dec;

    // Apfloat Constants
    protected final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    protected final Apfloat N0 = new Apfloat("0", super.precision);
    protected final Apfloat N1 = new Apfloat("1", super.precision);
    protected final Apfloat N2 = new Apfloat("2", super.precision);
    protected final Apfloat N3 = new Apfloat("3", super.precision);
    protected final Apfloat N4 = new Apfloat("4", super.precision);
    protected final Apfloat N5 = new Apfloat("5", super.precision);
    protected final Apfloat N6 = new Apfloat("6", super.precision);
    protected final Apfloat N7 = new Apfloat("7", super.precision);
    protected final Apfloat SQRT5 = ApfloatMath.sqrt(new Apfloat("5", super.precision));
    protected final Apfloat HALF = new Apfloat("0.5", super.precision);
    protected final Apfloat SQRT_2_TIMES_5_PLUS_SQRT5 = ApfloatMath.sqrt(N2.multiply(N5.add(SQRT5)));
    protected final Apfloat SQRT_5_PLUS_2_SQRT5 = ApfloatMath.sqrt(N5.add(N2.multiply(SQRT5)));
    protected final Apfloat SQRT3 = ApfloatMath.sqrt(new Apfloat("3", precision));
    protected final Apfloat SQRT15 = ApfloatMath.sqrt(new Apfloat("15", precision));
    protected final Apfloat SQRT_2_TIMES_5_MINUS_SQRT5 = ApfloatMath.sqrt(N2.multiply(N5.subtract(SQRT5)));


    public IcosidodecahedronTruncated(
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

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0 = mult(normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1 = mult(normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2 = mult(normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3 = mult(normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4 = mult(normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5 = mult(normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6 = mult(normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7 = mult(normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8 = mult(normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9 = mult(normalize(vB9), super.getRadius().toString());
        Triad<Apfloat> vC10 = mult(normalize(vB10), super.getRadius().toString());
        Triad<Apfloat> vC11 = mult(normalize(vB11), super.getRadius().toString());
        Triad<Apfloat> vC12 = mult(normalize(vB12), super.getRadius().toString());
        Triad<Apfloat> vC13 = mult(normalize(vB13), super.getRadius().toString());
        Triad<Apfloat> vC14 = mult(normalize(vB14), super.getRadius().toString());
        Triad<Apfloat> vC15 = mult(normalize(vB15), super.getRadius().toString());
        Triad<Apfloat> vC16 = mult(normalize(vB16), super.getRadius().toString());
        Triad<Apfloat> vC17 = mult(normalize(vB17), super.getRadius().toString());
        Triad<Apfloat> vC18 = mult(normalize(vB18), super.getRadius().toString());
        Triad<Apfloat> vC19 = mult(normalize(vB19), super.getRadius().toString());
        Triad<Apfloat> vC20 = mult(normalize(vB20), super.getRadius().toString());
        Triad<Apfloat> vC21 = mult(normalize(vB21), super.getRadius().toString());
        Triad<Apfloat> vC22 = mult(normalize(vB22), super.getRadius().toString());
        Triad<Apfloat> vC23 = mult(normalize(vB23), super.getRadius().toString());
        Triad<Apfloat> vC24 = mult(normalize(vB24), super.getRadius().toString());
        Triad<Apfloat> vC25 = mult(normalize(vB25), super.getRadius().toString());
        Triad<Apfloat> vC26 = mult(normalize(vB26), super.getRadius().toString());
        Triad<Apfloat> vC27 = mult(normalize(vB27), super.getRadius().toString());
        Triad<Apfloat> vC28 = mult(normalize(vB28), super.getRadius().toString());
        Triad<Apfloat> vC29 = mult(normalize(vB29), super.getRadius().toString());
        Triad<Apfloat> vC30 = mult(normalize(vB30), super.getRadius().toString());
        Triad<Apfloat> vC31 = mult(normalize(vB31), super.getRadius().toString());
        Triad<Apfloat> vC32 = mult(normalize(vB32), super.getRadius().toString());
        Triad<Apfloat> vC33 = mult(normalize(vB33), super.getRadius().toString());
        Triad<Apfloat> vC34 = mult(normalize(vB34), super.getRadius().toString());
        Triad<Apfloat> vC35 = mult(normalize(vB35), super.getRadius().toString());
        Triad<Apfloat> vC36 = mult(normalize(vB36), super.getRadius().toString());
        Triad<Apfloat> vC37 = mult(normalize(vB37), super.getRadius().toString());
        Triad<Apfloat> vC38 = mult(normalize(vB38), super.getRadius().toString());
        Triad<Apfloat> vC39 = mult(normalize(vB39), super.getRadius().toString());
        Triad<Apfloat> vC40 = mult(normalize(vB40), super.getRadius().toString());
        Triad<Apfloat> vC41 = mult(normalize(vB41), super.getRadius().toString());
        Triad<Apfloat> vC42 = mult(normalize(vB42), super.getRadius().toString());
        Triad<Apfloat> vC43 = mult(normalize(vB43), super.getRadius().toString());
        Triad<Apfloat> vC44 = mult(normalize(vB44), super.getRadius().toString());
        Triad<Apfloat> vC45 = mult(normalize(vB45), super.getRadius().toString());
        Triad<Apfloat> vC46 = mult(normalize(vB46), super.getRadius().toString());
        Triad<Apfloat> vC47 = mult(normalize(vB47), super.getRadius().toString());
        Triad<Apfloat> vC48 = mult(normalize(vB48), super.getRadius().toString());
        Triad<Apfloat> vC49 = mult(normalize(vB49), super.getRadius().toString());
        Triad<Apfloat> vC50 = mult(normalize(vB50), super.getRadius().toString());
        Triad<Apfloat> vC51 = mult(normalize(vB51), super.getRadius().toString());
        Triad<Apfloat> vC52 = mult(normalize(vB52), super.getRadius().toString());
        Triad<Apfloat> vC53 = mult(normalize(vB53), super.getRadius().toString());
        Triad<Apfloat> vC54 = mult(normalize(vB54), super.getRadius().toString());
        Triad<Apfloat> vC55 = mult(normalize(vB55), super.getRadius().toString());
        Triad<Apfloat> vC56 = mult(normalize(vB56), super.getRadius().toString());
        Triad<Apfloat> vC57 = mult(normalize(vB57), super.getRadius().toString());
        Triad<Apfloat> vC58 = mult(normalize(vB58), super.getRadius().toString());
        Triad<Apfloat> vC59 = mult(normalize(vB59), super.getRadius().toString());
        Triad<Apfloat> vC60 = mult(normalize(vB60), super.getRadius().toString());
        Triad<Apfloat> vC61 = mult(normalize(vB61), super.getRadius().toString());
        Triad<Apfloat> vC62 = mult(normalize(vB62), super.getRadius().toString());
        Triad<Apfloat> vC63 = mult(normalize(vB63), super.getRadius().toString());
        Triad<Apfloat> vC64 = mult(normalize(vB64), super.getRadius().toString());
        Triad<Apfloat> vC65 = mult(normalize(vB65), super.getRadius().toString());
        Triad<Apfloat> vC66 = mult(normalize(vB66), super.getRadius().toString());
        Triad<Apfloat> vC67 = mult(normalize(vB67), super.getRadius().toString());
        Triad<Apfloat> vC68 = mult(normalize(vB68), super.getRadius().toString());
        Triad<Apfloat> vC69 = mult(normalize(vB69), super.getRadius().toString());
        Triad<Apfloat> vC70 = mult(normalize(vB70), super.getRadius().toString());
        Triad<Apfloat> vC71 = mult(normalize(vB71), super.getRadius().toString());
        Triad<Apfloat> vC72 = mult(normalize(vB72), super.getRadius().toString());
        Triad<Apfloat> vC73 = mult(normalize(vB73), super.getRadius().toString());
        Triad<Apfloat> vC74 = mult(normalize(vB74), super.getRadius().toString());
        Triad<Apfloat> vC75 = mult(normalize(vB75), super.getRadius().toString());
        Triad<Apfloat> vC76 = mult(normalize(vB76), super.getRadius().toString());
        Triad<Apfloat> vC77 = mult(normalize(vB77), super.getRadius().toString());
        Triad<Apfloat> vC78 = mult(normalize(vB78), super.getRadius().toString());
        Triad<Apfloat> vC79 = mult(normalize(vB79), super.getRadius().toString());
        Triad<Apfloat> vC80 = mult(normalize(vB80), super.getRadius().toString());
        Triad<Apfloat> vC81 = mult(normalize(vB81), super.getRadius().toString());
        Triad<Apfloat> vC82 = mult(normalize(vB82), super.getRadius().toString());
        Triad<Apfloat> vC83 = mult(normalize(vB83), super.getRadius().toString());
        Triad<Apfloat> vC84 = mult(normalize(vB84), super.getRadius().toString());
        Triad<Apfloat> vC85 = mult(normalize(vB85), super.getRadius().toString());
        Triad<Apfloat> vC86 = mult(normalize(vB86), super.getRadius().toString());
        Triad<Apfloat> vC87 = mult(normalize(vB87), super.getRadius().toString());
        Triad<Apfloat> vC88 = mult(normalize(vB88), super.getRadius().toString());
        Triad<Apfloat> vC89 = mult(normalize(vB89), super.getRadius().toString());
        Triad<Apfloat> vC90 = mult(normalize(vB90), super.getRadius().toString());
        Triad<Apfloat> vC91 = mult(normalize(vB91), super.getRadius().toString());
        Triad<Apfloat> vC92 = mult(normalize(vB92), super.getRadius().toString());
        Triad<Apfloat> vC93 = mult(normalize(vB93), super.getRadius().toString());
        Triad<Apfloat> vC94 = mult(normalize(vB94), super.getRadius().toString());
        Triad<Apfloat> vC95 = mult(normalize(vB95), super.getRadius().toString());
        Triad<Apfloat> vC96 = mult(normalize(vB96), super.getRadius().toString());
        Triad<Apfloat> vC97 = mult(normalize(vB97), super.getRadius().toString());
        Triad<Apfloat> vC98 = mult(normalize(vB98), super.getRadius().toString());
        Triad<Apfloat> vC99 = mult(normalize(vB99), super.getRadius().toString());
        Triad<Apfloat> vC100 = mult(normalize(vB100), super.getRadius().toString());
        Triad<Apfloat> vC101 = mult(normalize(vB101), super.getRadius().toString());
        Triad<Apfloat> vC102 = mult(normalize(vB102), super.getRadius().toString());
        Triad<Apfloat> vC103 = mult(normalize(vB103), super.getRadius().toString());
        Triad<Apfloat> vC104 = mult(normalize(vB104), super.getRadius().toString());
        Triad<Apfloat> vC105 = mult(normalize(vB105), super.getRadius().toString());
        Triad<Apfloat> vC106 = mult(normalize(vB106), super.getRadius().toString());
        Triad<Apfloat> vC107 = mult(normalize(vB107), super.getRadius().toString());
        Triad<Apfloat> vC108 = mult(normalize(vB108), super.getRadius().toString());
        Triad<Apfloat> vC109 = mult(normalize(vB109), super.getRadius().toString());
        Triad<Apfloat> vC110 = mult(normalize(vB110), super.getRadius().toString());
        Triad<Apfloat> vC111 = mult(normalize(vB111), super.getRadius().toString());
        Triad<Apfloat> vC112 = mult(normalize(vB112), super.getRadius().toString());
        Triad<Apfloat> vC113 = mult(normalize(vB113), super.getRadius().toString());
        Triad<Apfloat> vC114 = mult(normalize(vB114), super.getRadius().toString());
        Triad<Apfloat> vC115 = mult(normalize(vB115), super.getRadius().toString());
        Triad<Apfloat> vC116 = mult(normalize(vB116), super.getRadius().toString());
        Triad<Apfloat> vC117 = mult(normalize(vB117), super.getRadius().toString());
        Triad<Apfloat> vC118 = mult(normalize(vB118), super.getRadius().toString());
        Triad<Apfloat> vC119 = mult(normalize(vB119), super.getRadius().toString());

        // ==== QUADRILATERAL FACES ====
        faces_sqr.add(new Tetrad<>(vC0, vC4, vC6, vC2));
        face_norms_sqr.add(normalQuad(vC0, vC4, vC6, vC2, true));
        faces_sqr.add(new Tetrad<>(vC1, vC3, vC7, vC5));
        face_norms_sqr.add(normalQuad(vC1, vC3, vC7, vC5, true));
        faces_sqr.add(new Tetrad<>(vC8, vC10, vC11, vC9));
        face_norms_sqr.add(normalQuad(vC8, vC10, vC11, vC9, true));         // { 8, 10, 11, 9 }
        faces_sqr.add(new Tetrad<>(vC12, vC13, vC15, vC14));
        face_norms_sqr.add(normalQuad(vC12, vC13, vC15, vC14, true));       // { 12, 13, 15, 14 }
        faces_sqr.add(new Tetrad<>(vC16, vC17, vC21, vC20));
        face_norms_sqr.add(normalQuad(vC16, vC17, vC21, vC20, true));       // { 16, 17, 21, 20 }
        faces_sqr.add(new Tetrad<>(vC18, vC22, vC23, vC19));
        face_norms_sqr.add(normalQuad(vC18, vC22, vC23, vC19, true));       // { 18, 22, 23, 19 }
        faces_sqr.add(new Tetrad<>(vC24, vC72, vC96, vC48));
        face_norms_sqr.add(normalQuad(vC24, vC72, vC96, vC48, true));       // { 24, 72, 96, 48 }
        faces_sqr.add(new Tetrad<>(vC25, vC49, vC97, vC73));
        face_norms_sqr.add(normalQuad(vC25, vC49, vC97, vC73, true));       // { 25, 49, 97, 73 }
        faces_sqr.add(new Tetrad<>(vC26, vC50, vC98, vC74));
        face_norms_sqr.add(normalQuad(vC26, vC50, vC98, vC74, true));       // { 26, 50, 98, 74 }
        faces_sqr.add(new Tetrad<>(vC27, vC75, vC99, vC51));
        face_norms_sqr.add(normalQuad(vC27, vC75, vC99, vC51, true));       // { 27, 75, 99, 51 }
        faces_sqr.add(new Tetrad<>(vC28, vC52, vC100, vC76));
        face_norms_sqr.add(normalQuad(vC28, vC52, vC100, vC76, true));      // { 28, 52, 100, 76 }
        faces_sqr.add(new Tetrad<>(vC29, vC77, vC101, vC53));
        face_norms_sqr.add(normalQuad(vC29, vC77, vC101, vC53, true));      // { 29, 77, 101, 53 }
        faces_sqr.add(new Tetrad<>(vC30, vC78, vC102, vC54));
        face_norms_sqr.add(normalQuad(vC30, vC78, vC102, vC54, true));      // { 30, 78, 102, 54 }
        faces_sqr.add(new Tetrad<>(vC31, vC55, vC103, vC79));
        face_norms_sqr.add(normalQuad(vC31, vC55, vC103, vC79, true));      // { 31, 55, 103, 79 }
        faces_sqr.add(new Tetrad<>(vC32, vC80, vC104, vC56));
        face_norms_sqr.add(normalQuad(vC32, vC80, vC104, vC56, true));      // { 32, 80, 104, 56 }
        faces_sqr.add(new Tetrad<>(vC33, vC57, vC105, vC81));
        face_norms_sqr.add(normalQuad(vC33, vC57, vC105, vC81, true));      // { 33, 57, 105, 81 }
        faces_sqr.add(new Tetrad<>(vC34, vC58, vC106, vC82));
        face_norms_sqr.add(normalQuad(vC34, vC58, vC106, vC82, true));      // { 34, 58, 106, 82 }
        faces_sqr.add(new Tetrad<>(vC35, vC83, vC107, vC59));
        face_norms_sqr.add(normalQuad(vC35, vC83, vC107, vC59, true));      // { 35, 83, 107, 59 }
        faces_sqr.add(new Tetrad<>(vC36, vC60, vC108, vC84));
        face_norms_sqr.add(normalQuad(vC36, vC60, vC108, vC84, true));      // { 36, 60, 108, 84 }
        faces_sqr.add(new Tetrad<>(vC37, vC85, vC109, vC61));
        face_norms_sqr.add(normalQuad(vC37, vC85, vC109, vC61, true));      // { 37, 85, 109, 61 }
        faces_sqr.add(new Tetrad<>(vC38, vC86, vC110, vC62));
        face_norms_sqr.add(normalQuad(vC38, vC86, vC110, vC62, true));      // { 38, 86, 110, 62 }
        faces_sqr.add(new Tetrad<>(vC39, vC63, vC111, vC87));
        face_norms_sqr.add(normalQuad(vC39, vC63, vC111, vC87, true));      // { 39, 63, 111, 87 }
        faces_sqr.add(new Tetrad<>(vC40, vC88, vC112, vC64));
        face_norms_sqr.add(normalQuad(vC40, vC88, vC112, vC64, true));      // { 40, 88, 112, 64 }
        faces_sqr.add(new Tetrad<>(vC41, vC65, vC113, vC89));
        face_norms_sqr.add(normalQuad(vC41, vC65, vC113, vC89, true));      // { 41, 65, 113, 89 }
        faces_sqr.add(new Tetrad<>(vC42, vC66, vC114, vC90));
        face_norms_sqr.add(normalQuad(vC42, vC66, vC114, vC90, true));      // { 42, 66, 114, 90 }
        faces_sqr.add(new Tetrad<>(vC43, vC91, vC115, vC67));
        face_norms_sqr.add(normalQuad(vC43, vC91, vC115, vC67, true));      // { 43, 91, 115, 67 }
        faces_sqr.add(new Tetrad<>(vC44, vC68, vC116, vC92));
        face_norms_sqr.add(normalQuad(vC44, vC68, vC116, vC92, true));      // { 44, 68, 116, 92 }
        faces_sqr.add(new Tetrad<>(vC45, vC93, vC117, vC69));
        face_norms_sqr.add(normalQuad(vC45, vC93, vC117, vC69, true));      // { 45, 93, 117, 69 }
        faces_sqr.add(new Tetrad<>(vC46, vC94, vC118, vC70));
        face_norms_sqr.add(normalQuad(vC46, vC94, vC118, vC70, true));      // { 46, 94, 118, 70 }
        faces_sqr.add(new Tetrad<>(vC47, vC71, vC119, vC95));
        face_norms_sqr.add(normalQuad(vC47, vC71, vC119, vC95, true));      // { 47, 71, 119, 95 }

        // ====  HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC0, vC24, vC48, vC52, vC28, vC4);
        Triad<Apfloat> hex0_norm = normalHex(vC0, vC24, vC48, vC52, vC28, vC4, true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC2, vC6, vC30, vC54, vC50, vC26);
        Triad<Apfloat> hex1_norm = normalHex(vC2, vC6, vC30, vC54, vC50, vC26, true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC3, vC27, vC51, vC55, vC31, vC7);
        Triad<Apfloat> hex2_norm = normalHex(vC3, vC27, vC51, vC55, vC31, vC7, true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC8, vC32, vC56, vC58, vC34, vC10);
        Triad<Apfloat> hex3_norm = normalHex(vC8, vC32, vC56, vC58, vC34, vC10, true);
        Hexad<Tuple<Apfloat>> hex4 = new Hexad<>(vC9, vC11, vC35, vC59, vC57, vC33);
        Triad<Apfloat> hex4_norm = normalHex(vC9, vC11, vC35, vC59, vC57, vC33, true);
        Hexad<Tuple<Apfloat>> hex5 = new Hexad<>(vC12, vC14, vC38, vC62, vC60, vC36);
        Triad<Apfloat> hex5_norm = normalHex(vC12, vC14, vC38, vC62, vC60, vC36, true);
        Hexad<Tuple<Apfloat>> hex6 = new Hexad<>(vC13, vC37, vC61, vC63, vC39, vC15);
        Triad<Apfloat> hex6_norm = normalHex(vC13, vC37, vC61, vC63, vC39, vC15, true);
        Hexad<Tuple<Apfloat>> hex7 = new Hexad<>(vC16, vC40, vC64, vC65, vC41, vC17);
        Triad<Apfloat> hex7_norm = normalHex(vC16, vC40, vC64, vC65, vC41, vC17, true);
        Hexad<Tuple<Apfloat>> hex8 = new Hexad<>(vC18, vC19, vC43, vC67, vC66, vC42);
        Triad<Apfloat> hex8_norm = normalHex(vC18, vC19, vC43, vC67, vC66, vC42, true);
        Hexad<Tuple<Apfloat>> hex9 = new Hexad<>(vC20, vC21, vC45, vC69, vC68, vC44);
        Triad<Apfloat> hex9_norm = normalHex(vC20, vC21, vC45, vC69, vC68, vC44, true);
        Hexad<Tuple<Apfloat>> hex10 = new Hexad<>(vC22, vC46, vC70, vC71, vC47, vC23);
        Triad<Apfloat> hex10_norm = normalHex(vC22, vC46, vC70, vC71, vC47, vC23, true);
        Hexad<Tuple<Apfloat>> hex11 = new Hexad<>(vC72, vC104, vC80, vC112, vC88, vC96);
        Triad<Apfloat> hex11_norm = normalHex(vC72, vC104, vC80, vC112, vC88, vC96, true);
        Hexad<Tuple<Apfloat>> hex12 = new Hexad<>(vC73, vC97, vC89, vC113, vC81, vC105);
        Triad<Apfloat> hex12_norm = normalHex(vC73, vC97, vC89, vC113, vC81, vC105, true);
        Hexad<Tuple<Apfloat>> hex13 = new Hexad<>(vC74, vC98, vC90, vC114, vC82, vC106);
        Triad<Apfloat> hex13_norm = normalHex(vC74, vC98, vC90, vC114, vC82, vC106, true);
        Hexad<Tuple<Apfloat>> hex14 = new Hexad<>(vC75, vC107, vC83, vC115, vC91, vC99);
        Triad<Apfloat> hex14_norm = normalHex(vC75, vC107, vC83, vC115, vC91, vC99, true);
        Hexad<Tuple<Apfloat>> hex15 = new Hexad<>(vC76, vC100, vC92, vC116, vC84, vC108);
        Triad<Apfloat> hex15_norm = normalHex(vC76, vC100, vC92, vC116, vC84, vC108, true);
        Hexad<Tuple<Apfloat>> hex16 = new Hexad<>(vC77, vC109, vC85, vC117, vC93, vC101);
        Triad<Apfloat> hex16_norm = normalHex(vC77, vC109, vC85, vC117, vC93, vC101, true);
        Hexad<Tuple<Apfloat>> hex17 = new Hexad<>(vC78, vC110, vC86, vC118, vC94, vC102);
        Triad<Apfloat> hex17_norm = normalHex(vC78, vC110, vC86, vC118, vC94, vC102, true);
        Hexad<Tuple<Apfloat>> hex18 = new Hexad<>(vC79, vC103, vC95, vC119, vC87, vC111);
        Triad<Apfloat> hex18_norm = normalHex(vC79, vC103, vC95, vC119, vC87, vC111, true);
        Hexad<Tuple<Apfloat>> hex19 = new Hexad<>(vC1, vC5, vC29, vC53, vC49, vC25);
        Triad<Apfloat> hex19_norm = normalHex(vC1, vC5, vC29, vC53, vC49, vC25, true);
        faces_hex = new Icosad<>(
        hex0, hex1, hex2, hex3,
        hex4, hex5, hex6, hex7,
        hex8, hex9, hex10, hex11,
        hex12, hex13, hex14, hex15,
        hex16, hex17, hex18, hex19
            );
        face_norms_hex = new Icosad<>(
        hex0_norm, hex1_norm, hex2_norm, hex3_norm,
        hex4_norm, hex5_norm, hex6_norm, hex7_norm,
        hex8_norm, hex9_norm, hex10_norm, hex11_norm,
        hex12_norm, hex13_norm, hex14_norm, hex15_norm,
        hex16_norm, hex17_norm, hex18_norm, hex19_norm
        );
        // ==== DODECAGONAL FACES ====
        Decad<Tuple<Apfloat>> dec0 = new Decad<>(vC0, vC2, vC26, vC74, vC106, vC58, vC56, vC104, vC72, vC24);
        Triad<Apfloat> dec0_norm = normalDeca(vC0, vC2, vC26, vC74, vC106, vC58, vC56, vC104, vC72, vC24, true);
        Decad<Tuple<Apfloat>> dec1 = new Decad<>(vC4, vC28, vC76, vC108, vC60, vC62, vC110, vC78, vC30, vC6);
        Triad<Apfloat> dec1_norm = normalDeca(vC4, vC28, vC76, vC108, vC60, vC62, vC110, vC78, vC30, vC6, true);
        Decad<Tuple<Apfloat>> dec2 = new Decad<>(vC5, vC7, vC31, vC79, vC111, vC63, vC61, vC109, vC77, vC29);
        Triad<Apfloat> dec2_norm = normalDeca(vC5, vC7, vC31, vC79, vC111, vC63, vC61, vC109, vC77, vC29, true);
        Decad<Tuple<Apfloat>> dec3 = new Decad<>(vC8, vC9, vC33, vC81, vC113, vC65, vC64, vC112, vC80, vC32);
        Triad<Apfloat> dec3_norm = normalDeca(vC8, vC9, vC33, vC81, vC113, vC65, vC64, vC112, vC80, vC32, true);
        Decad<Tuple<Apfloat>> dec4 = new Decad<>(vC10, vC34, vC82, vC114, vC66, vC67, vC115, vC83, vC35, vC11);
        Triad<Apfloat> dec4_norm = normalDeca(vC10, vC34, vC82, vC114, vC66, vC67, vC115, vC83, vC35, vC11, true);
        Decad<Tuple<Apfloat>> dec5 = new Decad<>(vC12, vC36, vC84, vC116, vC68, vC69, vC117, vC85, vC37, vC13);
        Triad<Apfloat> dec5_norm = normalDeca(vC12, vC36, vC84, vC116, vC68, vC69, vC117, vC85, vC37, vC13, true);
        Decad<Tuple<Apfloat>> dec6 = new Decad<>(vC14, vC15, vC39, vC87, vC119, vC71, vC70, vC118, vC86, vC38);
        Triad<Apfloat> dec6_norm = normalDeca(vC14, vC15, vC39, vC87, vC119, vC71, vC70, vC118, vC86, vC38, true);
        Decad<Tuple<Apfloat>> dec7 = new Decad<>(vC16, vC20, vC44, vC92, vC100, vC52, vC48, vC96, vC88, vC40);
        Triad<Apfloat> dec7_norm = normalDeca(vC16, vC20, vC44, vC92, vC100, vC52, vC48, vC96, vC88, vC40, true);
        Decad<Tuple<Apfloat>> dec8 = new Decad<>(vC17, vC41, vC89, vC97, vC49, vC53, vC101, vC93, vC45, vC21);
        Triad<Apfloat> dec8_norm = normalDeca(vC17, vC41, vC89, vC97, vC49, vC53, vC101, vC93, vC45, vC21, true);
        Decad<Tuple<Apfloat>> dec9 = new Decad<>(vC18, vC42, vC90, vC98, vC50, vC54, vC102, vC94, vC46, vC22);
        Triad<Apfloat> dec9_norm = normalDeca(vC18, vC42, vC90, vC98, vC50, vC54, vC102, vC94, vC46, vC22, true);
        Decad<Tuple<Apfloat>> dec10 = new Decad<>(vC19, vC23, vC47, vC95, vC103, vC55, vC51, vC99, vC91, vC43);
        Triad<Apfloat> dec10_norm = normalDeca(vC19, vC23, vC47, vC95, vC103, vC55, vC51, vC99, vC91, vC43, true);
        Decad<Tuple<Apfloat>> dec11 = new Decad<>(vC0, vC2, vC26, vC74, vC106, vC58, vC56, vC104, vC72, vC24);
        Triad<Apfloat> dec11_norm = normalDeca(vC0, vC2, vC26, vC74, vC106, vC58, vC56, vC104, vC72, vC24, true);
        faces_dec = new Dodecad<>(
                dec0, dec1, dec2, dec3,
                dec4, dec5, dec6, dec7,
                dec8, dec9, dec10, dec11
        );
        face_norms_dec = new Dodecad<>(
                dec0_norm, dec1_norm, dec2_norm, dec3_norm,
                dec4_norm, dec5_norm, dec6_norm, dec7_norm,
                dec8_norm, dec9_norm, dec10_norm, dec11_norm
        );
}

    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_sqr.size(); i++){

            // === Check dodecagonal faces ===
            if (i < faces_dec.fetchSize()) {
                Decad<Tuple<Apfloat>> face = (Decad<Tuple<Apfloat>>) faces_dec.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //  = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //  = face_norm_x*(p_x-vertA_x)
                //  + face_norm_y*(p_y-vertA_y)
                //  + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_dec.fetch(i);
                Apfloat d = dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check hexagonal faces ===
            if (i < faces_hex.fetchSize()) {
                Hexad<Tuple<Apfloat>> face = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);

                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //  = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //  = face_norm_x*(p_x-vertA_x)
                //  + face_norm_y*(p_y-vertA_y)
                //  + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_hex.fetch(i);
                Apfloat d = dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check square faces ===
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.get(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);

            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
            Triad<Apfloat> m = subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //  = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //  = face_norm_x*(p_x-vertA_x)
            //  + face_norm_y*(p_y-vertA_y)
            //  + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_sqr.get(i);
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }

        }
        return true;
    }
}
