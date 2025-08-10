package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.*;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b>Biscribed Truncated Cuboctahedron</b>
 * <p> A polyhedron where all faces
 * share the same inradius, meaning each face can be inscribed in a sphere of
 * equal radius.
 * <p> This shape is derived from the {@link CuboctahedronTruncated} by adjusting its
 * geometry so that all faces are tangent to a common sphere. It retains the
 * overall symmetry of the original while altering edge lengths to achieve the
 * biscribed property.
 * @see <a href="https://dmccooey.com/polyhedra/BiscribedTruncatedCuboctahedron.html">
 * Biscribed Truncated Cuboctahedron – dmccooey.com</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class CuboctahedronTruncatedBiscribed extends Shape {

    // 12 Quadrilateral Faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_sqr;
    private final Dodecad<Tuple<Apfloat>> face_norms_sqr;

    // 8 Hexagonal Faces
    private final Octad<Tuple<Tuple<Apfloat>>> faces_hex;
    private final Octad<Tuple<Apfloat>> face_norms_hex;

    // 6 octagonal Faces
    private final Hexad<Tuple<Tuple<Apfloat>>> faces_oct;
    private final Hexad<Tuple<Apfloat>> face_norms_oct;

    // Apfloat Constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat SQRT2 = ApfloatMath.sqrt(new Apfloat("2", super.precision));
    private final Apfloat SQRT3 = ApfloatMath.sqrt(new Apfloat("3", super.precision));

    // C0 = (1 + sqrt(2)) / 2
    private final Apfloat C0 = SQRT3.subtract(SQRT2);

    // C1 = (1 + 2 * sqrt(2)) / 2
    private final Apfloat C1 = SQRT2.subtract(N1);

    public CuboctahedronTruncatedBiscribed(
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
        Triad<Apfloat> vB0 = new Triad<>(C1, C0, N1);
        Triad<Apfloat> vB1 = new Triad<>(C1, C0, NEG_N1);
        Triad<Apfloat> vB2 = new Triad<>(C1, C0.multiply(NEG_N1), N1);
        Triad<Apfloat> vB3 = new Triad<>(C1, C0.multiply(NEG_N1), NEG_N1);
        Triad<Apfloat> vB4 = new Triad<>(C1.multiply(NEG_N1), C0, N1);
        Triad<Apfloat> vB5 = new Triad<>(C1.multiply(NEG_N1), C0, NEG_N1);
        Triad<Apfloat> vB6 = new Triad<>(C1.multiply(NEG_N1), C0.multiply(NEG_N1), N1);
        Triad<Apfloat> vB7 = new Triad<>(C1.multiply(NEG_N1), C0.multiply(NEG_N1), NEG_N1);
        Triad<Apfloat> vB8 = new Triad<>(N1, C1, C0);
        Triad<Apfloat> vB9 = new Triad<>(N1, C1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB10 = new Triad<>(N1, C1.multiply(NEG_N1), C0);
        Triad<Apfloat> vB11 = new Triad<>(N1, C1.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB12 = new Triad<>(NEG_N1, C1, C0);
        Triad<Apfloat> vB13 = new Triad<>(NEG_N1, C1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB14 = new Triad<>(NEG_N1, C1.multiply(NEG_N1), C0);
        Triad<Apfloat> vB15 = new Triad<>(NEG_N1, C1.multiply(NEG_N1), C0.multiply(NEG_N1));
        Triad<Apfloat> vB16 = new Triad<>(C0, N1, C1);
        Triad<Apfloat> vB17 = new Triad<>(C0, N1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB18 = new Triad<>(C0, NEG_N1, C1);
        Triad<Apfloat> vB19 = new Triad<>(C0, NEG_N1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB20 = new Triad<>(C0.multiply(NEG_N1), N1, C1);
        Triad<Apfloat> vB21 = new Triad<>(C0.multiply(NEG_N1), N1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB22 = new Triad<>(C0.multiply(NEG_N1), NEG_N1, C1);
        Triad<Apfloat> vB23 = new Triad<>(C0.multiply(NEG_N1), NEG_N1, C1.multiply(NEG_N1));
        Triad<Apfloat> vB24 = new Triad<>(C0, C1, N1);
        Triad<Apfloat> vB25 = new Triad<>(C0, C1, NEG_N1);
        Triad<Apfloat> vB26 = new Triad<>(C0, C1.multiply(NEG_N1), N1);
        Triad<Apfloat> vB27 = new Triad<>(C0, C1.multiply(NEG_N1), NEG_N1);
        Triad<Apfloat> vB28 = new Triad<>(C0.multiply(NEG_N1), C1, N1);
        Triad<Apfloat> vB29 = new Triad<>(C0.multiply(NEG_N1), C1, NEG_N1);
        Triad<Apfloat> vB30 = new Triad<>(C0.multiply(NEG_N1), C1.multiply(NEG_N1), N1);
        Triad<Apfloat> vB31 = new Triad<>(C0.multiply(NEG_N1), C1.multiply(NEG_N1), NEG_N1);
        Triad<Apfloat> vB32 = new Triad<>(N1, C0, C1);
        Triad<Apfloat> vB33 = new Triad<>(N1, C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB34 = new Triad<>(N1, C0.multiply(NEG_N1), C1);
        Triad<Apfloat> vB35 = new Triad<>(N1, C0.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB36 = new Triad<>(NEG_N1, C0, C1);
        Triad<Apfloat> vB37 = new Triad<>(NEG_N1, C0, C1.multiply(NEG_N1));
        Triad<Apfloat> vB38 = new Triad<>(NEG_N1, C0.multiply(NEG_N1), C1);
        Triad<Apfloat> vB39 = new Triad<>(NEG_N1, C0.multiply(NEG_N1), C1.multiply(NEG_N1));
        Triad<Apfloat> vB40 = new Triad<>(C1, N1, C0);
        Triad<Apfloat> vB41 = new Triad<>(C1, N1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB42 = new Triad<>(C1, NEG_N1, C0);
        Triad<Apfloat> vB43 = new Triad<>(C1, NEG_N1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB44 = new Triad<>(C1.multiply(NEG_N1), N1, C0);
        Triad<Apfloat> vB45 = new Triad<>(C1.multiply(NEG_N1), N1, C0.multiply(NEG_N1));
        Triad<Apfloat> vB46 = new Triad<>(C1.multiply(NEG_N1), NEG_N1, C0);
        Triad<Apfloat> vB47 = new Triad<>(C1.multiply(NEG_N1), NEG_N1, C0.multiply(NEG_N1));

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0  = mult(normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1  = mult(normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2  = mult(normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3  = mult(normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4  = mult(normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5  = mult(normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6  = mult(normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7  = mult(normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8  = mult(normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9  = mult(normalize(vB9), super.getRadius().toString());
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

        // ==== OCTAGONAL FACES ====
        Octad<Tuple<Apfloat>> oct0 = new Octad<>(vC0,vC24,vC28,vC4,vC6,vC30,vC26,vC2);
        Triad<Apfloat> oct0_norm = normalOct(vC0,vC24,vC28,vC4,vC6,vC30,vC26,vC2,true);
        Octad<Tuple<Apfloat>> oct1 = new Octad<>(vC1,vC3,vC27,vC31,vC7,vC5,vC29,vC25);
        Triad<Apfloat> oct1_norm = normalOct(vC1,vC3,vC27,vC31,vC7,vC5,vC29,vC25,true);
        Octad<Tuple<Apfloat>> oct2 = new Octad<>(vC8,vC32,vC34,vC10,vC11,vC35,vC33,vC9);
        Triad<Apfloat> oct2_norm = normalOct(vC8,vC32,vC34,vC10,vC11,vC35,vC33,vC9,true);
        Octad<Tuple<Apfloat>> oct3 = new Octad<>(vC12,vC13,vC37,vC39,vC15,vC14,vC38,vC36);
        Triad<Apfloat> oct3_norm = normalOct(vC12,vC13,vC37,vC39,vC15,vC14,vC38,vC36,true);
        Octad<Tuple<Apfloat>> oct4 = new Octad<>(vC16,vC40,vC41,vC17,vC21,vC45,vC44,vC20);
        Triad<Apfloat> oct4_norm = normalOct(vC16,vC40,vC41,vC17,vC21,vC45,vC44,vC20,true);
        Octad<Tuple<Apfloat>> oct5 = new Octad<>(vC18,vC22,vC46,vC47,vC23,vC19,vC43,vC42);
        Triad<Apfloat> oct5_norm = normalOct(vC18,vC22,vC46,vC47,vC23,vC19,vC43,vC42,true);
        faces_oct = new Hexad<>(oct0, oct1,oct2,oct3,oct4,oct5);
        face_norms_oct = new Hexad<>(oct0_norm,oct1_norm,oct2_norm,oct3_norm,oct4_norm,oct5_norm);

        // ==== HEXAGONAL FACES ====
        Hexad<Tuple<Apfloat>> hex0 = new Hexad<>(vC0, vC32, vC8, vC40, vC16, vC24);
        Triad<Apfloat> hex0_norm = normalHex(vC0, vC32, vC8, vC40, vC16, vC24, true);
        Hexad<Tuple<Apfloat>> hex1 = new Hexad<>(vC1, vC33, vC9, vC41, vC17, vC25);
        Triad<Apfloat> hex1_norm = normalHex(vC1, vC33, vC9, vC41, vC17, vC25, true);
        Hexad<Tuple<Apfloat>> hex2 = new Hexad<>(vC2, vC34, vC10, vC42, vC18, vC26);
        Triad<Apfloat> hex2_norm = normalHex(vC2, vC34, vC10, vC42, vC18, vC26, true);
        Hexad<Tuple<Apfloat>> hex3 = new Hexad<>(vC3, vC35, vC11, vC43, vC19, vC27);
        Triad<Apfloat> hex3_norm = normalHex(vC3, vC35, vC11, vC43, vC19, vC27, true);
        Hexad<Tuple<Apfloat>> hex4 = new Hexad<>(vC4, vC36, vC12, vC44, vC20, vC28);
        Triad<Apfloat> hex4_norm = normalHex(vC4, vC36, vC12, vC44, vC20, vC28, true);
        Hexad<Tuple<Apfloat>> hex5 = new Hexad<>(vC5, vC37, vC13, vC45, vC21, vC29);
        Triad<Apfloat> hex5_norm = normalHex(vC5, vC37, vC13, vC45, vC21, vC29, true);
        Hexad<Tuple<Apfloat>> hex6 = new Hexad<>(vC6, vC38, vC14, vC46, vC22, vC30);
        Triad<Apfloat> hex6_norm = normalHex(vC6, vC38, vC14, vC46, vC22, vC30, true);
        Hexad<Tuple<Apfloat>> hex7 = new Hexad<>(vC7, vC39, vC15, vC47, vC23, vC31);
        Triad<Apfloat> hex7_norm = normalHex(vC7, vC39, vC15, vC47, vC23, vC31, true);
        faces_hex = new Octad<>(hex0,hex1,hex2,hex3,hex4,hex5,hex6,hex7);
        face_norms_hex = new Octad<>(hex0_norm,hex1_norm,hex2_norm,hex3_norm,hex4_norm,hex5_norm,hex6_norm,hex7_norm);

        // ==== QUADRILATERAL FACES ====
        Tetrad<Tuple<Apfloat>> sqr0 = new Tetrad<>(vC0, vC2, vC34, vC3);
        Triad<Apfloat> sqr0_norm = normalQuad(vC0, vC2, vC34, vC3, true);
        Tetrad<Tuple<Apfloat>> sqr1 = new Tetrad<>(vC4, vC36, vC38, vC6);
        Triad<Apfloat> sqr1_norm = normalQuad(vC4, vC36, vC38, vC6, true);
        Tetrad<Tuple<Apfloat>> sqr2 = new Tetrad<>(vC5, vC7, vC39, vC37);
        Triad<Apfloat> sqr2_norm = normalQuad(vC5, vC7, vC39, vC37, true);
        Tetrad<Tuple<Apfloat>> sqr3 = new Tetrad<>(vC8, vC9, vC41, vC40);
        Triad<Apfloat> sqr3_norm = normalQuad(vC8, vC9, vC41, vC40, true);
        Tetrad<Tuple<Apfloat>> sqr4 = new Tetrad<>(vC10, vC42, vC43, vC11);
        Triad<Apfloat> sqr4_norm = normalQuad(vC10, vC42, vC43, vC11, true);
        Tetrad<Tuple<Apfloat>> sqr5 = new Tetrad<>(vC12, vC44, vC45, vC13);
        Triad<Apfloat> sqr5_norm = normalQuad(vC12, vC44, vC45, vC13, true);
        Tetrad<Tuple<Apfloat>> sqr6 = new Tetrad<>(vC14, vC15, vC47, vC46);
        Triad<Apfloat> sqr6_norm = normalQuad(vC14, vC15, vC47, vC46, true);
        Tetrad<Tuple<Apfloat>> sqr7 = new Tetrad<>(vC16, vC20, vC28, vC24);
        Triad<Apfloat> sqr7_norm = normalQuad(vC16, vC20, vC28, vC24, true);
        Tetrad<Tuple<Apfloat>> sqr8 = new Tetrad<>(vC17, vC25, vC29, vC21);
        Triad<Apfloat> sqr8_norm = normalQuad(vC17, vC25, vC29, vC21, true);
        Tetrad<Tuple<Apfloat>> sqr9 = new Tetrad<>(vC18, vC26, vC30, vC22);
        Triad<Apfloat> sqr9_norm = normalQuad(vC18, vC26, vC30, vC22, true);
        Tetrad<Tuple<Apfloat>> sqr10 = new Tetrad<>(vC19, vC23, vC31, vC27);
        Triad<Apfloat> sqr10_norm = normalQuad(vC19, vC23, vC31, vC27, true);
        Tetrad<Tuple<Apfloat>> sqr11 = new Tetrad<>(vC1, vC33, vC35, vC3);
        Triad<Apfloat> sqr11_norm = normalQuad(vC1, vC33, vC35, vC3, true);
        faces_sqr = new Dodecad<>(
                sqr0,sqr1,sqr2,sqr3,
                sqr4,sqr5,sqr6,sqr7,
                sqr8,sqr9,sqr10,sqr11
        );
        face_norms_sqr = new Dodecad<>(
                sqr0_norm,sqr1_norm,sqr2_norm,sqr3_norm,
                sqr4_norm,sqr5_norm,sqr6_norm,sqr7_norm,
                sqr8_norm,sqr9_norm,sqr10_norm,sqr11_norm
        );
    }


    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < faces_sqr.fetchSize(); i++){
            // === Check hexagonal faces ===
            if (i < 8) {
                Hexad<Tuple<Apfloat>> face = (Hexad<Tuple<Apfloat>>) faces_hex.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_hex.fetch(i);
                Apfloat d = dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check octagonal faces ===
            if (i < 6) {
                Octad<Tuple<Apfloat>> face = (Octad<Tuple<Apfloat>>) faces_oct.fetch(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


                // given point p and vertA, calculate vector from vertA -> p:
                // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
                Triad<Apfloat> m = subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_oct.fetch(i);
                Apfloat d = dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check square faces ===
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);
            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x -vertA_z)
            Triad<Apfloat> m = subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_sqr.fetch(i);
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }

        return true;
    }
}
