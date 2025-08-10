package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.platonic.Icosahedron;
import io.github.noshou.npg.shapes.archimedean.IcosidodecahedronTruncated;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.Polyad;
import io.github.noshou.tuple.Triad;
import io.github.noshou.tuple.Tuple;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;

import static io.github.noshou.npg.nputil.VectorMath.*;


/**
 * Represents a <b>Hexakis Icosahedron </b> (aka Disdyakis Triacontahedron).
 * <p>This polyhedron is constructed by augmenting each face of an
 * {@link Icosahedron} with a pyramid-like structure, resulting in a convex solid
 * with 80 vertices, 120 edges, and 60 triangular faces.
 * The hexakis icosahedron often appears in two common variants:
 * <ul>
 *   <li>{@link IcosahedronHexakisCanonical} is the regular, symmetric form with uniform face geometry.</li>
 *   <li>{@link IcosahedronHexakisBiscribed} is positioned between two concentric spheres—
 *       one contacting all vertices and the other contacting all face centers—
 *       producing slight variations in vertex placement and subtle geometric differences.</li>
 * </ul>
 * <p> It is the dual polyhedron of the {@link IcosidodecahedronTruncated}.
 */
@SuppressWarnings("FieldCanBeLocal")
public abstract class IcosahedronHexakis extends Shape {

    // 120 triangular faces
    private final ArrayList<Triad<Tuple<Apfloat>>> faces = new ArrayList<>();
    private final ArrayList<Triad<Apfloat>> face_norms = new ArrayList<>();

    // Apfloat constants
    protected final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    protected final Apfloat N0 = new Apfloat("0", super.precision);
    protected final Apfloat N1 = new Apfloat("-1", super.precision);
    protected final Apfloat N2 = new Apfloat("2", super.precision);
    protected final Apfloat N3 = new Apfloat("3", super.precision);
    protected final Apfloat N4 = new Apfloat("4", super.precision);
    protected final Apfloat N5 = new Apfloat("5", super.precision);
    protected final Apfloat N6 = new Apfloat("5", super.precision);
    protected final Apfloat N9 = new Apfloat("9", super.precision);
    protected final Apfloat N10 = new Apfloat("10", super.precision);
    protected final Apfloat N11 = new Apfloat("11", super.precision);
    protected final Apfloat N15 = new Apfloat("15", super.precision);
    protected final Apfloat N22 = new Apfloat("22", super.precision);
    protected final Apfloat N44 = new Apfloat("44", super.precision);
    protected final Apfloat N105 = new Apfloat("105", super.precision);
    protected final Apfloat SQRT5 = ApfloatMath.sqrt(N5);
    protected final Apfloat SQRT3 = ApfloatMath.sqrt(N3);
    protected final Apfloat SQRT15 = ApfloatMath.sqrt(N15);
    protected final Apfloat HALF = N1.divide(N2);

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

    public IcosahedronHexakis(
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
        Triad<Apfloat> vC0  = mult(normalize(vB0),  super.getRadius().toString());
        Triad<Apfloat> vC1  = mult(normalize(vB1),  super.getRadius().toString());
        Triad<Apfloat> vC2  = mult(normalize(vB2),  super.getRadius().toString());
        Triad<Apfloat> vC3  = mult(normalize(vB3),  super.getRadius().toString());
        Triad<Apfloat> vC4  = mult(normalize(vB4),  super.getRadius().toString());
        Triad<Apfloat> vC5  = mult(normalize(vB5),  super.getRadius().toString());
        Triad<Apfloat> vC6  = mult(normalize(vB6),  super.getRadius().toString());
        Triad<Apfloat> vC7  = mult(normalize(vB7),  super.getRadius().toString());
        Triad<Apfloat> vC8  = mult(normalize(vB8),  super.getRadius().toString());
        Triad<Apfloat> vC9  = mult(normalize(vB9),  super.getRadius().toString());
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

        // ==== TRIANGULAR FACES ====
        Triad<Tuple<Apfloat>> tri0 = new Triad<>(vC18, vC0, vC8);
        faces.add(tri0);
        Triad<Apfloat> tri0_norm = normalTriple(vC18, vC0, vC8, true);
        face_norms.add(tri0_norm);
        Triad<Tuple<Apfloat>> tri1 = new Triad<>(vC18, vC8, vC32);
        faces.add(tri1);
        Triad<Apfloat> tri1_norm = normalTriple(vC18, vC8, vC32, true);
        face_norms.add(tri1_norm);
        Triad<Tuple<Apfloat>> tri2 = new Triad<>(vC18, vC32, vC56);
        faces.add(tri2);
        Triad<Apfloat> tri2_norm = normalTriple(vC18, vC32, vC56, true);
        face_norms.add(tri2_norm);
        Triad<Tuple<Apfloat>> tri3 = new Triad<>(vC18, vC56, vC40);
        faces.add(tri3);
        Triad<Apfloat> tri3_norm = normalTriple(vC18, vC56, vC40, true);
        face_norms.add(tri3_norm);
        Triad<Tuple<Apfloat>> tri4 = new Triad<>(vC18, vC40, vC10);
        faces.add(tri4);
        Triad<Apfloat> tri4_norm = normalTriple(vC18, vC40, vC10, true);
        face_norms.add(tri4_norm);
        Triad<Tuple<Apfloat>> tri5 = new Triad<>(vC18, vC10, vC38);
        faces.add(tri5);
        Triad<Apfloat> tri5_norm = normalTriple(vC18, vC10, vC38, true);
        face_norms.add(tri5_norm);
        Triad<Tuple<Apfloat>> tri6 = new Triad<>(vC18, vC38, vC54);
        faces.add(tri6);
        Triad<Apfloat> tri6_norm = normalTriple(vC18, vC38, vC54, true);
        face_norms.add(tri6_norm);
        Triad<Tuple<Apfloat>> tri7 = new Triad<>(vC18, vC54, vC30);
        faces.add(tri7);
        Triad<Apfloat> tri7_norm = normalTriple(vC18, vC54, vC30, true);
        face_norms.add(tri7_norm);
        Triad<Tuple<Apfloat>> tri8 = new Triad<>(vC18, vC30, vC6);
        faces.add(tri8);
        Triad<Apfloat> tri8_norm = normalTriple(vC18, vC30, vC6, true);
        face_norms.add(tri8_norm);
        Triad<Tuple<Apfloat>> tri9 = new Triad<>(vC18, vC6, vC0);
        faces.add(tri9);
        Triad<Apfloat> tri9_norm = normalTriple(vC18, vC6, vC0, true);
        face_norms.add(tri9_norm);
        Triad<Tuple<Apfloat>> tri10 = new Triad<>(vC19, vC1, vC7);
        faces.add(tri10);
        Triad<Apfloat> tri10_norm = normalTriple(vC19, vC1, vC7, true);
        face_norms.add(tri10_norm);
        Triad<Tuple<Apfloat>> tri11 = new Triad<>(vC19, vC7, vC31);
        faces.add(tri11);
        Triad<Apfloat> tri11_norm = normalTriple(vC19, vC7, vC31, true);
        face_norms.add(tri11_norm);
        Triad<Tuple<Apfloat>> tri12 = new Triad<>(vC19, vC31, vC55);
        faces.add(tri12);
        Triad<Apfloat> tri12_norm = normalTriple(vC19, vC31, vC55, true);
        face_norms.add(tri12_norm);
        Triad<Tuple<Apfloat>> tri13 = new Triad<>(vC19, vC55, vC39);
        faces.add(tri13);
        Triad<Apfloat> tri13_norm = normalTriple(vC19, vC55, vC39, true);
        face_norms.add(tri13_norm);
        Triad<Tuple<Apfloat>> tri14 = new Triad<>(vC19, vC39, vC11);
        faces.add(tri14);
        Triad<Apfloat> tri14_norm = normalTriple(vC19, vC39, vC11, true);
        face_norms.add(tri14_norm);
        Triad<Tuple<Apfloat>> tri15 = new Triad<>(vC19, vC11, vC41);
        faces.add(tri15);
        Triad<Apfloat> tri15_norm = normalTriple(vC19, vC11, vC41, true);
        face_norms.add(tri15_norm);
        Triad<Tuple<Apfloat>> tri16 = new Triad<>(vC19, vC41, vC57);
        faces.add(tri16);
        Triad<Apfloat> tri16_norm = normalTriple(vC19, vC41, vC57, true);
        face_norms.add(tri16_norm);
        Triad<Tuple<Apfloat>> tri17 = new Triad<>(vC19, vC57, vC33);
        faces.add(tri17);
        Triad<Apfloat> tri17_norm = normalTriple(vC19, vC57, vC33, true);
        face_norms.add(tri17_norm);
        Triad<Tuple<Apfloat>> tri18 = new Triad<>(vC19, vC33, vC9);
        faces.add(tri18);
        Triad<Apfloat> tri18_norm = normalTriple(vC19, vC33, vC9, true);
        face_norms.add(tri18_norm);
        Triad<Tuple<Apfloat>> tri19 = new Triad<>(vC19, vC9, vC1);
        faces.add(tri19);
        Triad<Apfloat> tri19_norm = normalTriple(vC19, vC9, vC1, true);
        face_norms.add(tri19_norm);
        Triad<Tuple<Apfloat>> tri20 = new Triad<>(vC20, vC0, vC6);
        faces.add(tri20);
        Triad<Apfloat> tri20_norm = normalTriple(vC20, vC0, vC6, true);
        face_norms.add(tri20_norm);
        Triad<Tuple<Apfloat>> tri21 = new Triad<>(vC20, vC6, vC34);
        faces.add(tri21);
        Triad<Apfloat> tri21_norm = normalTriple(vC20, vC6, vC34, true);
        face_norms.add(tri21_norm);
        Triad<Tuple<Apfloat>> tri22 = new Triad<>(vC20, vC34, vC58);
        faces.add(tri22);
        Triad<Apfloat> tri22_norm = normalTriple(vC20, vC34, vC58, true);
        face_norms.add(tri22_norm);
        Triad<Tuple<Apfloat>> tri23 = new Triad<>(vC20, vC58, vC42);
        faces.add(tri23);
        Triad<Apfloat> tri23_norm = normalTriple(vC20, vC58, vC42, true);
        face_norms.add(tri23_norm);
        Triad<Tuple<Apfloat>> tri24 = new Triad<>(vC20, vC42, vC12);
        faces.add(tri24);
        Triad<Apfloat> tri24_norm = normalTriple(vC20, vC42, vC12, true);
        face_norms.add(tri24_norm);
        Triad<Tuple<Apfloat>> tri25 = new Triad<>(vC20, vC12, vC44);
        faces.add(tri25);
        Triad<Apfloat> tri25_norm = normalTriple(vC20, vC12, vC44, true);
        face_norms.add(tri25_norm);
        Triad<Tuple<Apfloat>> tri26 = new Triad<>(vC20, vC44, vC60);
        faces.add(tri26);
        Triad<Apfloat> tri26_norm = normalTriple(vC20, vC44, vC60, true);
        face_norms.add(tri26_norm);
        Triad<Tuple<Apfloat>> tri27 = new Triad<>(vC20, vC60, vC36);
        faces.add(tri27);
        Triad<Apfloat> tri27_norm = normalTriple(vC20, vC60, vC36, true);
        face_norms.add(tri27_norm);
        Triad<Tuple<Apfloat>> tri28 = new Triad<>(vC20, vC36, vC8);
        faces.add(tri28);
        Triad<Apfloat> tri28_norm = normalTriple(vC20, vC36, vC8, true);
        face_norms.add(tri28_norm);
        Triad<Tuple<Apfloat>> tri29 = new Triad<>(vC20, vC8, vC0);
        faces.add(tri29);
        Triad<Apfloat> tri29_norm = normalTriple(vC20, vC8, vC0, true);
        face_norms.add(tri29_norm);
        Triad<Tuple<Apfloat>> tri30 = new Triad<>(vC21, vC1, vC9);
        faces.add(tri30);
        Triad<Apfloat> tri30_norm = normalTriple(vC21, vC1, vC9, true);
        face_norms.add(tri30_norm);
        Triad<Tuple<Apfloat>> tri31 = new Triad<>(vC21, vC9, vC37);
        faces.add(tri31);
        Triad<Apfloat> tri31_norm = normalTriple(vC21, vC9, vC37, true);
        face_norms.add(tri31_norm);
        Triad<Tuple<Apfloat>> tri32 = new Triad<>(vC21, vC37, vC61);
        faces.add(tri32);
        Triad<Apfloat> tri32_norm = normalTriple(vC21, vC37, vC61, true);
        face_norms.add(tri32_norm);
        Triad<Tuple<Apfloat>> tri33 = new Triad<>(vC21, vC61, vC45);
        faces.add(tri33);
        Triad<Apfloat> tri33_norm = normalTriple(vC21, vC61, vC45, true);
        face_norms.add(tri33_norm);
        Triad<Tuple<Apfloat>> tri34 = new Triad<>(vC21, vC45, vC13);
        faces.add(tri34);
        Triad<Apfloat> tri34_norm = normalTriple(vC21, vC45, vC13, true);
        face_norms.add(tri34_norm);
        Triad<Tuple<Apfloat>> tri35 = new Triad<>(vC21, vC13, vC43);
        faces.add(tri35);
        Triad<Apfloat> tri35_norm = normalTriple(vC21, vC13, vC43, true);
        face_norms.add(tri35_norm);
        Triad<Tuple<Apfloat>> tri36 = new Triad<>(vC21, vC43, vC59);
        faces.add(tri36);
        Triad<Apfloat> tri36_norm = normalTriple(vC21, vC43, vC59, true);
        face_norms.add(tri36_norm);
        Triad<Tuple<Apfloat>> tri37 = new Triad<>(vC21, vC59, vC35);
        faces.add(tri37);
        Triad<Apfloat> tri37_norm = normalTriple(vC21, vC59, vC35, true);
        face_norms.add(tri37_norm);
        Triad<Tuple<Apfloat>> tri38 = new Triad<>(vC21, vC35, vC7);
        faces.add(tri38);
        Triad<Apfloat> tri38_norm = normalTriple(vC21, vC35, vC7, true);
        face_norms.add(tri38_norm);
        Triad<Tuple<Apfloat>> tri39 = new Triad<>(vC21, vC7, vC1);
        faces.add(tri39);
        Triad<Apfloat> tri39_norm = normalTriple(vC21, vC7, vC1, true);
        face_norms.add(tri39_norm);
        Triad<Tuple<Apfloat>> tri40 = new Triad<>(vC22, vC2, vC11);
        faces.add(tri40);
        Triad<Apfloat> tri40_norm = normalTriple(vC22, vC2, vC11, true);
        face_norms.add(tri40_norm);
        Triad<Tuple<Apfloat>> tri41 = new Triad<>(vC22, vC11, vC39);
        faces.add(tri41);
        Triad<Apfloat> tri41_norm = normalTriple(vC22, vC11, vC39, true);
        face_norms.add(tri41_norm);
        Triad<Tuple<Apfloat>> tri42 = new Triad<>(vC22, vC39, vC55);
        faces.add(tri42);
        Triad<Apfloat> tri42_norm = normalTriple(vC22, vC39, vC55, true);
        face_norms.add(tri42_norm);
        Triad<Tuple<Apfloat>> tri43 = new Triad<>(vC22, vC55, vC47);
        faces.add(tri43);
        Triad<Apfloat> tri43_norm = normalTriple(vC22, vC55, vC47, true);
        face_norms.add(tri43_norm);
        Triad<Tuple<Apfloat>> tri44 = new Triad<>(vC22, vC47, vC14);
        faces.add(tri44);
        Triad<Apfloat> tri44_norm = normalTriple(vC22, vC47, vC14, true);
        face_norms.add(tri44_norm);
        Triad<Tuple<Apfloat>> tri45 = new Triad<>(vC22, vC14, vC46);
        faces.add(tri45);
        Triad<Apfloat> tri45_norm = normalTriple(vC22, vC14, vC46, true);
        face_norms.add(tri45_norm);
        Triad<Tuple<Apfloat>> tri46 = new Triad<>(vC22, vC46, vC54);
        faces.add(tri46);
        Triad<Apfloat> tri46_norm = normalTriple(vC22, vC46, vC54, true);
        face_norms.add(tri46_norm);
        Triad<Tuple<Apfloat>> tri47 = new Triad<>(vC22, vC54, vC38);
        faces.add(tri47);
        Triad<Apfloat> tri47_norm = normalTriple(vC22, vC54, vC38, true);
        face_norms.add(tri47_norm);
        Triad<Tuple<Apfloat>> tri48 = new Triad<>(vC22, vC38, vC10);
        faces.add(tri48);
        Triad<Apfloat> tri48_norm = normalTriple(vC22, vC38, vC10, true);
        face_norms.add(tri48_norm);
        Triad<Tuple<Apfloat>> tri49 = new Triad<>(vC22, vC10, vC2);
        faces.add(tri49);
        Triad<Apfloat> tri49_norm = normalTriple(vC22, vC10, vC2, true);
        face_norms.add(tri49_norm);
        Triad<Tuple<Apfloat>> tri50 = new Triad<>(vC23, vC2, vC10);
        faces.add(tri50);
        Triad<Apfloat> tri50_norm = normalTriple(vC23, vC2, vC10, true);
        face_norms.add(tri50_norm);
        Triad<Tuple<Apfloat>> tri51 = new Triad<>(vC23, vC10, vC40);
        faces.add(tri51);
        Triad<Apfloat> tri51_norm = normalTriple(vC23, vC10, vC40, true);
        face_norms.add(tri51_norm);
        Triad<Tuple<Apfloat>> tri52 = new Triad<>(vC23, vC40, vC56);
        faces.add(tri52);
        Triad<Apfloat> tri52_norm = normalTriple(vC23, vC40, vC56, true);
        face_norms.add(tri52_norm);
        Triad<Tuple<Apfloat>> tri53 = new Triad<>(vC23, vC56, vC48);
        faces.add(tri53);
        Triad<Apfloat> tri53_norm = normalTriple(vC23, vC56, vC48, true);
        face_norms.add(tri53_norm);
        Triad<Tuple<Apfloat>> tri54 = new Triad<>(vC23, vC48, vC15);
        faces.add(tri54);
        Triad<Apfloat> tri54_norm = normalTriple(vC23, vC48, vC15, true);
        face_norms.add(tri54_norm);
        Triad<Tuple<Apfloat>> tri55 = new Triad<>(vC23, vC15, vC49);
        faces.add(tri55);
        Triad<Apfloat> tri55_norm = normalTriple(vC23, vC15, vC49, true);
        face_norms.add(tri55_norm);
        Triad<Tuple<Apfloat>> tri56 = new Triad<>(vC23, vC49, vC57);
        faces.add(tri56);
        Triad<Apfloat> tri56_norm = normalTriple(vC23, vC49, vC57, true);
        face_norms.add(tri56_norm);
        Triad<Tuple<Apfloat>> tri57 = new Triad<>(vC23, vC57, vC41);
        faces.add(tri57);
        Triad<Apfloat> tri57_norm = normalTriple(vC23, vC57, vC41, true);
        face_norms.add(tri57_norm);
        Triad<Tuple<Apfloat>> tri58 = new Triad<>(vC23, vC41, vC11);
        faces.add(tri58);
        Triad<Apfloat> tri58_norm = normalTriple(vC23, vC41, vC11, true);
        face_norms.add(tri58_norm);
        Triad<Tuple<Apfloat>> tri59 = new Triad<>(vC23, vC11, vC2);
        faces.add(tri59);
        Triad<Apfloat> tri59_norm = normalTriple(vC23, vC11, vC2, true);
        face_norms.add(tri59_norm);
        Triad<Tuple<Apfloat>> tri60 = new Triad<>(vC24, vC3, vC12);
        faces.add(tri60);
        Triad<Apfloat> tri60_norm = normalTriple(vC24, vC3, vC12, true);
        face_norms.add(tri60_norm);
        Triad<Tuple<Apfloat>> tri61 = new Triad<>(vC24, vC12, vC42);
        faces.add(tri61);
        Triad<Apfloat> tri61_norm = normalTriple(vC24, vC12, vC42, true);
        face_norms.add(tri61_norm);
        Triad<Tuple<Apfloat>> tri62 = new Triad<>(vC24, vC42, vC58);
        faces.add(tri62);
        Triad<Apfloat> tri62_norm = normalTriple(vC24, vC42, vC58, true);
        face_norms.add(tri62_norm);
        Triad<Tuple<Apfloat>> tri63 = new Triad<>(vC24, vC58, vC50);
        faces.add(tri63);
        Triad<Apfloat> tri63_norm = normalTriple(vC24, vC58, vC50, true);
        face_norms.add(tri63_norm);
        Triad<Tuple<Apfloat>> tri64 = new Triad<>(vC24, vC50, vC16);
        faces.add(tri64);
        Triad<Apfloat> tri64_norm = normalTriple(vC24, vC50, vC16, true);
        face_norms.add(tri64_norm);
        Triad<Tuple<Apfloat>> tri65 = new Triad<>(vC24, vC16, vC51);
        faces.add(tri65);
        Triad<Apfloat> tri65_norm = normalTriple(vC24, vC16, vC51, true);
        face_norms.add(tri65_norm);
        Triad<Tuple<Apfloat>> tri66 = new Triad<>(vC24, vC51, vC59);
        faces.add(tri66);
        Triad<Apfloat> tri66_norm = normalTriple(vC24, vC51, vC59, true);
        face_norms.add(tri66_norm);
        Triad<Tuple<Apfloat>> tri67 = new Triad<>(vC24, vC59, vC43);
        faces.add(tri67);
        Triad<Apfloat> tri67_norm = normalTriple(vC24, vC59, vC43, true);
        face_norms.add(tri67_norm);
        Triad<Tuple<Apfloat>> tri68 = new Triad<>(vC24, vC43, vC13);
        faces.add(tri68);
        Triad<Apfloat> tri68_norm = normalTriple(vC24, vC43, vC13, true);
        face_norms.add(tri68_norm);
        Triad<Tuple<Apfloat>> tri69 = new Triad<>(vC24, vC13, vC3);
        faces.add(tri69);
        Triad<Apfloat> tri69_norm = normalTriple(vC24, vC13, vC3, true);
        face_norms.add(tri69_norm);
        Triad<Tuple<Apfloat>> tri70 = new Triad<>(vC25, vC3, vC13);
        faces.add(tri70);
        Triad<Apfloat> tri70_norm = normalTriple(vC25, vC3, vC13, true);
        face_norms.add(tri70_norm);
        Triad<Tuple<Apfloat>> tri71 = new Triad<>(vC25, vC13, vC45);
        faces.add(tri71);
        Triad<Apfloat> tri71_norm = normalTriple(vC25, vC13, vC45, true);
        face_norms.add(tri71_norm);
        Triad<Tuple<Apfloat>> tri72 = new Triad<>(vC25, vC45, vC61);
        faces.add(tri72);
        Triad<Apfloat> tri72_norm = normalTriple(vC25, vC45, vC61, true);
        face_norms.add(tri72_norm);
        Triad<Tuple<Apfloat>> tri73 = new Triad<>(vC25, vC61, vC53);
        faces.add(tri73);
        Triad<Apfloat> tri73_norm = normalTriple(vC25, vC61, vC53, true);
        face_norms.add(tri73_norm);
        Triad<Tuple<Apfloat>> tri74 = new Triad<>(vC25, vC53, vC17);
        faces.add(tri74);
        Triad<Apfloat> tri74_norm = normalTriple(vC25, vC53, vC17, true);
        face_norms.add(tri74_norm);
        Triad<Tuple<Apfloat>> tri75 = new Triad<>(vC25, vC17, vC52);
        faces.add(tri75);
        Triad<Apfloat> tri75_norm = normalTriple(vC25, vC17, vC52, true);
        face_norms.add(tri75_norm);
        Triad<Tuple<Apfloat>> tri76 = new Triad<>(vC25, vC52, vC60);
        faces.add(tri76);
        Triad<Apfloat> tri76_norm = normalTriple(vC25, vC52, vC60, true);
        face_norms.add(tri76_norm);
        Triad<Tuple<Apfloat>> tri77 = new Triad<>(vC25, vC60, vC44);
        faces.add(tri77);
        Triad<Apfloat> tri77_norm = normalTriple(vC25, vC60, vC44, true);
        face_norms.add(tri77_norm);
        Triad<Tuple<Apfloat>> tri78 = new Triad<>(vC25, vC44, vC12);
        faces.add(tri78);
        Triad<Apfloat> tri78_norm = normalTriple(vC25, vC44, vC12, true);
        face_norms.add(tri78_norm);
        Triad<Tuple<Apfloat>> tri79 = new Triad<>(vC25, vC12, vC3);
        faces.add(tri79);
        Triad<Apfloat> tri79_norm = normalTriple(vC25, vC12, vC3, true);
        face_norms.add(tri79_norm);
        Triad<Tuple<Apfloat>> tri80 = new Triad<>(vC26, vC4, vC16);
        faces.add(tri80);
        Triad<Apfloat> tri80_norm = normalTriple(vC26, vC4, vC16, true);
        face_norms.add(tri80_norm);
        Triad<Tuple<Apfloat>> tri81 = new Triad<>(vC26, vC16, vC50);
        faces.add(tri81);
        Triad<Apfloat> tri81_norm = normalTriple(vC26, vC16, vC50, true);
        face_norms.add(tri81_norm);
        Triad<Tuple<Apfloat>> tri82 = new Triad<>(vC26, vC50, vC58);
        faces.add(tri82);
        Triad<Apfloat> tri82_norm = normalTriple(vC26, vC50, vC58, true);
        face_norms.add(tri82_norm);
        Triad<Tuple<Apfloat>> tri83 = new Triad<>(vC26, vC58, vC34);
        faces.add(tri83);
        Triad<Apfloat> tri83_norm = normalTriple(vC26, vC58, vC34, true);
        face_norms.add(tri83_norm);
        Triad<Tuple<Apfloat>> tri84 = new Triad<>(vC26, vC34, vC6);
        faces.add(tri84);
        Triad<Apfloat> tri84_norm = normalTriple(vC26, vC34, vC6, true);
        face_norms.add(tri84_norm);
        Triad<Tuple<Apfloat>> tri85 = new Triad<>(vC26, vC6, vC30);
        faces.add(tri85);
        Triad<Apfloat> tri85_norm = normalTriple(vC26, vC6, vC30, true);
        face_norms.add(tri85_norm);
        Triad<Tuple<Apfloat>> tri86 = new Triad<>(vC26, vC30, vC54);
        faces.add(tri86);
        Triad<Apfloat> tri86_norm = normalTriple(vC26, vC30, vC54, true);
        face_norms.add(tri86_norm);
        Triad<Tuple<Apfloat>> tri87 = new Triad<>(vC26, vC54, vC46);
        faces.add(tri87);
        Triad<Apfloat> tri87_norm = normalTriple(vC26, vC54, vC46, true);
        face_norms.add(tri87_norm);
        Triad<Tuple<Apfloat>> tri88 = new Triad<>(vC26, vC46, vC14);
        faces.add(tri88);
        Triad<Apfloat> tri88_norm = normalTriple(vC26, vC46, vC14, true);
        face_norms.add(tri88_norm);
        Triad<Tuple<Apfloat>> tri89 = new Triad<>(vC26, vC14, vC4);
        faces.add(tri89);
        Triad<Apfloat> tri89_norm = normalTriple(vC26, vC14, vC4, true);
        face_norms.add(tri89_norm);
        Triad<Tuple<Apfloat>> tri90 = new Triad<>(vC27, vC4, vC14);
        faces.add(tri90);
        Triad<Apfloat> tri90_norm = normalTriple(vC27, vC4, vC14, true);
        face_norms.add(tri90_norm);
        Triad<Tuple<Apfloat>> tri91 = new Triad<>(vC27, vC14, vC47);
        faces.add(tri91);
        Triad<Apfloat> tri91_norm = normalTriple(vC27, vC14, vC47, true);
        face_norms.add(tri91_norm);
        Triad<Tuple<Apfloat>> tri92 = new Triad<>(vC27, vC47, vC55);
        faces.add(tri92);
        Triad<Apfloat> tri92_norm = normalTriple(vC27, vC47, vC55, true);
        face_norms.add(tri92_norm);
        Triad<Tuple<Apfloat>> tri93 = new Triad<>(vC27, vC55, vC31);
        faces.add(tri93);
        Triad<Apfloat> tri93_norm = normalTriple(vC27, vC55, vC31, true);
        face_norms.add(tri93_norm);
        Triad<Tuple<Apfloat>> tri94 = new Triad<>(vC27, vC31, vC7);
        faces.add(tri94);
        Triad<Apfloat> tri94_norm = normalTriple(vC27, vC31, vC7, true);
        face_norms.add(tri94_norm);
        Triad<Tuple<Apfloat>> tri95 = new Triad<>(vC27, vC7, vC35);
        faces.add(tri95);
        Triad<Apfloat> tri95_norm = normalTriple(vC27, vC7, vC35, true);
        face_norms.add(tri95_norm);
        Triad<Tuple<Apfloat>> tri96 = new Triad<>(vC27, vC35, vC59);
        faces.add(tri96);
        Triad<Apfloat> tri96_norm = normalTriple(vC27, vC35, vC59, true);
        face_norms.add(tri96_norm);
        Triad<Tuple<Apfloat>> tri97 = new Triad<>(vC27, vC59, vC51);
        faces.add(tri97);
        Triad<Apfloat> tri97_norm = normalTriple(vC27, vC59, vC51, true);
        face_norms.add(tri97_norm);
        Triad<Tuple<Apfloat>> tri98 = new Triad<>(vC27, vC51, vC16);
        faces.add(tri98);
        Triad<Apfloat> tri98_norm = normalTriple(vC27, vC51, vC16, true);
        face_norms.add(tri98_norm);
        Triad<Tuple<Apfloat>> tri99 = new Triad<>(vC27, vC16, vC4);
        faces.add(tri99);
        Triad<Apfloat> tri99_norm = normalTriple(vC27, vC16, vC4, true);
        face_norms.add(tri99_norm);
        Triad<Tuple<Apfloat>> tri100 = new Triad<>(vC28, vC5, vC15);
        faces.add(tri100);
        Triad<Apfloat> tri100_norm = normalTriple(vC28, vC5, vC15, true);
        face_norms.add(tri100_norm);
        Triad<Tuple<Apfloat>> tri101 = new Triad<>(vC28, vC15, vC48);
        faces.add(tri101);
        Triad<Apfloat> tri101_norm = normalTriple(vC28, vC15, vC48, true);
        face_norms.add(tri101_norm);
        Triad<Tuple<Apfloat>> tri102 = new Triad<>(vC28, vC48, vC56);
        faces.add(tri102);
        Triad<Apfloat> tri102_norm = normalTriple(vC28, vC48, vC56, true);
        face_norms.add(tri102_norm);
        Triad<Tuple<Apfloat>> tri103 = new Triad<>(vC28, vC56, vC32);
        faces.add(tri103);
        Triad<Apfloat> tri103_norm = normalTriple(vC28, vC56, vC32, true);
        face_norms.add(tri103_norm);
        Triad<Tuple<Apfloat>> tri104 = new Triad<>(vC28, vC32, vC8);
        faces.add(tri104);
        Triad<Apfloat> tri104_norm = normalTriple(vC28, vC32, vC8, true);
        face_norms.add(tri104_norm);
        Triad<Tuple<Apfloat>> tri105 = new Triad<>(vC28, vC8, vC36);
        faces.add(tri105);
        Triad<Apfloat> tri105_norm = normalTriple(vC28, vC8, vC36, true);
        face_norms.add(tri105_norm);
        Triad<Tuple<Apfloat>> tri106 = new Triad<>(vC28, vC36, vC60);
        faces.add(tri106);
        Triad<Apfloat> tri106_norm = normalTriple(vC28, vC36, vC60, true);
        face_norms.add(tri106_norm);
        Triad<Tuple<Apfloat>> tri107 = new Triad<>(vC28, vC60, vC52);
        faces.add(tri107);
        Triad<Apfloat> tri107_norm = normalTriple(vC28, vC60, vC52, true);
        face_norms.add(tri107_norm);
        Triad<Tuple<Apfloat>> tri108 = new Triad<>(vC28, vC52, vC17);
        faces.add(tri108);
        Triad<Apfloat> tri108_norm = normalTriple(vC28, vC52, vC17, true);
        face_norms.add(tri108_norm);
        Triad<Tuple<Apfloat>> tri109 = new Triad<>(vC28, vC17, vC5);
        faces.add(tri109);
        Triad<Apfloat> tri109_norm = normalTriple(vC28, vC17, vC5, true);
        face_norms.add(tri109_norm);
        Triad<Tuple<Apfloat>> tri110 = new Triad<>(vC29, vC5, vC17);
        faces.add(tri110);
        Triad<Apfloat> tri110_norm = normalTriple(vC29, vC5, vC17, true);
        face_norms.add(tri110_norm);
        Triad<Tuple<Apfloat>> tri111 = new Triad<>(vC29, vC17, vC53);
        faces.add(tri111);
        Triad<Apfloat> tri111_norm = normalTriple(vC29, vC17, vC53, true);
        face_norms.add(tri111_norm);
        Triad<Tuple<Apfloat>> tri112 = new Triad<>(vC29, vC53, vC61);
        faces.add(tri112);
        Triad<Apfloat> tri112_norm = normalTriple(vC29, vC53, vC61, true);
        face_norms.add(tri112_norm);
        Triad<Tuple<Apfloat>> tri113 = new Triad<>(vC29, vC61, vC37);
        faces.add(tri113);
        Triad<Apfloat> tri113_norm = normalTriple(vC29, vC61, vC37, true);
        face_norms.add(tri113_norm);
        Triad<Tuple<Apfloat>> tri114 = new Triad<>(vC29, vC37, vC9);
        faces.add(tri114);
        Triad<Apfloat> tri114_norm = normalTriple(vC29, vC37, vC9, true);
        face_norms.add(tri114_norm);
        Triad<Tuple<Apfloat>> tri115 = new Triad<>(vC29, vC9, vC33);
        faces.add(tri115);
        Triad<Apfloat> tri115_norm = normalTriple(vC29, vC9, vC33, true);
        face_norms.add(tri115_norm);
        Triad<Tuple<Apfloat>> tri116 = new Triad<>(vC29, vC33, vC57);
        faces.add(tri116);
        Triad<Apfloat> tri116_norm = normalTriple(vC29, vC33, vC57, true);
        face_norms.add(tri116_norm);
        Triad<Tuple<Apfloat>> tri117 = new Triad<>(vC29, vC57, vC49);
        faces.add(tri117);
        Triad<Apfloat> tri117_norm = normalTriple(vC29, vC57, vC49, true);
        face_norms.add(tri117_norm);
        Triad<Tuple<Apfloat>> tri118 = new Triad<>(vC29, vC49, vC15);
        faces.add(tri118);
        Triad<Apfloat> tri118_norm = normalTriple(vC29, vC49, vC15, true);
        face_norms.add(tri118_norm);
        Triad<Tuple<Apfloat>> tri119 = new Triad<>(vC29, vC15, vC5);
        faces.add(tri119);
        Triad<Apfloat> tri119_norm = normalTriple(vC29, vC15, vC5, true);
        face_norms.add(tri119_norm);
    }

    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        // check if each point lies within bounds of each face
        for (int i = 0; i < faces.size(); i++) {
            // get vertices of vertA
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces.get(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x - vertA_z)
            Triad<Apfloat> m = subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> face_norm = (Triad<Apfloat>) face_norms.get(i);
            Apfloat d = dot_prod(face_norm, m);

            // if d < 0,  point lies behind face (within bounds)
            // if d == 0, point lies on face (within bounds)
            // if d > 0,  point lies in front of face (out of bounds)
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
