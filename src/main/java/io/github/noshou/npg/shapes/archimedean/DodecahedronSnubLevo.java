package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.catalan.*;
import io.github.noshou.tuple.Pentad;
import io.github.noshou.tuple.Polyad;
import io.github.noshou.tuple.Triad;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;

/**
 * Represents a <b><i>levo</i>-Snub Dodecahedron</b>
 * <p> The right-handed (<i>levo</i>) enantiomorph of the {@link DodecahedronSnub}.
 * <p> It is the dual of the {@link HexecontahedronPentagonalDextro}.
 * @see <a href="https://dmccooey.com/polyhedra/LsnubDodecahedron.html">
 *      <i>levo</i>-Snub Dodecahedron (David McCooey)</a>
 */
public class DodecahedronSnubLevo extends DodecahedronSnub {

    /**
     * Constructs a new shape instance, resolving units and initializing the lattice.
     *
     * @param radius           the radius of the atom as a string (interpreted using {@code radius_type}), non-null
     * @param radius_type      the unit of the radius (must be "pm", "A", "Å", or "nm"), non-null
     * @param lattice_type     the lattice type (currently only {@code Lattice.LatticeType.FCC} is supported), non-null
     * @param precision        the numeric precision for Apfloat operations
     * @param basis            the atom basis for the unit cell (must contain exactly four atoms for FCC), non-null
     * @param lattice_constant the lattice constant (edge length of unit cell) as a string, non-null
     * @param file_name        the base name for any file output operations, non-null
     * @param structure_name   a user-defined name for the structure, non-null
     * @param structure_index  a unique structure ID used for tracking, non-null
     * @throws IllegalArgumentException if the radius unit or lattice type is not supported
     */
    public DodecahedronSnubLevo(
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

        // define vertices
        Triad<Apfloat> vC0  = mult(normalize(new Triad<>(C2(), NEG_C1(), C14())), super.getRadius().toString());
        Triad<Apfloat> vC1  = mult(normalize(new Triad<>(C2(), C1(), NEG_C14())), super.getRadius().toString());
        Triad<Apfloat> vC2  = mult(normalize(new Triad<>(NEG_C2(), C1(), C14())), super.getRadius().toString());
        Triad<Apfloat> vC3  = mult(normalize(new Triad<>(NEG_C2(), NEG_C1(), NEG_C14())), super.getRadius().toString());
        Triad<Apfloat> vC4  = mult(normalize(new Triad<>(C14(), NEG_C2(), C1())), super.getRadius().toString());
        Triad<Apfloat> vC5  = mult(normalize(new Triad<>(C14(), C2(), NEG_C1())), super.getRadius().toString());
        Triad<Apfloat> vC6  = mult(normalize(new Triad<>(NEG_C14(), C2(), C1())), super.getRadius().toString());
        Triad<Apfloat> vC7  = mult(normalize(new Triad<>(NEG_C14(), NEG_C2(), NEG_C1())), super.getRadius().toString());
        Triad<Apfloat> vC8  = mult(normalize(new Triad<>(C1(), NEG_C14(), C2())), super.getRadius().toString());
        Triad<Apfloat> vC9  = mult(normalize(new Triad<>(C1(), C14(), NEG_C2())), super.getRadius().toString());
        Triad<Apfloat> vC10 = mult(normalize(new Triad<>(NEG_C1(), C14(), C2())), super.getRadius().toString());
        Triad<Apfloat> vC11 = mult(normalize(new Triad<>(NEG_C1(), NEG_C14(), NEG_C2())), super.getRadius().toString());
        Triad<Apfloat> vC12 = mult(normalize(new Triad<>(C3(), C4(), C13())), super.getRadius().toString());
        Triad<Apfloat> vC13 = mult(normalize(new Triad<>(C3(), NEG_C4(), NEG_C13())), super.getRadius().toString());
        Triad<Apfloat> vC14 = mult(normalize(new Triad<>(NEG_C3(), NEG_C4(), C13())), super.getRadius().toString());
        Triad<Apfloat> vC15 = mult(normalize(new Triad<>(NEG_C3(), C4(), NEG_C13())), super.getRadius().toString());
        Triad<Apfloat> vC16 = mult(normalize(new Triad<>(C13(), C3(), C4())), super.getRadius().toString());
        Triad<Apfloat> vC17 = mult(normalize(new Triad<>(C13(), NEG_C3(), NEG_C4())), super.getRadius().toString());
        Triad<Apfloat> vC18 = mult(normalize(new Triad<>(NEG_C13(), NEG_C3(), C4())), super.getRadius().toString());
        Triad<Apfloat> vC19 = mult(normalize(new Triad<>(NEG_C13(), C3(), NEG_C4())), super.getRadius().toString());
        Triad<Apfloat> vC20 = mult(normalize(new Triad<>(C4(), C13(), C3())), super.getRadius().toString());
        Triad<Apfloat> vC21 = mult(normalize(new Triad<>(C4(), NEG_C13(), NEG_C3())), super.getRadius().toString());
        Triad<Apfloat> vC22 = mult(normalize(new Triad<>(NEG_C4(), NEG_C13(), C3())), super.getRadius().toString());
        Triad<Apfloat> vC23 = mult(normalize(new Triad<>(NEG_C4(), C13(), NEG_C3())), super.getRadius().toString());
        Triad<Apfloat> vC24 = mult(normalize(new Triad<>(C0(), NEG_C8(), C12())), super.getRadius().toString());
        Triad<Apfloat> vC25 = mult(normalize(new Triad<>(C0(), C8(), NEG_C12())), super.getRadius().toString());
        Triad<Apfloat> vC26 = mult(normalize(new Triad<>(NEG_C0(), C8(), C12())), super.getRadius().toString());
        Triad<Apfloat> vC27 = mult(normalize(new Triad<>(NEG_C0(), NEG_C8(), NEG_C12())), super.getRadius().toString());
        Triad<Apfloat> vC28 = mult(normalize(new Triad<>(C12(), NEG_C0(), C8())), super.getRadius().toString());
        Triad<Apfloat> vC29 = mult(normalize(new Triad<>(C12(), C0(), NEG_C8())), super.getRadius().toString());
        Triad<Apfloat> vC30 = mult(normalize(new Triad<>(NEG_C12(), C0(), C8())), super.getRadius().toString());
        Triad<Apfloat> vC31 = mult(normalize(new Triad<>(NEG_C12(), NEG_C0(), NEG_C8())), super.getRadius().toString());
        Triad<Apfloat> vC32 = mult(normalize(new Triad<>(C8(), NEG_C12(), C0())), super.getRadius().toString());
        Triad<Apfloat> vC33 = mult(normalize(new Triad<>(C8(), C12(), NEG_C0())), super.getRadius().toString());
        Triad<Apfloat> vC34 = mult(normalize(new Triad<>(NEG_C8(), C12(), C0())), super.getRadius().toString());
        Triad<Apfloat> vC35 = mult(normalize(new Triad<>(NEG_C8(), NEG_C12(), NEG_C0())), super.getRadius().toString());
        Triad<Apfloat> vC36 = mult(normalize(new Triad<>(C7(), NEG_C6(), C11())), super.getRadius().toString());
        Triad<Apfloat> vC37 = mult(normalize(new Triad<>(C7(), C6(), NEG_C11())), super.getRadius().toString());
        Triad<Apfloat> vC38 = mult(normalize(new Triad<>(NEG_C7(), C6(), C11())), super.getRadius().toString());
        Triad<Apfloat> vC39 = mult(normalize(new Triad<>(NEG_C7(), NEG_C6(), NEG_C11())), super.getRadius().toString());
        Triad<Apfloat> vC40 = mult(normalize(new Triad<>(C11(), NEG_C7(), C6())), super.getRadius().toString());
        Triad<Apfloat> vC41 = mult(normalize(new Triad<>(C11(), C7(), NEG_C6())), super.getRadius().toString());
        Triad<Apfloat> vC42 = mult(normalize(new Triad<>(NEG_C11(), C7(), C6())), super.getRadius().toString());
        Triad<Apfloat> vC43 = mult(normalize(new Triad<>(NEG_C11(), NEG_C7(), NEG_C6())), super.getRadius().toString());
        Triad<Apfloat> vC44 = mult(normalize(new Triad<>(C6(), NEG_C11(), C7())), super.getRadius().toString());
        Triad<Apfloat> vC45 = mult(normalize(new Triad<>(C6(), C11(), NEG_C7())), super.getRadius().toString());
        Triad<Apfloat> vC46 = mult(normalize(new Triad<>(NEG_C6(), C11(), C7())), super.getRadius().toString());
        Triad<Apfloat> vC47 = mult(normalize(new Triad<>(NEG_C6(), NEG_C11(), NEG_C7())), super.getRadius().toString());
        Triad<Apfloat> vC48 = mult(normalize(new Triad<>(C9(), C5(), C10())), super.getRadius().toString());
        Triad<Apfloat> vC49 = mult(normalize(new Triad<>(C9(), NEG_C5(), NEG_C10())), super.getRadius().toString());
        Triad<Apfloat> vC50 = mult(normalize(new Triad<>(NEG_C9(), NEG_C5(), C10())), super.getRadius().toString());
        Triad<Apfloat> vC51 = mult(normalize(new Triad<>(NEG_C9(), C5(), NEG_C10())), super.getRadius().toString());
        Triad<Apfloat> vC52 = mult(normalize(new Triad<>(C10(), C9(), C5())), super.getRadius().toString());
        Triad<Apfloat> vC53 = mult(normalize(new Triad<>(C10(), NEG_C9(), NEG_C5())), super.getRadius().toString());
        Triad<Apfloat> vC54 = mult(normalize(new Triad<>(NEG_C10(), NEG_C9(), C5())), super.getRadius().toString());
        Triad<Apfloat> vC55 = mult(normalize(new Triad<>(NEG_C10(), C9(), NEG_C5())), super.getRadius().toString());
        Triad<Apfloat> vC56 = mult(normalize(new Triad<>(C5(), C10(), C9())), super.getRadius().toString());
        Triad<Apfloat> vC57 = mult(normalize(new Triad<>(C5(), NEG_C10(), NEG_C9())), super.getRadius().toString());
        Triad<Apfloat> vC58 = mult(normalize(new Triad<>(NEG_C5(), NEG_C10(), C9())), super.getRadius().toString());
        Triad<Apfloat> vC59 = mult(normalize(new Triad<>(NEG_C5(), C10(), NEG_C9())), super.getRadius().toString());

        // ==== PENTAGONAL FACES ====

        faces_pnt.add(new Pentad<>(vC0, vC36, vC28, vC48, vC12));
        face_norms_pnt.add(normalPent(vC0, vC36, vC28, vC48, vC12, true));

        faces_pnt.add(new Pentad<>(vC1, vC37, vC29, vC49, vC13));
        face_norms_pnt.add(normalPent(vC1, vC37, vC29, vC49, vC13, true));

        faces_pnt.add(new Pentad<>(vC2, vC38, vC30, vC50, vC14));
        face_norms_pnt.add(normalPent(vC2, vC38, vC30, vC50, vC14, true));

        faces_pnt.add(new Pentad<>(vC3, vC39, vC31, vC51, vC15));
        face_norms_pnt.add(normalPent(vC3, vC39, vC31, vC51, vC15, true));

        faces_pnt.add(new Pentad<>(vC4, vC40, vC32, vC53, vC17));
        face_norms_pnt.add(normalPent(vC4, vC40, vC32, vC53, vC17, true));

        faces_pnt.add(new Pentad<>(vC5, vC41, vC33, vC52, vC16));
        face_norms_pnt.add(normalPent(vC5, vC41, vC33, vC52, vC16, true));

        faces_pnt.add(new Pentad<>(vC6, vC42, vC34, vC55, vC19));
        face_norms_pnt.add(normalPent(vC6, vC42, vC34, vC55, vC19, true));

        faces_pnt.add(new Pentad<>(vC7, vC43, vC35, vC54, vC18));
        face_norms_pnt.add(normalPent(vC7, vC43, vC35, vC54, vC18, true));

        faces_pnt.add(new Pentad<>(vC8, vC44, vC24, vC58, vC22));
        face_norms_pnt.add(normalPent(vC8, vC44, vC24, vC58, vC22, true));

        faces_pnt.add(new Pentad<>(vC9, vC45, vC25, vC59, vC23));
        face_norms_pnt.add(normalPent(vC9, vC45, vC25, vC59, vC23, true));

        faces_pnt.add(new Pentad<>(vC10, vC46, vC26, vC56, vC20));
        face_norms_pnt.add(normalPent(vC10, vC46, vC26, vC56, vC20, true));

        faces_pnt.add(new Pentad<>(vC11, vC47, vC27, vC57, vC21));
        face_norms_pnt.add(normalPent(vC11, vC47, vC27, vC57, vC21, true));


        // ==== TRIANGULAR FACES ====

        faces_tri.add(new Triad<>(vC0, vC2, vC14));
        face_norms_tri.add(normalTriple(vC0, vC2, vC14, true));
        faces_tri.add(new Triad<>(vC1, vC3, vC15));
        face_norms_tri.add(normalTriple(vC1, vC3, vC15, true));

        faces_tri.add(new Triad<>(vC2, vC0, vC12));
        face_norms_tri.add(normalTriple(vC2, vC0, vC12, true));

        faces_tri.add(new Triad<>(vC3, vC1, vC13));
        face_norms_tri.add(normalTriple(vC3, vC1, vC13, true));

        faces_tri.add(new Triad<>(vC4, vC5, vC16));
        face_norms_tri.add(normalTriple(vC4, vC5, vC16, true));

        faces_tri.add(new Triad<>(vC5, vC4, vC17));
        face_norms_tri.add(normalTriple(vC5, vC4, vC17, true));

        faces_tri.add(new Triad<>(vC6, vC7, vC18));
        face_norms_tri.add(normalTriple(vC6, vC7, vC18, true));

        faces_tri.add(new Triad<>(vC7, vC6, vC19));
        face_norms_tri.add(normalTriple(vC7, vC6, vC19, true));

        faces_tri.add(new Triad<>(vC8, vC11, vC21));
        face_norms_tri.add(normalTriple(vC8, vC11, vC21, true));

        faces_tri.add(new Triad<>(vC9, vC10, vC20));
        face_norms_tri.add(normalTriple(vC9, vC10, vC20, true));

        faces_tri.add(new Triad<>(vC10, vC9, vC23));
        face_norms_tri.add(normalTriple(vC10, vC9, vC23, true));

        faces_tri.add(new Triad<>(vC11, vC8, vC22));
        face_norms_tri.add(normalTriple(vC11, vC8, vC22, true));

        faces_tri.add(new Triad<>(vC12, vC48, vC56));
        face_norms_tri.add(normalTriple(vC12, vC48, vC56, true));

        faces_tri.add(new Triad<>(vC13, vC49, vC57));
        face_norms_tri.add(normalTriple(vC13, vC49, vC57, true));

        faces_tri.add(new Triad<>(vC14, vC50, vC58));
        face_norms_tri.add(normalTriple(vC14, vC50, vC58, true));

        faces_tri.add(new Triad<>(vC15, vC51, vC59));
        face_norms_tri.add(normalTriple(vC15, vC51, vC59, true));

        faces_tri.add(new Triad<>(vC16, vC52, vC48));
        face_norms_tri.add(normalTriple(vC16, vC52, vC48, true));

        faces_tri.add(new Triad<>(vC17, vC53, vC49));
        face_norms_tri.add(normalTriple(vC17, vC53, vC49, true));

        faces_tri.add(new Triad<>(vC18, vC54, vC50));
        face_norms_tri.add(normalTriple(vC18, vC54, vC50, true));

        faces_tri.add(new Triad<>(vC19, vC55, vC51));
        face_norms_tri.add(normalTriple(vC19, vC55, vC51, true));

        faces_tri.add(new Triad<>(vC20, vC56, vC52));
        face_norms_tri.add(normalTriple(vC20, vC56, vC52, true));

        faces_tri.add(new Triad<>(vC21, vC57, vC53));
        face_norms_tri.add(normalTriple(vC21, vC57, vC53, true));

        faces_tri.add(new Triad<>(vC22, vC58, vC54));
        face_norms_tri.add(normalTriple(vC22, vC58, vC54, true));

        faces_tri.add(new Triad<>(vC23, vC59, vC55));
        face_norms_tri.add(normalTriple(vC23, vC59, vC55, true));

        faces_tri.add(new Triad<>(vC24, vC44, vC36));
        face_norms_tri.add(normalTriple(vC24, vC44, vC36, true));

        faces_tri.add(new Triad<>(vC25, vC45, vC37));
        face_norms_tri.add(normalTriple(vC25, vC45, vC37, true));

        faces_tri.add(new Triad<>(vC26, vC46, vC38));
        face_norms_tri.add(normalTriple(vC26, vC46, vC38, true));

        faces_tri.add(new Triad<>(vC27, vC47, vC39));
        face_norms_tri.add(normalTriple(vC27, vC47, vC39, true));

        faces_tri.add(new Triad<>(vC28, vC36, vC40));
        face_norms_tri.add(normalTriple(vC28, vC36, vC40, true));

        faces_tri.add(new Triad<>(vC29, vC37, vC41));
        face_norms_tri.add(normalTriple(vC29, vC37, vC41, true));

        faces_tri.add(new Triad<>(vC30, vC38, vC42));
        face_norms_tri.add(normalTriple(vC30, vC38, vC42, true));

        faces_tri.add(new Triad<>(vC31, vC39, vC43));
        face_norms_tri.add(normalTriple(vC31, vC39, vC43, true));

        faces_tri.add(new Triad<>(vC32, vC40, vC44));
        face_norms_tri.add(normalTriple(vC32, vC40, vC44, true));

        faces_tri.add(new Triad<>(vC33, vC41, vC45));
        face_norms_tri.add(normalTriple(vC33, vC41, vC45, true));

        faces_tri.add(new Triad<>(vC34, vC42, vC46));
        face_norms_tri.add(normalTriple(vC34, vC42, vC46, true));

        faces_tri.add(new Triad<>(vC35, vC43, vC47));
        face_norms_tri.add(normalTriple(vC35, vC43, vC47, true));

        faces_tri.add(new Triad<>(vC36, vC0, vC24));
        face_norms_tri.add(normalTriple(vC36, vC0, vC24, true));

        faces_tri.add(new Triad<>(vC37, vC1, vC25));
        face_norms_tri.add(normalTriple(vC37, vC1, vC25, true));

        faces_tri.add(new Triad<>(vC38, vC2, vC26));
        face_norms_tri.add(normalTriple(vC38, vC2, vC26, true));

        faces_tri.add(new Triad<>(vC39, vC3, vC27));
        face_norms_tri.add(normalTriple(vC39, vC3, vC27, true));

        faces_tri.add(new Triad<>(vC40, vC4, vC28));
        face_norms_tri.add(normalTriple(vC40, vC4, vC28, true));

        faces_tri.add(new Triad<>(vC41, vC5, vC29));
        face_norms_tri.add(normalTriple(vC41, vC5, vC29, true));

        faces_tri.add(new Triad<>(vC42, vC6, vC30));
        face_norms_tri.add(normalTriple(vC42, vC6, vC30, true));

        faces_tri.add(new Triad<>(vC43, vC7, vC31));
        face_norms_tri.add(normalTriple(vC43, vC7, vC31, true));

        faces_tri.add(new Triad<>(vC44, vC8, vC32));
        face_norms_tri.add(normalTriple(vC44, vC8, vC32, true));

        faces_tri.add(new Triad<>(vC45, vC9, vC33));
        face_norms_tri.add(normalTriple(vC45, vC9, vC33, true));

        faces_tri.add(new Triad<>(vC46, vC10, vC34));
        face_norms_tri.add(normalTriple(vC46, vC10, vC34, true));

        faces_tri.add(new Triad<>(vC47, vC11, vC35));
        face_norms_tri.add(normalTriple(vC47, vC11, vC35, true));

        faces_tri.add(new Triad<>(vC48, vC28, vC16));
        face_norms_tri.add(normalTriple(vC48, vC28, vC16, true));

        faces_tri.add(new Triad<>(vC49, vC29, vC17));
        face_norms_tri.add(normalTriple(vC49, vC29, vC17, true));

        faces_tri.add(new Triad<>(vC50, vC30, vC18));
        face_norms_tri.add(normalTriple(vC50, vC30, vC18, true));

        faces_tri.add(new Triad<>(vC51, vC31, vC19));
        face_norms_tri.add(normalTriple(vC51, vC31, vC19, true));

        faces_tri.add(new Triad<>(vC52, vC33, vC20));
        face_norms_tri.add(normalTriple(vC52, vC33, vC20, true));

        faces_tri.add(new Triad<>(vC53, vC32, vC21));
        face_norms_tri.add(normalTriple(vC53, vC32, vC21, true));

        faces_tri.add(new Triad<>(vC54, vC35, vC22));
        face_norms_tri.add(normalTriple(vC54, vC35, vC22, true));

        faces_tri.add(new Triad<>(vC55, vC34, vC23));
        face_norms_tri.add(normalTriple(vC55, vC34, vC23, true));

        faces_tri.add(new Triad<>(vC56, vC26, vC12));
        face_norms_tri.add(normalTriple(vC56, vC26, vC12, true));

        faces_tri.add(new Triad<>(vC57, vC27, vC13));
        face_norms_tri.add(normalTriple(vC57, vC27, vC13, true));

        faces_tri.add(new Triad<>(vC58, vC24, vC14));
        face_norms_tri.add(normalTriple(vC58, vC24, vC14, true));

        faces_tri.add(new Triad<>(vC59, vC25, vC15));
        face_norms_tri.add(normalTriple(vC59, vC25, vC15, true));

        faces_tri.add(new Triad<>(vC24, vC0, vC14));
        face_norms_tri.add(normalTriple(vC24, vC0, vC14, true));

        faces_tri.add(new Triad<>(vC25, vC1, vC15));
        face_norms_tri.add(normalTriple(vC25, vC1, vC15, true));

        faces_tri.add(new Triad<>(vC26, vC2, vC12));
        face_norms_tri.add(normalTriple(vC26, vC2, vC12, true));

        faces_tri.add(new Triad<>(vC27, vC3, vC13));
        face_norms_tri.add(normalTriple(vC27, vC3, vC13, true));

        faces_tri.add(new Triad<>(vC28, vC4, vC16));
        face_norms_tri.add(normalTriple(vC28, vC4, vC16, true));

        faces_tri.add(new Triad<>(vC29, vC5, vC17));
        face_norms_tri.add(normalTriple(vC29, vC5, vC17, true));

        faces_tri.add(new Triad<>(vC30, vC6, vC18));
        face_norms_tri.add(normalTriple(vC30, vC6, vC18, true));

        faces_tri.add(new Triad<>(vC31, vC7, vC19));
        face_norms_tri.add(normalTriple(vC31, vC7, vC19, true));

        faces_tri.add(new Triad<>(vC32, vC8, vC21));
        face_norms_tri.add(normalTriple(vC32, vC8, vC21, true));

        faces_tri.add(new Triad<>(vC33, vC9, vC20));
        face_norms_tri.add(normalTriple(vC33, vC9, vC20, true));

        faces_tri.add(new Triad<>(vC34, vC10, vC23));
        face_norms_tri.add(normalTriple(vC34, vC10, vC23, true));

        faces_tri.add(new Triad<>(vC35, vC11, vC22));
        face_norms_tri.add(normalTriple(vC35, vC11, vC22, true));

        faces_tri.add(new Triad<>(vC36, vC44, vC40));
        face_norms_tri.add(normalTriple(vC36, vC44, vC40, true));

        faces_tri.add(new Triad<>(vC37, vC45, vC41));
        face_norms_tri.add(normalTriple(vC37, vC45, vC41, true));

        faces_tri.add(new Triad<>(vC38, vC46, vC42));
        face_norms_tri.add(normalTriple(vC38, vC46, vC42, true));

        faces_tri.add(new Triad<>(vC39, vC47, vC43));
        face_norms_tri.add(normalTriple(vC39, vC47, vC43, true));

        faces_tri.add(new Triad<>(vC48, vC52, vC56));
        face_norms_tri.add(normalTriple(vC48, vC52, vC56, true));

        faces_tri.add(new Triad<>(vC49, vC53, vC57));
        face_norms_tri.add(normalTriple(vC49, vC53, vC57, true));

        faces_tri.add(new Triad<>(vC50, vC54, vC58));
        face_norms_tri.add(normalTriple(vC50, vC54, vC58, true));

        faces_tri.add(new Triad<>(vC51, vC55, vC59));
        face_norms_tri.add(normalTriple(vC51, vC55, vC59, true));
    }
}

