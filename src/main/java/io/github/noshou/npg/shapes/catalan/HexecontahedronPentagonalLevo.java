package io.github.noshou.npg.shapes.catalan;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;

/**
 * Represents a <b><i>levo</i>-Pentagonal Hexecontahedron</b>
 * <p>
 * The left-handed (<i>levo</i>) enantiomorph of the {@link HexecontahedronPentagonal}
 * @see <a href="https://dmccooey.com/polyhedra/LpentagonalHexecontahedron.html">
 *      <i>levo</i>-Pentagonal Hexecontahedron (David McCooey)</a>
 */
@SuppressWarnings("FieldCanBeLocal")
public class HexecontahedronPentagonalLevo extends HexecontahedronPentagonal {

    public HexecontahedronPentagonalLevo(
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
        ArrayList<Triad<Apfloat>> vBase = new ArrayList<>();
        vBase.add(new Triad<>(NEG_C0(), NEG_C1(), NEG_C19())); // V0
        vBase.add(new Triad<>(NEG_C0(), C1(), C19()));         // V1
        vBase.add(new Triad<>(C0(), C1(), NEG_C19()));         // V2
        vBase.add(new Triad<>(C0(), NEG_C1(), C19()));         // V3
        vBase.add(new Triad<>(NEG_C19(), NEG_C0(), NEG_C1())); // V4
        vBase.add(new Triad<>(NEG_C19(), C0(), C1()));         // V5
        vBase.add(new Triad<>(C19(), C0(), NEG_C1()));         // V6
        vBase.add(new Triad<>(C19(), NEG_C0(), C1()));         // V7
        vBase.add(new Triad<>(NEG_C1(), NEG_C19(), NEG_C0())); // V8
        vBase.add(new Triad<>(NEG_C1(), C19(), C0()));         // V9
        vBase.add(new Triad<>(C1(), C19(), NEG_C0()));         // V10
        vBase.add(new Triad<>(C1(), NEG_C19(), C0()));         // V11
        vBase.add(new Triad<>(N0(), NEG_C5(), NEG_C18()));       // V12
        vBase.add(new Triad<>(N0(), NEG_C5(), C18()));           // V13
        vBase.add(new Triad<>(N0(), C5(), NEG_C18()));           // V14
        vBase.add(new Triad<>(N0(), C5(), C18()));               // V15
        vBase.add(new Triad<>(NEG_C18(), N0(), NEG_C5()));       // V16
        vBase.add(new Triad<>(NEG_C18(), N0(), C5()));           // V17
        vBase.add(new Triad<>(C18(), N0(), NEG_C5()));           // V18
        vBase.add(new Triad<>(C18(), N0(), C5()));               // V19
        vBase.add(new Triad<>(NEG_C5(), NEG_C18(), N0()));       // V20
        vBase.add(new Triad<>(NEG_C5(), C18(), N0()));           // V21
        vBase.add(new Triad<>(C5(), NEG_C18(), N0()));           // V22
        vBase.add(new Triad<>(C5(), C18(), N0()));               // V23
        vBase.add(new Triad<>(NEG_C10(), N0(), NEG_C17()));      // V24
        vBase.add(new Triad<>(NEG_C10(), N0(), C17()));          // V25
        vBase.add(new Triad<>(C10(), N0(), NEG_C17()));          // V26
        vBase.add(new Triad<>(C10(), N0(), C17()));              // V27
        vBase.add(new Triad<>(NEG_C17(), NEG_C10(), N0()));      // V28
        vBase.add(new Triad<>(NEG_C17(), C10(), N0()));          // V29
        vBase.add(new Triad<>(C17(), NEG_C10(), N0()));          // V30
        vBase.add(new Triad<>(C17(), C10(), N0()));              // V31
        vBase.add(new Triad<>(N0(), NEG_C17(), NEG_C10()));      // V32
        vBase.add(new Triad<>(N0(), NEG_C17(), C10()));          // V33
        vBase.add(new Triad<>(N0(), C17(), NEG_C10()));          // V34
        vBase.add(new Triad<>(N0(), C17(), C10()));              // V35
        vBase.add(new Triad<>(NEG_C3(), C6(), NEG_C16()));      // V36
        vBase.add(new Triad<>(NEG_C3(), NEG_C6(), C16()));      // V37
        vBase.add(new Triad<>(C3(), NEG_C6(), NEG_C16()));      // V38
        vBase.add(new Triad<>(C3(), C6(), C16()));              // V39
        vBase.add(new Triad<>(NEG_C16(), C3(), NEG_C6()));      // V40
        vBase.add(new Triad<>(NEG_C16(), NEG_C3(), C6()));      // V41
        vBase.add(new Triad<>(C16(), NEG_C3(), NEG_C6()));      // V42
        vBase.add(new Triad<>(C16(), C3(), C6()));              // V43
        vBase.add(new Triad<>(NEG_C6(), C16(), NEG_C3()));      // V44
        vBase.add(new Triad<>(NEG_C6(), NEG_C16(), C3()));      // V45
        vBase.add(new Triad<>(C6(), NEG_C16(), NEG_C3()));      // V46
        vBase.add(new Triad<>(C6(), C16(), C3()));              // V47
        vBase.add(new Triad<>(NEG_C2(), NEG_C9(), NEG_C15()));  // V48
        vBase.add(new Triad<>(NEG_C2(), C9(), C15()));          // V49
        vBase.add(new Triad<>(C2(), C9(), NEG_C15()));          // V50
        vBase.add(new Triad<>(C2(), NEG_C9(), C15()));          // V51
        vBase.add(new Triad<>(NEG_C15(), NEG_C2(), NEG_C9()));  // V52
        vBase.add(new Triad<>(NEG_C15(), C2(), C9()));          // V53
        vBase.add(new Triad<>(C15(), C2(), NEG_C9()));          // V54
        vBase.add(new Triad<>(C15(), NEG_C2(), C9()));          // V55
        vBase.add(new Triad<>(NEG_C9(), NEG_C15(), NEG_C2()));  // V56
        vBase.add(new Triad<>(NEG_C9(), C15(), C2()));          // V57
        vBase.add(new Triad<>(C9(), C15(), NEG_C2()));          // V58
        vBase.add(new Triad<>(C9(), NEG_C15(), C2()));          // V59
        vBase.add(new Triad<>(NEG_C7(), NEG_C8(), NEG_C14()));  // V60
        vBase.add(new Triad<>(NEG_C7(), C8(), C14()));          // V61
        vBase.add(new Triad<>(C7(), C8(), NEG_C14()));          // V62
        vBase.add(new Triad<>(C7(), NEG_C8(), C14()));          // V63
        vBase.add(new Triad<>(NEG_C14(), NEG_C7(), NEG_C8()));  // V64
        vBase.add(new Triad<>(NEG_C14(), C7(), C8()));          // V65
        vBase.add(new Triad<>(C14(), C7(), NEG_C8()));          // V66
        vBase.add(new Triad<>(C14(), NEG_C7(), C8()));          // V67
        vBase.add(new Triad<>(NEG_C8(), NEG_C14(), NEG_C7()));  // V68
        vBase.add(new Triad<>(NEG_C8(), C14(), C7()));          // V69
        vBase.add(new Triad<>(C8(), C14(), NEG_C7()));          // V70
        vBase.add(new Triad<>(C8(), NEG_C14(), C7()));          // V71
        vBase.add(new Triad<>(NEG_C4(), C12(), NEG_C13()));      // V72
        vBase.add(new Triad<>(NEG_C4(), NEG_C12(), C13()));      // V73
        vBase.add(new Triad<>(C4(), NEG_C12(), NEG_C13()));      // V74
        vBase.add(new Triad<>(C4(), C12(), C13()));              // V75
        vBase.add(new Triad<>(NEG_C13(), C4(), NEG_C12()));      // V76
        vBase.add(new Triad<>(NEG_C13(), NEG_C4(), C12()));      // V77
        vBase.add(new Triad<>(C13(), NEG_C4(), NEG_C12()));      // V78
        vBase.add(new Triad<>(C13(), C4(), C12()));              // V79
        vBase.add(new Triad<>(NEG_C12(), C13(), NEG_C4()));      // V80
        vBase.add(new Triad<>(NEG_C12(), NEG_C13(), C4()));      // V81
        vBase.add(new Triad<>(C12(), NEG_C13(), NEG_C4()));      // V82
        vBase.add(new Triad<>(C12(), C13(), C4()));              // V83
        vBase.add(new Triad<>(NEG_C11(), NEG_C11(), NEG_C11())); // V84
        vBase.add(new Triad<>(NEG_C11(), NEG_C11(), C11()));     // V85
        vBase.add(new Triad<>(NEG_C11(), C11(), NEG_C11()));     // V86
        vBase.add(new Triad<>(NEG_C11(), C11(), C11()));         // V87
        vBase.add(new Triad<>(C11(), NEG_C11(), NEG_C11()));     // V88
        vBase.add(new Triad<>(C11(), NEG_C11(), C11()));         // V89
        vBase.add(new Triad<>(C11(), C11(), NEG_C11()));         // V90
        vBase.add(new Triad<>(C11(), C11(), C11()));             // V91

        // ==== SCALE VERTICES ====
        scaledVert(vBase);

        // ==== PENTAGONAL FACES ====
        ArrayList<Pentad<Tuple<Apfloat>>> pnt_faces = new ArrayList<>();

        pnt_faces.add(new Pentad<>(vertices.get(24), vertices.get(36), vertices.get(14), vertices.get(2), vertices.get(0)));
        pnt_faces.add(new Pentad<>(vertices.get(24), vertices.get(76), vertices.get(86), vertices.get(72), vertices.get(36)));
        pnt_faces.add(new Pentad<>(vertices.get(24), vertices.get(52), vertices.get(16), vertices.get(40), vertices.get(76)));
        pnt_faces.add(new Pentad<>(vertices.get(24), vertices.get(60), vertices.get(84), vertices.get(64), vertices.get(52)));
        pnt_faces.add(new Pentad<>(vertices.get(24), vertices.get(0), vertices.get(12), vertices.get(48), vertices.get(60)));

        pnt_faces.add(new Pentad<>(vertices.get(25), vertices.get(37), vertices.get(13), vertices.get(3), vertices.get(1)));
        pnt_faces.add(new Pentad<>(vertices.get(25), vertices.get(77), vertices.get(85), vertices.get(73), vertices.get(37)));
        pnt_faces.add(new Pentad<>(vertices.get(25), vertices.get(53), vertices.get(17), vertices.get(41), vertices.get(77)));
        pnt_faces.add(new Pentad<>(vertices.get(25), vertices.get(61), vertices.get(87), vertices.get(65), vertices.get(53)));
        pnt_faces.add(new Pentad<>(vertices.get(25), vertices.get(1), vertices.get(15), vertices.get(49), vertices.get(61)));

        pnt_faces.add(new Pentad<>(vertices.get(26), vertices.get(38), vertices.get(12), vertices.get(0), vertices.get(2)));
        pnt_faces.add(new Pentad<>(vertices.get(26), vertices.get(78), vertices.get(88), vertices.get(74), vertices.get(38)));
        pnt_faces.add(new Pentad<>(vertices.get(26), vertices.get(54), vertices.get(18), vertices.get(42), vertices.get(78)));
        pnt_faces.add(new Pentad<>(vertices.get(26), vertices.get(62), vertices.get(90), vertices.get(66), vertices.get(54)));
        pnt_faces.add(new Pentad<>(vertices.get(26), vertices.get(2), vertices.get(14), vertices.get(50), vertices.get(62)));

        pnt_faces.add(new Pentad<>(vertices.get(27), vertices.get(39), vertices.get(15), vertices.get(1), vertices.get(3)));
        pnt_faces.add(new Pentad<>(vertices.get(27), vertices.get(79), vertices.get(91), vertices.get(75), vertices.get(39)));
        pnt_faces.add(new Pentad<>(vertices.get(27), vertices.get(55), vertices.get(19), vertices.get(43), vertices.get(79)));
        pnt_faces.add(new Pentad<>(vertices.get(27), vertices.get(63), vertices.get(89), vertices.get(67), vertices.get(55)));
        pnt_faces.add(new Pentad<>(vertices.get(27), vertices.get(3), vertices.get(13), vertices.get(51), vertices.get(63)));

        pnt_faces.add(new Pentad<>(vertices.get(28), vertices.get(41), vertices.get(17), vertices.get(5), vertices.get(4)));
        pnt_faces.add(new Pentad<>(vertices.get(28), vertices.get(81), vertices.get(85), vertices.get(77), vertices.get(41)));
        pnt_faces.add(new Pentad<>(vertices.get(28), vertices.get(56), vertices.get(20), vertices.get(45), vertices.get(81)));
        pnt_faces.add(new Pentad<>(vertices.get(28), vertices.get(64), vertices.get(84), vertices.get(68), vertices.get(56)));
        pnt_faces.add(new Pentad<>(vertices.get(28), vertices.get(4), vertices.get(16), vertices.get(52), vertices.get(64)));

        pnt_faces.add(new Pentad<>(vertices.get(29), vertices.get(40), vertices.get(16), vertices.get(4), vertices.get(5)));
        pnt_faces.add(new Pentad<>(vertices.get(29), vertices.get(80), vertices.get(86), vertices.get(76), vertices.get(40)));
        pnt_faces.add(new Pentad<>(vertices.get(29), vertices.get(57), vertices.get(21), vertices.get(44), vertices.get(80)));
        pnt_faces.add(new Pentad<>(vertices.get(29), vertices.get(65), vertices.get(87), vertices.get(69), vertices.get(57)));
        pnt_faces.add(new Pentad<>(vertices.get(29), vertices.get(5), vertices.get(17), vertices.get(53), vertices.get(65)));

        pnt_faces.add(new Pentad<>(vertices.get(30), vertices.get(42), vertices.get(18), vertices.get(6), vertices.get(7)));
        pnt_faces.add(new Pentad<>(vertices.get(30), vertices.get(82), vertices.get(88), vertices.get(78), vertices.get(42)));
        pnt_faces.add(new Pentad<>(vertices.get(30), vertices.get(59), vertices.get(22), vertices.get(46), vertices.get(82)));
        pnt_faces.add(new Pentad<>(vertices.get(30), vertices.get(67), vertices.get(89), vertices.get(71), vertices.get(59)));
        pnt_faces.add(new Pentad<>(vertices.get(30), vertices.get(7), vertices.get(19), vertices.get(55), vertices.get(67)));

        pnt_faces.add(new Pentad<>(vertices.get(31), vertices.get(43), vertices.get(19), vertices.get(7), vertices.get(6)));
        pnt_faces.add(new Pentad<>(vertices.get(31), vertices.get(83), vertices.get(91), vertices.get(79), vertices.get(43)));
        pnt_faces.add(new Pentad<>(vertices.get(31), vertices.get(58), vertices.get(23), vertices.get(47), vertices.get(83)));
        pnt_faces.add(new Pentad<>(vertices.get(31), vertices.get(66), vertices.get(90), vertices.get(70), vertices.get(58)));
        pnt_faces.add(new Pentad<>(vertices.get(31), vertices.get(6), vertices.get(18), vertices.get(54), vertices.get(66)));

        pnt_faces.add(new Pentad<>(vertices.get(32), vertices.get(46), vertices.get(22), vertices.get(11), vertices.get(8)));
        pnt_faces.add(new Pentad<>(vertices.get(32), vertices.get(74), vertices.get(88), vertices.get(82), vertices.get(46)));
        pnt_faces.add(new Pentad<>(vertices.get(32), vertices.get(48), vertices.get(12), vertices.get(38), vertices.get(74)));
        pnt_faces.add(new Pentad<>(vertices.get(32), vertices.get(68), vertices.get(84), vertices.get(60), vertices.get(48)));
        pnt_faces.add(new Pentad<>(vertices.get(32), vertices.get(8), vertices.get(20), vertices.get(56), vertices.get(68)));

        pnt_faces.add(new Pentad<>(vertices.get(33), vertices.get(45), vertices.get(20), vertices.get(8), vertices.get(11)));
        pnt_faces.add(new Pentad<>(vertices.get(33), vertices.get(73), vertices.get(85), vertices.get(81), vertices.get(45)));
        pnt_faces.add(new Pentad<>(vertices.get(33), vertices.get(51), vertices.get(13), vertices.get(37), vertices.get(73)));
        pnt_faces.add(new Pentad<>(vertices.get(33), vertices.get(71), vertices.get(89), vertices.get(63), vertices.get(51)));
        pnt_faces.add(new Pentad<>(vertices.get(33), vertices.get(11), vertices.get(22), vertices.get(59), vertices.get(71)));

        pnt_faces.add(new Pentad<>(vertices.get(34), vertices.get(44), vertices.get(21), vertices.get(9), vertices.get(10)));
        pnt_faces.add(new Pentad<>(vertices.get(34), vertices.get(72), vertices.get(86), vertices.get(80), vertices.get(44)));
        pnt_faces.add(new Pentad<>(vertices.get(34), vertices.get(50), vertices.get(14), vertices.get(36), vertices.get(72)));
        pnt_faces.add(new Pentad<>(vertices.get(34), vertices.get(70), vertices.get(90), vertices.get(62), vertices.get(50)));
        pnt_faces.add(new Pentad<>(vertices.get(34), vertices.get(10), vertices.get(23), vertices.get(58), vertices.get(70)));

        pnt_faces.add(new Pentad<>(vertices.get(35), vertices.get(47), vertices.get(23), vertices.get(10), vertices.get(9)));
        pnt_faces.add(new Pentad<>(vertices.get(35), vertices.get(75), vertices.get(91), vertices.get(83), vertices.get(47)));
        pnt_faces.add(new Pentad<>(vertices.get(35), vertices.get(49), vertices.get(15), vertices.get(39), vertices.get(75)));
        pnt_faces.add(new Pentad<>(vertices.get(35), vertices.get(69), vertices.get(87), vertices.get(61), vertices.get(49)));
        pnt_faces.add(new Pentad<>(vertices.get(35), vertices.get(9), vertices.get(21), vertices.get(57), vertices.get(69)));

        pntFaces(pnt_faces);
    }

}
