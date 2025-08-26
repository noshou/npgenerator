package io.github.noshou.npg.shapes.johnson;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;
import io.github.noshou.npg.shapes.archimedean.*;
import java.util.ArrayList;

/**
 * <p> Represents a <b>Parabigyrate Rhombicosidodecahedron</b>
 * <p>This polyhedron is formed by taking a {@link Rhombicosidodecahedron} and rotating two
 * opposing pentagonal cupola by 36 degrees.
 * <p> It consists of:
 * <ul>
 *     <li>60 vertices</li>
 *     <li>120 edges</li>
 *     <li>62 faces: 20 triangles,30 squares,and 12 pentagons</li>
 * </ul>
 * <p> It is a member of the {@link RhombicoisododecahedronGyrateAbstract}.
 * @see <a href="https://dmccooey.com/polyhedra/ParabigyrateRhombicosidodecahedron.html">
 *      Parabigyrate Rhombicosidodecahedron (David McCooey)</a>
 */
public class RhombicoisododecahedronParaBiGyrate extends RhombicoisododecahedronGyrateAbstract {
    public RhombicoisododecahedronParaBiGyrate (
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
        vB0 = new Triad<>(HALF(),HALF(),C9());
        vB1 = new Triad<>(HALF(),HALF(),NEG_C9());
        vB2 = new Triad<>(HALF(),NEG_HALF(),C9());
        vB3  = new Triad<>(HALF(),NEG_HALF(),NEG_C9());
        vB4  = new Triad<>(NEG_HALF(),HALF(),C9());
        vB5  = new Triad<>(NEG_HALF(),HALF(),NEG_C9());
        vB6  = new Triad<>(NEG_HALF(),NEG_HALF(),C9());
        vB7  = new Triad<>(NEG_HALF(),NEG_HALF(),NEG_C9());
        vB8  = new Triad<>(C9(),HALF(),HALF());
        vB9  = new Triad<>(C9(),HALF(),NEG_HALF());
        vB10 = new Triad<>(C9(),NEG_HALF(),HALF());
        vB11 = new Triad<>(C9(),NEG_HALF(),NEG_HALF());
        vB12 = new Triad<>(NEG_C9(),HALF(),HALF());
        vB13 = new Triad<>(NEG_C9(),HALF(),NEG_HALF());
        vB14 = new Triad<>(NEG_C9(),NEG_HALF(),HALF());
        vB15 = new Triad<>(NEG_C9(),NEG_HALF(),NEG_HALF());
        vB16 = new Triad<>(HALF(),C4(),C6());
        vB17 = new Triad<>(HALF(),C9(),NEG_HALF());
        vB18 = new Triad<>(HALF(),NEG_C9(),HALF());
        vB19 = new Triad<>(HALF(),NEG_C4(),NEG_C6());
        vB20 = new Triad<>(NEG_HALF(),C4(),C6());
        vB21 = new Triad<>(NEG_HALF(),C9(),NEG_HALF());
        vB22 = new Triad<>(NEG_HALF(),NEG_C9(),HALF());
        vB23 = new Triad<>(NEG_HALF(),NEG_C4(),NEG_C6());
        vB24 = new Triad<>(N0(),C10(),C0());
        vB25 = new Triad<>(N0(),C3(),NEG_C7());
        vB26 = new Triad<>(N0(),NEG_C3(),C7());
        vB27 = new Triad<>(N0(),NEG_C10(),NEG_C0());
        vB28 = new Triad<>(C7(),N0(),C3());
        vB29 = new Triad<>(C7(),N0(),NEG_C3());
        vB30 = new Triad<>(NEG_C7(),N0(),C3());
        vB31 = new Triad<>(NEG_C7(),N0(),NEG_C3());
        vB32 = new Triad<>(C3(),C7(),N0());
        vB33 = new Triad<>(C3(),NEG_C7(),N0());
        vB34 = new Triad<>(NEG_C3(),C7(),N0());
        vB35 = new Triad<>(NEG_C3(),NEG_C7(),N0());
        vB36 = new Triad<>(C3(),C1(),C5());
        vB37 = new Triad<>(C3(),C1(),NEG_C5());
        vB38 = new Triad<>(C3(),NEG_C1(),C5());
        vB39 = new Triad<>(C3(),NEG_C1(),NEG_C5());
        vB40 = new Triad<>(NEG_C3(),C1(),C5());
        vB41 = new Triad<>(NEG_C3(),C1(),NEG_C5());
        vB42 = new Triad<>(NEG_C3(),NEG_C1(),C5());
        vB43 = new Triad<>(NEG_C3(),NEG_C1(),NEG_C5());
        vB44 = new Triad<>(C5(),C3(),C1());
        vB45 = new Triad<>(C5(),C3(),NEG_C1());
        vB46 = new Triad<>(C5(),NEG_C3(),C1());
        vB47 = new Triad<>(C5(),NEG_C3(),NEG_C1());
        vB48 = new Triad<>(NEG_C5(),C3(),C1());
        vB49 = new Triad<>(NEG_C5(),C3(),NEG_C1());
        vB50 = new Triad<>(NEG_C5(),NEG_C3(),C1());
        vB51 = new Triad<>(NEG_C5(),NEG_C3(),NEG_C1());
        vB52 = new Triad<>(C1(),C8(),C2());
        vB53 = new Triad<>(C1(),C5(),NEG_C3());
        vB54 = new Triad<>(C1(),NEG_C5(),C3());
        vB55 = new Triad<>(C1(),NEG_C8(),NEG_C2());
        vB56 = new Triad<>(NEG_C1(),C8(),C2());
        vB57 = new Triad<>(NEG_C1(),C5(),NEG_C3());
        vB58 = new Triad<>(NEG_C1(),NEG_C5(),C3());
        vB59 = new Triad<>(NEG_C1(),NEG_C8(),NEG_C2());

        // ==== SCALED VERTICES ====
        scaledVert();

        // === SQUARE FACES ===
        ArrayList<Tetrad<Tuple<Apfloat>>> squares = new ArrayList<>();
        squares.add(new Tetrad<>(vertices.get(32),vertices.get(17),vertices.get(24),vertices.get(52)));
        squares.add(new Tetrad<>(vertices.get(1),vertices.get(25),vertices.get(53),vertices.get(37)));
        squares.add(new Tetrad<>(vertices.get(2),vertices.get(26),vertices.get(54),vertices.get(38)));
        squares.add(new Tetrad<>(vertices.get(33),vertices.get(18),vertices.get(27),vertices.get(55)));
        squares.add(new Tetrad<>(vertices.get(34),vertices.get(56),vertices.get(24),vertices.get(21)));
        squares.add(new Tetrad<>(vertices.get(5),vertices.get(41),vertices.get(57),vertices.get(25)));
        squares.add(new Tetrad<>(vertices.get(6),vertices.get(42),vertices.get(58),vertices.get(26)));
        squares.add(new Tetrad<>(vertices.get(35),vertices.get(59),vertices.get(27),vertices.get(22)));
        squares.add(new Tetrad<>(vertices.get(8),vertices.get(44),vertices.get(36),vertices.get(28)));
        squares.add(new Tetrad<>(vertices.get(9),vertices.get(29),vertices.get(37),vertices.get(45)));
        squares.add(new Tetrad<>(vertices.get(10),vertices.get(28),vertices.get(38),vertices.get(46)));
        squares.add(new Tetrad<>(vertices.get(11),vertices.get(47),vertices.get(39),vertices.get(29)));
        squares.add(new Tetrad<>(vertices.get(12),vertices.get(30),vertices.get(40),vertices.get(48)));
        squares.add(new Tetrad<>(vertices.get(13),vertices.get(49),vertices.get(41),vertices.get(31)));
        squares.add(new Tetrad<>(vertices.get(14),vertices.get(50),vertices.get(42),vertices.get(30)));
        squares.add(new Tetrad<>(vertices.get(15),vertices.get(31),vertices.get(43),vertices.get(51)));
        squares.add(new Tetrad<>(vertices.get(36),vertices.get(44),vertices.get(52),vertices.get(16)));
        squares.add(new Tetrad<>(vertices.get(17),vertices.get(32),vertices.get(45),vertices.get(53)));
        squares.add(new Tetrad<>(vertices.get(18),vertices.get(33),vertices.get(46),vertices.get(54)));
        squares.add(new Tetrad<>(vertices.get(39),vertices.get(47),vertices.get(55),vertices.get(19)));
        squares.add(new Tetrad<>(vertices.get(40),vertices.get(20),vertices.get(56),vertices.get(48)));
        squares.add(new Tetrad<>(vertices.get(21),vertices.get(57),vertices.get(49),vertices.get(34)));
        squares.add(new Tetrad<>(vertices.get(22),vertices.get(58),vertices.get(50),vertices.get(35)));
        squares.add(new Tetrad<>(vertices.get(43),vertices.get(23),vertices.get(59),vertices.get(51)));
        squares.add(new Tetrad<>(vertices.get(0),vertices.get(4),vertices.get(6),vertices.get(2)));
        squares.add(new Tetrad<>(vertices.get(1),vertices.get(3),vertices.get(7),vertices.get(5)));
        squares.add(new Tetrad<>(vertices.get(8),vertices.get(10),vertices.get(11),vertices.get(9)));
        squares.add(new Tetrad<>(vertices.get(12),vertices.get(13),vertices.get(15),vertices.get(14)));
        squares.add(new Tetrad<>(vertices.get(20),vertices.get(4),vertices.get(0),vertices.get(16)));
        squares.add(new Tetrad<>(vertices.get(7),vertices.get(3),vertices.get(19),vertices.get(23)));

        // === PENTAGONAL FACES ===
        ArrayList<Pentad<Tuple<Apfloat>>> pentagons = new ArrayList<>();
        pentagons.add(new Pentad<>(vertices.get(24),vertices.get(56),vertices.get(20),vertices.get(16),vertices.get(52)));
        pentagons.add(new Pentad<>(vertices.get(25),vertices.get(57),vertices.get(21),vertices.get(17),vertices.get(53)));
        pentagons.add(new Pentad<>(vertices.get(26),vertices.get(58),vertices.get(22),vertices.get(18),vertices.get(54)));
        pentagons.add(new Pentad<>(vertices.get(27),vertices.get(59),vertices.get(23),vertices.get(19),vertices.get(55)));
        pentagons.add(new Pentad<>(vertices.get(28),vertices.get(36),vertices.get(0),vertices.get(2),vertices.get(38)));
        pentagons.add(new Pentad<>(vertices.get(29),vertices.get(39),vertices.get(3),vertices.get(1),vertices.get(37)));
        pentagons.add(new Pentad<>(vertices.get(30),vertices.get(42),vertices.get(6),vertices.get(4),vertices.get(40)));
        pentagons.add(new Pentad<>(vertices.get(31),vertices.get(41),vertices.get(5),vertices.get(7),vertices.get(43)));
        pentagons.add(new Pentad<>(vertices.get(32),vertices.get(44),vertices.get(8),vertices.get(9),vertices.get(45)));
        pentagons.add(new Pentad<>(vertices.get(33),vertices.get(47),vertices.get(11),vertices.get(10),vertices.get(46)));
        pentagons.add(new Pentad<>(vertices.get(34),vertices.get(49),vertices.get(13),vertices.get(12),vertices.get(48)));
        pentagons.add(new Pentad<>(vertices.get(35),vertices.get(50),vertices.get(14),vertices.get(15),vertices.get(51)));

        // === TRIANGULAR FACES ===
        ArrayList<Triad<Tuple<Apfloat>>> triangles = new ArrayList<>();
        triangles.add(new Triad<>(vertices.get(24),vertices.get(17),vertices.get(21)));
        triangles.add(new Triad<>(vertices.get(25),vertices.get(1),vertices.get(5)));
        triangles.add(new Triad<>(vertices.get(26),vertices.get(2),vertices.get(6)));
        triangles.add(new Triad<>(vertices.get(27),vertices.get(18),vertices.get(22)));
        triangles.add(new Triad<>(vertices.get(28),vertices.get(10),vertices.get(8)));
        triangles.add(new Triad<>(vertices.get(29),vertices.get(9),vertices.get(11)));
        triangles.add(new Triad<>(vertices.get(30),vertices.get(12),vertices.get(14)));
        triangles.add(new Triad<>(vertices.get(31),vertices.get(15),vertices.get(13)));
        triangles.add(new Triad<>(vertices.get(32),vertices.get(52),vertices.get(44)));
        triangles.add(new Triad<>(vertices.get(33),vertices.get(55),vertices.get(47)));
        triangles.add(new Triad<>(vertices.get(34),vertices.get(48),vertices.get(56)));
        triangles.add(new Triad<>(vertices.get(35),vertices.get(51),vertices.get(59)));
        triangles.add(new Triad<>(vertices.get(36),vertices.get(16),vertices.get(0)));
        triangles.add(new Triad<>(vertices.get(37),vertices.get(53),vertices.get(45)));
        triangles.add(new Triad<>(vertices.get(38),vertices.get(54),vertices.get(46)));
        triangles.add(new Triad<>(vertices.get(39),vertices.get(19),vertices.get(3)));
        triangles.add(new Triad<>(vertices.get(40),vertices.get(4),vertices.get(20)));
        triangles.add(new Triad<>(vertices.get(41),vertices.get(49),vertices.get(57)));
        triangles.add(new Triad<>(vertices.get(42),vertices.get(50),vertices.get(58)));
        triangles.add(new Triad<>(vertices.get(43),vertices.get(7),vertices.get(23)));

        setFaces(pentagons,squares,triangles);
    }
}
