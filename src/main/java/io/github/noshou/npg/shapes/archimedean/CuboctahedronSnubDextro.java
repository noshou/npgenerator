package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;

/**
 * Represents a <b><i>dextro</i>-Snub Cuboctahedron</b>
 * <p>
 * The right-handed (<i>dextro</i>) enantiomorph of the {@link CuboctahedronSnub}
 * @see <a href="https://dmccooey.com/polyhedra/LsnubCube.html">
 *      <i>dextro</i>-Snub Cuboctahedron (David McCooey)</a>
 */
public class CuboctahedronSnubDextro extends CuboctahedronSnub {

    public CuboctahedronSnubDextro(
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
        Triad<Apfloat> vB0 = new Triad<>(C1(), NEG_C0(), C2());
        vBase.add(vB0);
        Triad<Apfloat> vB1 = new Triad<>(C1(), C0(), NEG_C2());
        vBase.add(vB1);
        Triad<Apfloat> vB2 = new Triad<>(NEG_C1(), C0(), C2());
        vBase.add(vB2);
        Triad<Apfloat> vB3 = new Triad<>(NEG_C1(), NEG_C0(), NEG_C2());
        vBase.add(vB3);
        Triad<Apfloat> vB4 = new Triad<>(C2(), NEG_C1(), C0());
        vBase.add(vB4);
        Triad<Apfloat> vB5 = new Triad<>(C2(), C1(), NEG_C0());
        vBase.add(vB5);
        Triad<Apfloat> vB6 = new Triad<>(NEG_C2(), C1(), C0());
        vBase.add(vB6);
        Triad<Apfloat> vB7 = new Triad<>(NEG_C2(), NEG_C1(), NEG_C0());
        vBase.add(vB7);
        Triad<Apfloat> vB8 = new Triad<>(C0(), NEG_C2(), C1());
        vBase.add(vB8);
        Triad<Apfloat> vB9 = new Triad<>(C0(), C2(), NEG_C1());
        vBase.add(vB9);
        Triad<Apfloat> vB10 = new Triad<>(NEG_C0(), C2(), C1());
        vBase.add(vB10);
        Triad<Apfloat> vB11 = new Triad<>(NEG_C0(), NEG_C2(), NEG_C1());
        vBase.add(vB11);
        Triad<Apfloat> vB12 = new Triad<>(C0(), C1(), C2());
        vBase.add(vB12);
        Triad<Apfloat> vB13 = new Triad<>(C0(), NEG_C1(), NEG_C2());
        vBase.add(vB13);
        Triad<Apfloat> vB14 = new Triad<>(NEG_C0(), NEG_C1(), C2());
        vBase.add(vB14);
        Triad<Apfloat> vB15 = new Triad<>(NEG_C0(), C1(), NEG_C2());
        vBase.add(vB15);
        Triad<Apfloat> vB16 = new Triad<>(C2(), C0(), C1());
        vBase.add(vB16);
        Triad<Apfloat> vB17 = new Triad<>(C2(), NEG_C0(), NEG_C1());
        vBase.add(vB17);
        Triad<Apfloat> vB18 = new Triad<>(NEG_C2(), NEG_C0(), C1());
        vBase.add(vB18);
        Triad<Apfloat> vB19 = new Triad<>(NEG_C2(), C0(), NEG_C1());
        vBase.add(vB19);
        Triad<Apfloat> vB20 = new Triad<>(C1(), C2(), C0());
        vBase.add(vB20);
        Triad<Apfloat> vB21 = new Triad<>(C1(), NEG_C2(), NEG_C0());
        vBase.add(vB21);
        Triad<Apfloat> vB22 = new Triad<>(NEG_C1(), NEG_C2(), C0());
        vBase.add(vB22);
        Triad<Apfloat> vB23 = new Triad<>(NEG_C1(), C2(), NEG_C0());
        vBase.add(vB23);

        // ==== SCALE VERTICES ====
        scaledVert(vBase);

        // ==== TRIANGULAR FACES ====
        ArrayList<Triad<Tuple<Apfloat>>> tri_faces = new ArrayList<>();
        tri_faces.add(new Triad<>(vertices.get(0),vertices.get(14),vertices.get(8)));
        tri_faces.add(new Triad<>(vertices.get(1),vertices.get(15),vertices.get(9)));
        tri_faces.add(new Triad<>(vertices.get(2),vertices.get(12),vertices.get(10)));
        tri_faces.add(new Triad<>(vertices.get(3),vertices.get(13),vertices.get(11)));
        tri_faces.add(new Triad<>(vertices.get(4),vertices.get(16),vertices.get(0)));
        tri_faces.add(new Triad<>(vertices.get(5),vertices.get(17),vertices.get(1)));
        tri_faces.add(new Triad<>(vertices.get(6),vertices.get(18),vertices.get(2)));
        tri_faces.add(new Triad<>(vertices.get(7),vertices.get(19),vertices.get(3)));
        tri_faces.add(new Triad<>(vertices.get(8),vertices.get(21),vertices.get(4)));
        tri_faces.add(new Triad<>(vertices.get(9),vertices.get(20),vertices.get(5)));
        tri_faces.add(new Triad<>(vertices.get(10),vertices.get(23),vertices.get(6)));
        tri_faces.add(new Triad<>(vertices.get(11),vertices.get(22),vertices.get(7)));
        tri_faces.add(new Triad<>(vertices.get(12),vertices.get(0),vertices.get(16)));
        tri_faces.add(new Triad<>(vertices.get(13),vertices.get(1),vertices.get(17)));
        tri_faces.add(new Triad<>(vertices.get(14),vertices.get(2),vertices.get(18)));
        tri_faces.add(new Triad<>(vertices.get(15),vertices.get(3),vertices.get(19)));
        tri_faces.add(new Triad<>(vertices.get(16),vertices.get(5),vertices.get(20)));
        tri_faces.add(new Triad<>(vertices.get(17),vertices.get(4),vertices.get(21)));
        tri_faces.add(new Triad<>(vertices.get(18),vertices.get(7),vertices.get(22)));
        tri_faces.add(new Triad<>(vertices.get(19),vertices.get(6),vertices.get(23)));
        tri_faces.add(new Triad<>(vertices.get(20),vertices.get(10),vertices.get(12)));
        tri_faces.add(new Triad<>(vertices.get(21),vertices.get(11),vertices.get(13)));
        tri_faces.add(new Triad<>(vertices.get(22),vertices.get(8),vertices.get(14)));
        tri_faces.add(new Triad<>(vertices.get(23),vertices.get(9),vertices.get(15)));
        tri_faces.add(new Triad<>(vertices.get(8),vertices.get(4),vertices.get(0)));
        tri_faces.add(new Triad<>(vertices.get(9),vertices.get(5),vertices.get(1)));
        tri_faces.add(new Triad<>(vertices.get(10),vertices.get(6),vertices.get(2)));
        tri_faces.add(new Triad<>(vertices.get(11),vertices.get(7),vertices.get(3)));
        tri_faces.add(new Triad<>(vertices.get(12),vertices.get(16),vertices.get(20)));
        tri_faces.add(new Triad<>(vertices.get(13),vertices.get(17),vertices.get(21)));
        tri_faces.add(new Triad<>(vertices.get(14),vertices.get(18),vertices.get(22)));
        tri_faces.add(new Triad<>(vertices.get(15),vertices.get(19),vertices.get(23)));
        triFaces(tri_faces);

        // ==== SQUARE FACES ====
        ArrayList<Tetrad<Tuple<Apfloat>>> sqr_faces = new ArrayList<>();
        sqr_faces.add(new Tetrad<>(vertices.get(2),vertices.get(14),vertices.get(0),vertices.get(12)));
        sqr_faces.add(new Tetrad<>(vertices.get(3),vertices.get(15),vertices.get(1),vertices.get(13)));
        sqr_faces.add(new Tetrad<>(vertices.get(4),vertices.get(17),vertices.get(5),vertices.get(16)));
        sqr_faces.add(new Tetrad<>(vertices.get(7),vertices.get(18),vertices.get(6),vertices.get(19)));
        sqr_faces.add(new Tetrad<>(vertices.get(8),vertices.get(22),vertices.get(11),vertices.get(21)));
        sqr_faces.add(new Tetrad<>(vertices.get(9),vertices.get(23),vertices.get(10),vertices.get(20)));
        sqrFaces(sqr_faces);
    }
}
