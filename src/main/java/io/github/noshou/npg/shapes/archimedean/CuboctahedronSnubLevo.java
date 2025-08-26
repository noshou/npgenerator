package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;
import java.util.ArrayList;
import io.github.noshou.npg.shapes.catalan.*;

/**
 * Represents a <b><i>levo</i>-Snub Cuboctahedron</b>
 * <p> The left-handed (<i>levo</i>) enantiomorph of the {@link CuboctahedronSnub}.
 * <p> It is the dual of the {@link IcositetrahedronPentagonalLevo}. </p>
 * <p> May be circumscribed ({@link CuboctahedronSnubLevoCanonical}) or
 * circumscribed and inscribed ({@link CuboctahedronSnubLevoBiscribed}).
 */
public abstract class CuboctahedronSnubLevo extends CuboctahedronSnub {

    // list of vertices
    private ArrayList<Triad<Apfloat>> vBase;

    // ==== BASIS VERTICES ====
    // MUST be filled by children
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

    /*
     * Child classes must implement this method in their constructors
     */
    protected void setVerts() {

        vBase.add(vB0);
        vBase.add(vB1);
        vBase.add(vB2);
        vBase.add(vB3);
        vBase.add(vB4);
        vBase.add(vB5);
        vBase.add(vB6);
        vBase.add(vB7);
        vBase.add(vB8);
        vBase.add(vB9);
        vBase.add(vB10);
        vBase.add(vB11);
        vBase.add(vB12);
        vBase.add(vB13);
        vBase.add(vB14);
        vBase.add(vB16);
        vBase.add(vB17);
        vBase.add(vB18);
        vBase.add(vB19);
        vBase.add(vB20);
        vBase.add(vB21);
        vBase.add(vB22);
        vBase.add(vB23);
        
        // ==== SCALED VERTICES ====
        scaledVert(vBase);

        // ==== TRIANGULAR FACES ====
        ArrayList<Triad<Tuple<Apfloat>>> tri_faces = new ArrayList<>();
        tri_faces.add(new Triad<>(vertices.get(0),vertices.get(8),vertices.get(14)));
        tri_faces.add(new Triad<>(vertices.get(1),vertices.get(9),vertices.get(15)));
        tri_faces.add(new Triad<>(vertices.get(2),vertices.get(10),vertices.get(12)));
        tri_faces.add(new Triad<>(vertices.get(3),vertices.get(11),vertices.get(13)));
        tri_faces.add(new Triad<>(vertices.get(4),vertices.get(0),vertices.get(16)));
        tri_faces.add(new Triad<>(vertices.get(5),vertices.get(1),vertices.get(17)));
        tri_faces.add(new Triad<>(vertices.get(6),vertices.get(2),vertices.get(18)));
        tri_faces.add(new Triad<>(vertices.get(7),vertices.get(3),vertices.get(19)));
        tri_faces.add(new Triad<>(vertices.get(8),vertices.get(4),vertices.get(21)));
        tri_faces.add(new Triad<>(vertices.get(9),vertices.get(5),vertices.get(20)));
        tri_faces.add(new Triad<>(vertices.get(10),vertices.get(6),vertices.get(23)));
        tri_faces.add(new Triad<>(vertices.get(11),vertices.get(7),vertices.get(22)));
        tri_faces.add(new Triad<>(vertices.get(12),vertices.get(16),vertices.get(0)));
        tri_faces.add(new Triad<>(vertices.get(13),vertices.get(17),vertices.get(1)));
        tri_faces.add(new Triad<>(vertices.get(14),vertices.get(18),vertices.get(2)));
        tri_faces.add(new Triad<>(vertices.get(15),vertices.get(19),vertices.get(3)));
        tri_faces.add(new Triad<>(vertices.get(16),vertices.get(20),vertices.get(5)));
        tri_faces.add(new Triad<>(vertices.get(17),vertices.get(21),vertices.get(4)));
        tri_faces.add(new Triad<>(vertices.get(18),vertices.get(22),vertices.get(7)));
        tri_faces.add(new Triad<>(vertices.get(19),vertices.get(23),vertices.get(6)));
        tri_faces.add(new Triad<>(vertices.get(20),vertices.get(12),vertices.get(10)));
        tri_faces.add(new Triad<>(vertices.get(21),vertices.get(13),vertices.get(11)));
        tri_faces.add(new Triad<>(vertices.get(22),vertices.get(14),vertices.get(8)));
        tri_faces.add(new Triad<>(vertices.get(23),vertices.get(15),vertices.get(9)));
        tri_faces.add(new Triad<>(vertices.get(8),vertices.get(0),vertices.get(4)));
        tri_faces.add(new Triad<>(vertices.get(9),vertices.get(1),vertices.get(5)));
        tri_faces.add(new Triad<>(vertices.get(10),vertices.get(2),vertices.get(6)));
        tri_faces.add(new Triad<>(vertices.get(11),vertices.get(3),vertices.get(7)));
        tri_faces.add(new Triad<>(vertices.get(12),vertices.get(20),vertices.get(16)));
        tri_faces.add(new Triad<>(vertices.get(13),vertices.get(21),vertices.get(17)));
        tri_faces.add(new Triad<>(vertices.get(14),vertices.get(22),vertices.get(18)));
        tri_faces.add(new Triad<>(vertices.get(15),vertices.get(23),vertices.get(19)));
        triFaces(tri_faces);

        // ==== SQUARE FACES ====
        ArrayList<Tetrad<Tuple<Apfloat>>> sqr_faces = new ArrayList<>();
        sqr_faces.add(new Tetrad<>(vertices.get(2), vertices.get(12), vertices.get(0), vertices.get(14)));
        sqr_faces.add(new Tetrad<>(vertices.get(3), vertices.get(13), vertices.get(1), vertices.get(15)));
        sqr_faces.add(new Tetrad<>(vertices.get(4), vertices.get(16), vertices.get(5), vertices.get(17)));
        sqr_faces.add(new Tetrad<>(vertices.get(7), vertices.get(19), vertices.get(6), vertices.get(18)));
        sqr_faces.add(new Tetrad<>(vertices.get(8), vertices.get(21), vertices.get(11), vertices.get(22)));
        sqr_faces.add(new Tetrad<>(vertices.get(9), vertices.get(20), vertices.get(10), vertices.get(23)));
        sqrFaces(sqr_faces);

    }
    public CuboctahedronSnubLevo(
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
    }
}
