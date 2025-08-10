package io.github.noshou.npg.shapes.platonic;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.NotNull;
import static io.github.noshou.npg.nputil.VectorMath.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents an <b>Octahedron</b>.
 * <p>An octahedron has 8 triangular faces, 12 edges, and 6 vertices.
 */
@SuppressWarnings("FieldCanBeLocal")
public class Octahedron extends Shape {

    private final Octad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Octad<Tuple<Apfloat>> face_norms_tri;

    public Octahedron(
            @NotNull String radius, @NotNull String radius_type, @NotNull LatticeType lattice_type,
            int precision, @NotNull Polyad<Atom> basis, @NotNull String lattice_constant,
            @NotNull String file_name, @NotNull String structure_name, @NotNull String structure_index
    ) {
        super(radius, radius_type, lattice_type, precision, basis, lattice_constant, file_name, structure_name, structure_index);

        Apfloat N0 = new Apfloat("0", precision);
        Apfloat N1 = new Apfloat("1", precision);
        Apfloat NEG_N1 = new Apfloat("-1", precision);

        // === Canonical vertices ===
        ArrayList<Triad<Apfloat>> v = new ArrayList<>();
        v.add(new Triad<>( N1,  N0, N0));
        v.add(new Triad<>(NEG_N1,  N0, N0));
        v.add(new Triad<>( N0,  N1, N0));
        v.add(new Triad<>( N0, NEG_N1, N0));
        v.add(new Triad<>( N0, N0,  N1));
        v.add(new Triad<>( N0, N0, NEG_N1));

        // Scale (already unit length → radius)
        Apfloat scale = super.getRadius();
        List<Triad<Apfloat>> vScaled = v.stream().map(v_i -> mult(v_i, scale.toString())).toList();

        // === Define 8 triangular faces ===
        ArrayList<Triad<Tuple<Apfloat>>> tr = new ArrayList<>();
        tr.add(new Triad<>(vScaled.get(0), vScaled.get(2), vScaled.get(4)));
        tr.add(new Triad<>(vScaled.get(2), vScaled.get(1), vScaled.get(4)));
        tr.add(new Triad<>(vScaled.get(1), vScaled.get(3), vScaled.get(4)));
        tr.add(new Triad<>(vScaled.get(3), vScaled.get(0), vScaled.get(4)));
        tr.add(new Triad<>(vScaled.get(2), vScaled.get(0), vScaled.get(5)));
        tr.add(new Triad<>(vScaled.get(1), vScaled.get(2), vScaled.get(5)));
        tr.add(new Triad<>(vScaled.get(3), vScaled.get(1), vScaled.get(5)));
        tr.add(new Triad<>(vScaled.get(0), vScaled.get(3), vScaled.get(5)));

        List<Triad<Apfloat>> trn = tr.stream().map(
                f -> normalTriple(
                        (Triad<Apfloat>) f.fetch(0),
                        (Triad<Apfloat>) f.fetch(1),
                        (Triad<Apfloat>) f.fetch(2),
                        true
                )
        ).toList();

        faces_tri = new Octad<>(tr.get(0), tr.get(1), tr.get(2), tr.get(3),
                tr.get(4), tr.get(5), tr.get(6), tr.get(7));
        face_norms_tri = new Octad<>(trn.get(0), trn.get(1), trn.get(2), trn.get(3),
                trn.get(4), trn.get(5), trn.get(6), trn.get(7));
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        Apfloat N0 = new Apfloat("0", super.precision);
        Apfloat THREE = new Apfloat("3", super.precision);

        for (int i = 0; i < faces_tri.fetchSize(); i++) {
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.fetch(i);
            Triad<Apfloat> v0 = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> v1 = (Triad<Apfloat>) face.fetch(1);
            Triad<Apfloat> v2 = (Triad<Apfloat>) face.fetch(2);

            Triad<Apfloat> centroid = new Triad<>(
                    v0.fetch(0).add(v1.fetch(0)).add(v2.fetch(0)).divide(THREE),
                    v0.fetch(1).add(v1.fetch(1)).add(v2.fetch(1)).divide(THREE),
                    v0.fetch(2).add(v1.fetch(2)).add(v2.fetch(2)).divide(THREE)
            );
            Triad<Apfloat> m = subs(point_cart, centroid);
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.fetch(i);
            Apfloat d = dot_prod(norm, m);
            if (d.compareTo(N0) > 0) return false;
        }
        return true;
    }
}