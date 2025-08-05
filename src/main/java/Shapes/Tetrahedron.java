package Shapes;

import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;
import java.util.*;
import Lattice.*;
import Atom.*;

/**
 * Represents a regular tetrahedron inscribed in a sphere of given radius.
 *
 * <p>A tetrahedron has 4 triangular faces, 6 edges, and 4 vertices.</p>
 */
public class Tetrahedron extends Shape {

    /** 4 triangular faces */
    private final Tetrad<Tuple<Tuple<Apfloat>>> faces_tri;

    /** Normals of the 4 triangular faces */
    private final Tetrad<Tuple<Apfloat>> face_norms_tri;

    // === Constants ===
    private final Apfloat ONE = new Apfloat("1", precision);
    private final Apfloat NEG = new Apfloat("-1", precision);
    private final Apfloat ZERO = new Apfloat("0", super.precision);
    private final Apfloat THREE = new Apfloat("3", super.precision);

    /**
     * Constructs a tetrahedron.
     */
    public Tetrahedron(
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
        super(radius, radius_type, lattice_type, precision, basis, lattice_constant, file_name, structure_name, structure_index);


        // === Canonical vertices (tetrahedron centered at origin) ===
        // Using coordinates: (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
        ArrayList<Triad<Apfloat>> v = new ArrayList<>();
        v.add(new Triad<>( ONE,  ONE,  ONE));
        v.add(new Triad<>( ONE, NEG, NEG));
        v.add(new Triad<>(NEG,  ONE, NEG));
        v.add(new Triad<>(NEG, NEG,  ONE));

        // Scale to circumsphere radius
        Apfloat dist = ApfloatMath.sqrt(new Apfloat("3", precision)); // |(1,1,1)|
        Apfloat scale = super.getRadius().divide(dist);
        List<Triad<Apfloat>> vScaled = v.stream().map(v_i -> VectorMath.mult(v_i, scale.toString())).toList();

        // === Define 4 triangular faces ===
        ArrayList<Triad<Tuple<Apfloat>>> tr = new ArrayList<>();
        tr.add(new Triad<>(vScaled.get(0), vScaled.get(1), vScaled.get(2)));
        tr.add(new Triad<>(vScaled.get(0), vScaled.get(3), vScaled.get(1)));
        tr.add(new Triad<>(vScaled.get(0), vScaled.get(2), vScaled.get(3)));
        tr.add(new Triad<>(vScaled.get(1), vScaled.get(3), vScaled.get(2)));

        // Compute normals
        List<Triad<Apfloat>> trn = tr.stream().map(
                tr_i -> VectorMath.normalTriple(
                        (Triad<Apfloat>) tr_i.fetch(0),
                        (Triad<Apfloat>) tr_i.fetch(1),
                        (Triad<Apfloat>) tr_i.fetch(2),
                        true
                )
        ).toList();

        faces_tri = new Tetrad<>(tr.get(0), tr.get(1), tr.get(2), tr.get(3));
        face_norms_tri = new Tetrad<>(trn.get(0), trn.get(1), trn.get(2), trn.get(3));
    }

    @Override
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

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
            Triad<Apfloat> m = VectorMath.subs(point_cart, centroid);
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.fetch(i);
            Apfloat d = VectorMath.dot_prod(norm, m);
            if (d.compareTo(ZERO) > 0) return false;
        }
        return true;
    }
}
