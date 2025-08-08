package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.NotNull;

/**
 * {@link:<a href=https://dmccooey.com/polyhedra/RhombicDodecahedron.html>...</a>}
 */
@SuppressWarnings("FieldCanBeLocal")
public class DodecahedronRhombic extends Shape {

    //  12 rhombohedral faces
    private final Dodecad<Tuple<Tuple<Apfloat>>> faces_rho;
    private final Dodecad<Tuple<Apfloat>> face_norms_rho;


    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N4 = new Apfloat("3", super.precision);
    private final Apfloat N3 = new Apfloat("4", super.precision);
    private final Apfloat N8 = new Apfloat("8", super.precision);
    private final Apfloat FRAC_N3mulSQRT2_over_N8 = N3.multiply(ApfloatMath.sqrt(N2)).divide(N8);
    private final Apfloat FRAC_N3mulSQRT2_over_N4 = N3.multiply(ApfloatMath.sqrt(N4)).divide(N8);

    public DodecahedronRhombic(
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
        Triad<Apfloat> vB0 = new Triad<>(
                N0,
                N0,
                FRAC_N3mulSQRT2_over_N4
        );
        Triad<Apfloat> vB1 = new Triad<>(
                N0,
                N0,
                FRAC_N3mulSQRT2_over_N4.multiply(NEG_N1)
        );
        Triad<Apfloat> vB2 = new Triad<>(
                FRAC_N3mulSQRT2_over_N4,
                N0,
                N0
        );
        Triad<Apfloat> vB3 = new Triad<>(
                FRAC_N3mulSQRT2_over_N4.multiply(NEG_N1),
                N0,
                N0
        );
        Triad<Apfloat> vB4 = new Triad<>(
                N0,
                FRAC_N3mulSQRT2_over_N4,
                N0
        );
        Triad<Apfloat> vB5 = new Triad<>(
                N0,
                FRAC_N3mulSQRT2_over_N4.multiply(NEG_N1),
                N0);
        Triad<Apfloat> vB6 = new Triad<>(
                FRAC_N3mulSQRT2_over_N8,
                FRAC_N3mulSQRT2_over_N8,
                FRAC_N3mulSQRT2_over_N8
        );
        Triad<Apfloat> vB7 = new Triad<>(
                FRAC_N3mulSQRT2_over_N8,
                FRAC_N3mulSQRT2_over_N8,
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1)
        );
        Triad<Apfloat> vB8 = new Triad<>(
                FRAC_N3mulSQRT2_over_N8,
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1),
                FRAC_N3mulSQRT2_over_N8
        );
        Triad<Apfloat> vB9 = new Triad<>(
                FRAC_N3mulSQRT2_over_N8,
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1),
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1)
        );
        Triad<Apfloat> vB10 = new Triad<>(
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1),
                FRAC_N3mulSQRT2_over_N8,
                FRAC_N3mulSQRT2_over_N8
        );
        Triad<Apfloat> vB11 = new Triad<>(
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1),
                FRAC_N3mulSQRT2_over_N8,
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1)
        );
        Triad<Apfloat> vB12 = new Triad<>(
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1),
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1),
                FRAC_N3mulSQRT2_over_N8
        );
        Triad<Apfloat> vB13 = new Triad<>(
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1),
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1),
                FRAC_N3mulSQRT2_over_N8.multiply(NEG_N1)
        );

        // ==== SCALED VERTICES ====
        Triad<Apfloat> vC0 = VectorMath.mult(VectorMath.normalize(vB0), super.getRadius().toString());
        Triad<Apfloat> vC1 = VectorMath.mult(VectorMath.normalize(vB1), super.getRadius().toString());
        Triad<Apfloat> vC2 = VectorMath.mult(VectorMath.normalize(vB2), super.getRadius().toString());
        Triad<Apfloat> vC3  = VectorMath.mult(VectorMath.normalize(vB3), super.getRadius().toString());
        Triad<Apfloat> vC4  = VectorMath.mult(VectorMath.normalize(vB4), super.getRadius().toString());
        Triad<Apfloat> vC5  = VectorMath.mult(VectorMath.normalize(vB5), super.getRadius().toString());
        Triad<Apfloat> vC6  = VectorMath.mult(VectorMath.normalize(vB6), super.getRadius().toString());
        Triad<Apfloat> vC7  = VectorMath.mult(VectorMath.normalize(vB7), super.getRadius().toString());
        Triad<Apfloat> vC8  = VectorMath.mult(VectorMath.normalize(vB8), super.getRadius().toString());
        Triad<Apfloat> vC9  = VectorMath.mult(VectorMath.normalize(vB9), super.getRadius().toString());
        Triad<Apfloat> vC10 = VectorMath.mult(VectorMath.normalize(vB10), super.getRadius().toString());
        Triad<Apfloat> vC11 = VectorMath.mult(VectorMath.normalize(vB11), super.getRadius().toString());
        Triad<Apfloat> vC12 = VectorMath.mult(VectorMath.normalize(vB12), super.getRadius().toString());
        Triad<Apfloat> vC13 = VectorMath.mult(VectorMath.normalize(vB13), super.getRadius().toString());

        // ==== RHOMBOHEDRAL FACES ====
        Tetrad<Tuple<Apfloat>> rho0 = new Tetrad<>(vC6, vC0, vC8, vC2);
        Triad<Apfloat> rho0_norm = VectorMath.normalQuad(vC6, vC0, vC8, vC2, true);
        Tetrad<Tuple<Apfloat>> rho1 = new Tetrad<>(vC6, vC2, vC7, vC4);
        Triad<Apfloat> rho1_norm = VectorMath.normalQuad(vC6, vC2, vC7, vC4, true);
        Tetrad<Tuple<Apfloat>> rho2 = new Tetrad<>(vC6, vC4, vC10, vC0);
        Triad<Apfloat> rho2_norm = VectorMath.normalQuad(vC6, vC4, vC10, vC0, true);
        Tetrad<Tuple<Apfloat>> rho3 = new Tetrad<>(vC9, vC1, vC7, vC2);
        Triad<Apfloat> rho3_norm = VectorMath.normalQuad(vC9, vC1, vC7, vC2, true);
        Tetrad<Tuple<Apfloat>> rho4 = new Tetrad<>(vC9, vC2, vC8, vC5);
        Triad<Apfloat> rho4_norm = VectorMath.normalQuad(vC9, vC2, vC8, vC5, true);
        Tetrad<Tuple<Apfloat>> rho5 = new Tetrad<>(vC9, vC5, vC13, vC1);
        Triad<Apfloat> rho5_norm = VectorMath.normalQuad(vC9, vC2, vC8, vC5, true);
        Tetrad<Tuple<Apfloat>> rho6 = new Tetrad<>(vC11, vC1, vC13, vC3);
        Triad<Apfloat> rho6_norm = VectorMath.normalQuad(vC11, vC1, vC13, vC3, true);
        Tetrad<Tuple<Apfloat>> rho7 = new Tetrad<>(vC11, vC3, vC10, vC4);
        Triad<Apfloat> rho7_norm = VectorMath.normalQuad(vC11, vC3, vC10, vC4, true);
        Tetrad<Tuple<Apfloat>> rho8 = new Tetrad<>(vC11, vC4, vC7, vC1);
        Triad<Apfloat> rho8_norm = VectorMath.normalQuad(vC11, vC4, vC7, vC1, true);
        Tetrad<Tuple<Apfloat>> rho9 = new Tetrad<>(vC12, vC0, vC10, vC3);
        Triad<Apfloat> rho9_norm = VectorMath.normalQuad(vC12, vC0, vC10, vC3, true);
        Tetrad<Tuple<Apfloat>> rho10 = new Tetrad<>(vC12, vC3, vC13, vC5);
        Triad<Apfloat> rho10_norm = VectorMath.normalQuad(vC12, vC3, vC13, vC5, true);
        Tetrad<Tuple<Apfloat>> rho11 = new Tetrad<>(vC12, vC5, vC8, vC0);
        Triad<Apfloat> rho11_norm = VectorMath.normalQuad(vC12, vC5, vC8, vC0, true);
        faces_rho = new Dodecad<>(
                rho0, rho1, rho2, rho3,
                rho4, rho5, rho6, rho7,
                rho8, rho9, rho10,rho11
        );
        face_norms_rho = new Dodecad<>(
                rho0_norm, rho1_norm, rho2_norm, rho3_norm,
                rho4_norm, rho5_norm, rho6_norm, rho7_norm,
                rho8_norm, rho9_norm, rho10_norm,rho11_norm
        );
    }

    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < faces_rho.fetchSize(); i++){

            // === Check rhombohedral faces ===
            // check if each point lies within bounds of each face

            // get vertices of vertA
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_rho.fetch(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);


            // given point p and vertA, calculate vector from vertA -> p:
            // m = p - vertA = (p_x - vertA_x, p_x - vertA_y, p_x - vertA_z)
            Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> face_norm = (Triad<Apfloat>) face_norms_rho.fetch(i);
            Apfloat d = VectorMath.dot_prod(face_norm, m);

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
