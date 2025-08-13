package io.github.noshou.npg.shapes.archimedean;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.shapes.catalan.*;
import io.github.noshou.npg.shapes.Shape;
import io.github.noshou.tuple.Pentad;
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
 * Represents a <b>Snub Dodecahedron</b>
 * <p> An Archimedean solid with 60 vertices, 92 faces
 * (comprised of 80 triangles and 12 pentagons), and 150 edges.
 * <p> It is a chiral polyhedron, meaning it exists in two mirror-image forms:
 * the {@link DodecahedronSnubLevo} (left-handed) and {@link DodecahedronSnubDextro} variants.
 * <p> It is the dual of the {@link HexecontahedronPentagonal}.
 */
public abstract class DodecahedronSnub extends Shape {

    // 12 pentagonal faces
    protected ArrayList<Pentad<Tuple<Apfloat>>> faces_pnt= new ArrayList<>();
    protected ArrayList<Triad<Apfloat>> face_norms_pnt = new ArrayList<>();

    // 80 triangular faces
    protected ArrayList<Triad<Tuple<Apfloat>>> faces_tri = new ArrayList<>();
    protected ArrayList<Triad<Apfloat>> face_norms_tri = new ArrayList<>();

    // vertices
    protected ArrayList<Triad<Apfloat>> vertices = new ArrayList<>();


    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N6 = new Apfloat("6", super.precision);
    private final Apfloat N27 = new Apfloat("27", super.precision);
    private final Apfloat SQRT5 = ApfloatMath.sqrt(N5);

    // φ = (1 + sqrt(5) / 2)
    private final Apfloat PHI = (N1.add(SQRT5)).divide(N2);

    // φ² (1 + sqrt(5) / 2)² = 6 + 2 * sqrt(5)
    private final Apfloat PHISQR = N6.add(N2.multiply(SQRT5));

    // FRAC_PHI_over_N2 = φ/2
    private final Apfloat FRAC_PHI_over_N2 = PHI.divide(N2);

    // SQRT_FRAC_N5_over_N27_over_2 = sqrt(5/27)/2
    private final Apfloat SQRT_FRAC_PHI_sub_N5_over_N27_over_2 = ApfloatMath.sqrt(PHI.subtract(N5.divide(N27))).divide(N2);

    // CBRT_01 = cbrt(FRAC_PHI_over_N2 + SQRT_FRAC_N5_over_N27_over_2)
    private final Apfloat CBRT_01 = ApfloatMath.cbrt(FRAC_PHI_over_N2.add(SQRT_FRAC_PHI_sub_N5_over_N27_over_2));

    // CBRT_02 = cbrt(FRAC_PHI_over_N2 - SQRT_FRAC_N5_over_N27_over_2)
    private final Apfloat CBRT_02 = ApfloatMath.cbrt(FRAC_PHI_over_N2.subtract(SQRT_FRAC_PHI_sub_N5_over_N27_over_2));

    // CBRT = CBRT_01 - CBRT_02
    private final Apfloat CBRT = CBRT_01.add(CBRT_02);
    private final Apfloat CBRTSQR = ApfloatMath.pow(CBRT, 2);

    // PHI_mul_FRAC_SQRT_3_sub_CBRTSQR_over_N2 = φ * sqrt(3 - (CBRT²)) / 2
    private final Apfloat PHI_mul_FRAC_SQRT_3_sub_CBRTSQR_over_N2 = PHI.multiply(
            ApfloatMath.sqrt(N3.subtract(CBRTSQR))
    ).divide(N2);

    // PHI_mul_SQRT_CBRT_sub_1_sub_FRAC_N1_over_CBRT_mul_PHI_over_N2 = φ * sqrt((CBRT - 1 - (1/CBRT)) * φ) / 2
    private final Apfloat PHI_mul_SQRT_CBRT_sub_1_sub_FRAC_N1_over_CBRT_mul_PHI_over_N2 = PHI.multiply(
            ApfloatMath.sqrt(CBRT.subtract(N1).subtract(N1.divide(CBRT))).multiply(PHI)
    ).divide(N2);

    // FRAC_PHI_mul_SQRT_CBRT_mul_CBRT_add_PHI_add_N1_over_N2 = φ * sqrt(CBRT * (CBRT + φ) + 1) / 2
    private final Apfloat FRAC_PHI_mul_SQRT_CBRT_mul_CBRT_add_PHI_add_N1_over_N2 = PHI.multiply(
            ApfloatMath.sqrt(CBRT.multiply(CBRT.add(PHI)).add(N1)).divide(N2)
    );

    /**
     * @return φ * sqrt(3 - (CBRT²)) / 2
     */
    public Apfloat C0() {
        return PHI_mul_FRAC_SQRT_3_sub_CBRTSQR_over_N2;
    }

    /**
     * @return CBRT * φ * sqrt(3 - (CBRT²)) / 2
     */
    public Apfloat C1() {
        return CBRT.multiply(PHI_mul_FRAC_SQRT_3_sub_CBRTSQR_over_N2);
    }

    /**
     * @return φ * sqrt((CBRT - 1 - (1/CBRT)) * φ) / 2
     */
    public Apfloat C2() {
        return PHI_mul_SQRT_CBRT_sub_1_sub_FRAC_N1_over_CBRT_mul_PHI_over_N2;
    }

    /**
     * @return (CBRT²) * φ * sqrt(3 - (CBRT²)) / 2
     */
    public Apfloat C3() {
        return CBRTSQR.multiply(PHI_mul_FRAC_SQRT_3_sub_CBRTSQR_over_N2);
    }

    /**
     * @return CBRT * φ * sqrt((CBRT - 1 - (1/CBRT)) * φ) / 2
     */
    public Apfloat C4() {
        return CBRT.multiply(PHI_mul_SQRT_CBRT_sub_1_sub_FRAC_N1_over_CBRT_mul_PHI_over_N2);
    }

    /**
     * @return φ * sqrt(1 - CBRT + (1 + φ) / CBRT) / 2
     */
    public Apfloat C5() {
        return PHI.multiply(
                ApfloatMath.sqrt(N1.subtract(CBRT).add((N1.add(PHI)).divide(CBRT)))
        ).divide(N2);
    }

    /**
     * @return φ * sqrt(CBRT + 1 - φ) / 2
     */
    public Apfloat C6() {
        return PHI.multiply(
                ApfloatMath.sqrt(CBRT.add(N1).subtract(PHI))
        ).divide(N2);
    }

    /**
     * @return (CBRT²) * φ * sqrt((CBRT - 1 - (1/CBRT)) * φ) / 2
     */
    public Apfloat C7() {
        return CBRTSQR.multiply(PHI_mul_SQRT_CBRT_sub_1_sub_FRAC_N1_over_CBRT_mul_PHI_over_N2);
    }

    /**
     * @return CBRT * φ * sqrt(1 - CBRT + (1 + φ) / CBRT) / 2
     */
    public Apfloat C8() {
        return CBRT.multiply(C5());
    }

    /**
     * @return sqrt((CBRT + 2) * φ + 2) / 2
     */
    public Apfloat C9() {
        return ApfloatMath.sqrt(CBRT.add(N2).multiply(PHI).add(N2)).divide(N2);
    }

    /**
     * @return CBRT * sqrt(CBRT * (1 + φ) - φ) / 2
     */
    public Apfloat C10() {
        return CBRT.multiply(
                ApfloatMath.sqrt(CBRT.multiply(N1.add(PHI)).subtract(PHI))
        ).divide(N2);
    }

    /**
     * @return sqrt((CBRT²) * (1 + 2 * φ) - φ) / 2
     */
    public Apfloat C11() {
        return ApfloatMath.sqrt(
                CBRTSQR.multiply(N1.add(N2.multiply(PHI))).subtract(PHI)
        ).divide(N2);
    }

    /**
     * @return φ * sqrt(CBRT² + CBRT) / 2
     */
    public Apfloat C12() {
        return PHI.multiply(
                ApfloatMath.sqrt(CBRTSQR.add(CBRT))
        ).divide(N2);
    }

    /**
     * @return (φ²) * sqrt(CBRT * (CBRT + φ) + 1) / (2 * CBRT)
     */
    public Apfloat C13() {
        return PHISQR.multiply(
                ApfloatMath.sqrt(CBRT.multiply(CBRT.add(PHI)).add(N1))
        ).divide(N2.multiply(CBRT));
    }

    /**
     * @return φ * sqrt(CBRT * (CBRT + φ) + 1) / 2
     */
    public Apfloat C14() {
        return PHI.multiply(
                ApfloatMath.sqrt(CBRT.multiply(CBRT.add(PHI)).add(N1))
        ).divide(N2);
    }

    /**
     * @return -( φ * sqrt(3 - (CBRT²)) / 2)
     */
    public Apfloat NEG_C0() {
        return C0().multiply(NEG_N1);
    }

    /**
     * @return -( CBRT * φ * sqrt(3 - (CBRT²)) / 2)
     */
    public Apfloat NEG_C1() {
        return C1().multiply(NEG_N1);
    }

    /**
     * @return -( φ * sqrt((CBRT - 1 - (1/CBRT)) * φ) / 2)
     */
    public Apfloat NEG_C2() {
        return C2().multiply(NEG_N1);
    }

    /**
     * @return -((CBRT²) * φ * sqrt(3 - (CBRT²)) / 2)
     */
    public Apfloat NEG_C3() {
        return C3().multiply(NEG_N1);
    }

    /**
     * @return -( CBRT * φ * sqrt((CBRT - 1 - (1/CBRT)) * φ) / 2)
     */
    public Apfloat NEG_C4() {
        return C4().multiply(NEG_N1);
    }

    /**
     * @return -(φ * sqrt(1 - CBRT + (1 + φ) / CBRT) / 2)
     */
    public Apfloat NEG_C5() {
        return C5().multiply(NEG_N1);
    }

    /**
     * @return -(φ * sqrt(CBRT + 1 - φ) / 2)
     */
    public Apfloat NEG_C6() {
        return C6().multiply(NEG_N1);
    }

    /**
     * @return -( CBRT² * φ * sqrt((CBRT - 1 - (1/CBRT)) * φ) / 2)
     */
    public Apfloat NEG_C7() {
        return C7().multiply(NEG_N1);
    }

    /**
     * @return -(CBRT * φ * sqrt(1 - CBRT + (1 + φ) / CBRT) / 2)
     */
    public Apfloat NEG_C8 () {
        return C8().multiply(NEG_N1);
    }

    /**
     * @return  -(sqrt((CBRT+2)*φ+2) / 2)
     */
    public Apfloat NEG_C9() {
        return C9().multiply(NEG_N1);
    }

    /**
     * @return -(CBRT * sqrt(CBRT * (1 + φ) - φ) / 2)
     */
    public Apfloat NEG_C10() {
        return C10().multiply(NEG_N1);
    }

    /**
     * @return -(sqrt(CBRT² * (1 + 2 * φ) - φ) / 2)
     */
    public Apfloat NEG_C11() {
        return C11().multiply(NEG_N1);
    }

    /**
     * @return  -(φ* sqrt(CBRT² + CBRT) / 2)
     */
    public Apfloat NEG_C12() {
        return C12().multiply(NEG_N1);
    }

    /**
     * @return   -(φ² * sqrt(CBRT * (CBRT + φ) + 1) / 2 * 1/CBRT)
     */
    public Apfloat NEG_C13() {
        return C13().multiply(NEG_N1);
    }

    /**
     * @return  -( φ * sqrt(CBRT * (CBRT + φ) + 1) / 2)
     */
    public Apfloat NEG_C14() {
        return C14().multiply(NEG_N1);
    }

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
    public DodecahedronSnub(
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


    /**
     * Tests whether the given Cartesian point lies inside or on the boundary
     * of the snub dodecahedron.
     *
     * <p>The test is performed by checking dot products against all
     * face normals (pentagons and triangles).
     *
     * @param point_cart the Cartesian coordinate point
     * @return {@code true} if inside or on surface, {@code false} if outside
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_tri.size(); i++) {

            // === Check pentagonal faces ===
            if (i < faces_pnt.size()) {
                Pentad<Tuple<Apfloat>> face = (Pentad<Tuple<Apfloat>>) faces_pnt.get(i);
                Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);

                // given point p and vertA, calculate vector from centroid -> p:
                // m = p - centroid = (p_x - centroid_x, p_x - centroid_y, p_x -centroid_z)
                Triad<Apfloat> m = subs(point_cart, vertA);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_pnt.get(i);
                Apfloat d = dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(N0) > 0) {
                    return false;
                }
            }

            // === Check triangular faces ===
            Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.get(i);
            Triad<Apfloat> vert0 = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> vert1 = (Triad<Apfloat>) face.fetch(1);
            Triad<Apfloat> vert2 = (Triad<Apfloat>) face.fetch(2);

            // centroid of triangle
            Triad<Apfloat> centroid = new Triad<>(
                    vert0.fetch(0).add(vert1.fetch(0)).add(vert2.fetch(0)).divide(N3),
                    vert0.fetch(1).add(vert1.fetch(1)).add(vert2.fetch(1)).divide(N3),
                    vert0.fetch(2).add(vert1.fetch(2)).add(vert2.fetch(2)).divide(N3)
            );

            // given point p and vertA, calculate vector from centroid -> p:
            // m = p - centroid = (p_x - centroid_x, p_x - centroid_y, p_x -centroid_z)
            Triad<Apfloat> m = subs(point_cart, centroid);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.get(i);
            Apfloat d = dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
