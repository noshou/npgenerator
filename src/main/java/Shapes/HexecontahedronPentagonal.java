package Shapes;

import Atom.Atom;
import Lattice.LatticeType;
import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;

/**
 * Represents a Pentagonal Hexecontahedron, a convex polyhedron composed of 60 pentagonal faces.
 * This Catalan solid is the dual of the truncated icosidodecahedron and exhibits icosahedral symmetry.
 * It has 92 vertices and 150 edges, with coordinates based on the golden ratio (phi) and a special
 * parameter x, calculated here with high-precision Apfloat arithmetic.
 *
 * The shape can exist in two chiral forms: dextro (right-handed) and levo (left-handed), which are
 * mirror images of each other but not superimposable, similar to left and right hands.
 */
@SuppressWarnings("FieldCanBeLocal")
public abstract class HexecontahedronPentagonal extends Shape {

    // children must fill in these faces!

    // 60 pentagonal faces
    private ArrayList<Pentad<Tuple<Apfloat>>> faces_pnt = new ArrayList<>();
    private final ArrayList<Triad<Apfloat>> face_norms_pnt= new ArrayList<>();

    // 92 vertices
    protected ArrayList<Triad<Apfloat>> vertices = new ArrayList<>();

    // Apfloat constants
    private final Apfloat NEG_N1 = new Apfloat("-1", super.precision);
    private final Apfloat N0 = new Apfloat("0", super.precision);
    private final Apfloat N1 = new Apfloat("1", super.precision);
    private final Apfloat N2 = new Apfloat("2", super.precision);
    private final Apfloat N3 = new Apfloat("3", super.precision);
    private final Apfloat N4 = new Apfloat("4", super.precision);
    private final Apfloat N5 = new Apfloat("5", super.precision);
    private final Apfloat N62 = new Apfloat("62", super.precision);
    private final Apfloat N157 = new Apfloat("157", super.precision);
    private final Apfloat N225 = new Apfloat("225", super.precision);
    private final Apfloat N249 = new Apfloat("249", super.precision);
    private final Apfloat N392 = new Apfloat("392", super.precision);
    private final Apfloat N470 = new Apfloat("470", super.precision);
    private final Apfloat N627 = new Apfloat("627", super.precision);
    private final Apfloat N617 = new Apfloat("617", super.precision);
    private final Apfloat N670 = new Apfloat("670", super.precision);
    private final Apfloat N784 = new Apfloat("784", super.precision);
    private final Apfloat N842 = new Apfloat("842", super.precision);
    private final Apfloat N919 = new Apfloat("919", super.precision);
    private final Apfloat N1589 = new Apfloat("1589", super.precision);
    private final Apfloat FRAC_N5_over_N27 = (N5).divide(new Apfloat("27", super.precision));
    private final Apfloat PHI = (N1.add(ApfloatMath.sqrt(N5))).divide(N2); // golden ratio = (1 + sqrt(5))/2
    private final Apfloat PHISQR = ApfloatMath.pow(PHI, 2);
    private final Apfloat PHICUB = ApfloatMath.pow(PHI, 3);
    private final Apfloat FRAC_PHI_over_N2 = PHI.divide(N2);
    private final Apfloat FRAC_PHIminFRAC_5overN27_over_N2 = ApfloatMath.sqrt(
            PHI.subtract(FRAC_N5_over_N27)
    ).divide(N2); // sqrt(phi-(5/27)/2
    private final Apfloat x_01 = ApfloatMath.cbrt(
            FRAC_PHI_over_N2.add(FRAC_PHIminFRAC_5overN27_over_N2)
    ); // x_01 = cbrt(FRAC_PHI_over_N2 + FRAC_PHIminFRAC_5overN27_over_N2)
    private final Apfloat x_02 = ApfloatMath.cbrt(
            FRAC_PHI_over_N2.subtract(FRAC_PHIminFRAC_5overN27_over_N2)
    ); // x_02 = cbrt(FRAC_PHI_over_N2 - FRAC_PHIminFRAC_5overN27_over_N2)
    private final Apfloat X = x_01.add(x_02);   // x = x_01 + x_02
    private final Apfloat XSQR = ApfloatMath.pow(X, 2);
    private final Apfloat SQRT_3minXSQR = ApfloatMath.sqrt(N3.subtract(XSQR)); // sqrt(3 - x²)
    private final Apfloat SQRT_xminN1minFRAC_N1overX_mulPHI = ApfloatMath.sqrt(
            (X.subtract(N1).subtract(N1.divide(X))).multiply(PHI)
    ); // sqrt((x - 1 - 1/x) * phi)
    private final Apfloat SQRT_N1minXplusFRAC_N1plusPHI_over_X = ApfloatMath.sqrt(
            N1.subtract(X).add((N1.add(PHI)).divide(X))
    ); // sqrt(1 - x + (1 + phi) / x)
    private final Apfloat SQRT_x_mul_xplusPHI_plus1 = ApfloatMath.sqrt(
            X.multiply(X.add(PHI)).add(N1)
    ); // sqrt(x(x + phi) + 1)

    private final Apfloat SQRT_xplus2_mul_PHI_plus2 = ApfloatMath.sqrt(
            (X.add(N2)).multiply(PHI).add(N2)
    ); // sqrt((x + 2) * phi + 2)

    private final Apfloat SQRT_negXSQR_mul_2plusPHI_plus_x_mul_1plus3PHI_plus4 = ApfloatMath.sqrt(
            XSQR.negate().multiply(N2.add(PHI)).add(
                    X.multiply(N1.add(N3.multiply(PHI))).add(N4)
            )
    ); // sqrt(-x²*(2 + phi) + x*(1 + 3*phi) + 4)

    private final Apfloat SQRT_1plusFRAC_N1overX = ApfloatMath.sqrt(
            N1.add(N1.divide(X))
    ); // sqrt(1 + (1/x))

    private final Apfloat SQRT_2plus3PHI_min2X_plus_3overX = ApfloatMath.sqrt(
            N2.add(N3.multiply(PHI)).subtract(N2.multiply(X)).add(N3.divide(X))
    ); // sqrt(2 + 3*phi - 2*x + (3/x))

    private final Apfloat SQRT_poly_C10 = ApfloatMath.sqrt(
            XSQR.multiply(N392.add(N225.multiply(PHI)))
                    .add(X.multiply(N249.add(N670.multiply(PHI))))
                    .add(N470.add(N157.multiply(PHI)))
    ); // sqrt(x²*(392 + 225*phi) + x*(249 + 670*phi) + (470 + 157*phi))

    private final Apfloat SQRT_xsqr_plus_x_plus_1_plus_phi = ApfloatMath.sqrt(
            XSQR.add(X).add(N1).add(PHI)
    ); // sqrt(x² + x + 1 + phi)

    private final Apfloat SQRT_xsqr_plus_2xphi_plus_2 = ApfloatMath.sqrt(
            XSQR.add(X.multiply(N2).multiply(PHI)).add(N2)
    ); // sqrt(x² + 2*x*phi + 2)

    private final Apfloat SQRT_xsqr_mul_1plus2PHI_minPHI = ApfloatMath.sqrt(
            XSQR.multiply(N1.add(N2.multiply(PHI))).subtract(PHI)
    ); // sqrt(x²*(1 + 2*phi) - phi)

    private final Apfloat SQRT_xsqr_plus_x = ApfloatMath.sqrt(
            XSQR.add(X)
    ); // sqrt(x² + x)

    private final Apfloat SQRT_poly_C17 = ApfloatMath.sqrt(
            XSQR.multiply(N617.add(N842.multiply(PHI)))
                    .add(X.multiply(N919.add(N1589.multiply(PHI))))
                    .add(N627.add(N784.multiply(PHI)))
    ); // sqrt(x²*(617 + 842*phi) + x*(919 + 1589*phi) + (627 + 784*phi))

    /**
     * C0 = phi * sqrt(3 - x²) / 2
     */
    public Apfloat C0() {
        return PHI.multiply(SQRT_3minXSQR).divide(N2);
    }
    /** Negation of C0 */
    public Apfloat NEG_C0() {
        return C0().multiply(NEG_N1);
    }

    /**
     * returns 0
     */
    public Apfloat N0() {return N0;}
    /**
     * C1 = phi * sqrt((x - 1 - 1/x) * phi) / (2 * x)
     */
    public Apfloat C1() {
        return PHI.multiply(SQRT_xminN1minFRAC_N1overX_mulPHI).divide(N2.multiply(X));
    }
    /** Negation of C1 */
    public Apfloat NEG_C1() {
        return C1().multiply(NEG_N1);
    }

    /**
     * C2 = phi * sqrt((x - 1 - 1/x) * phi) / 2
     */
    public Apfloat C2() {
        return PHI.multiply(SQRT_xminN1minFRAC_N1overX_mulPHI).divide(N2);
    }
    /** Negation of C2 */
    public Apfloat NEG_C2() {
        return C2().multiply(NEG_N1);
    }

    /**
     * C3 = x² * phi * sqrt(3 - x²) / 2
     */
    public Apfloat C3() {
        return XSQR.multiply(PHI).multiply(SQRT_3minXSQR).divide(N2);
    }
    /** Negation of C3 */
    public Apfloat NEG_C3() {
        return C3().multiply(NEG_N1);
    }

    /**
     * C4 = phi * sqrt(1 - x + (1 + phi)/x) / 2
     */
    public Apfloat C4() {
        return PHI.multiply(SQRT_N1minXplusFRAC_N1plusPHI_over_X).divide(N2);
    }
    /** Negation of C4 */
    public Apfloat NEG_C4() {
        return C4().multiply(NEG_N1);
    }

    /**
     * C5 = sqrt(x * (x + phi) + 1) / (2 * x)
     */
    public Apfloat C5() {
        return SQRT_x_mul_xplusPHI_plus1.divide(N2.multiply(X));
    }
    /** Negation of C5 */
    public Apfloat NEG_C5() {
        return C5().multiply(NEG_N1);
    }

    /**
     * C6 = sqrt((x + 2) * phi + 2) / (2 * x)
     */
    public Apfloat C6() {
        return SQRT_xplus2_mul_PHI_plus2.divide(N2.multiply(X));
    }
    /** Negation of C6 */
    public Apfloat NEG_C6() {
        return C6().multiply(NEG_N1);
    }

    /**
     * C7 = sqrt(-x²*(2 + phi) + x*(1 + 3*phi) + 4) / 2
     */
    public Apfloat C7() {
        return SQRT_negXSQR_mul_2plusPHI_plus_x_mul_1plus3PHI_plus4.divide(N2);
    }
    /** Negation of C7 */
    public Apfloat NEG_C7() {
        return C7().multiply(NEG_N1);
    }

    /**
     * C8 = (1 + phi) * sqrt(1 + 1/x) / (2 * x)
     */
    public Apfloat C8() {
        return (N1.add(PHI)).multiply(SQRT_1plusFRAC_N1overX).divide(N2.multiply(X));
    }
    /** Negation of C8 */
    public Apfloat NEG_C8() {
        return C8().multiply(NEG_N1);
    }

    /**
     * C9 = sqrt(2 + 3*phi - 2*x + 3/x) / 2
     */
    public Apfloat C9() {
        return SQRT_2plus3PHI_min2X_plus_3overX.divide(N2);
    }
    /** Negation of C9 */
    public Apfloat NEG_C9() {
        return C9().multiply(NEG_N1);
    }

    /**
     * C10 = sqrt(x²*(392 + 225*phi) + x*(249 + 670*phi) + (470 + 157*phi)) / 62
     */
    public Apfloat C10() {
        return ApfloatMath.sqrt(
                XSQR.multiply(N392.add(N225.multiply(PHI)))
                        .add(X.multiply(N249.add(N670.multiply(PHI))))
                        .add(N470.add(N157.multiply(PHI)))
        ).divide(N62);
    }
    /** Negation of C10 */
    public Apfloat NEG_C10() {
        return C10().multiply(NEG_N1);
    }

    /**
     * C11 = phi * sqrt(x * (x + phi) + 1) / (2 * x)
     */
    public Apfloat C11() {
        return PHI.multiply(SQRT_x_mul_xplusPHI_plus1).divide(N2.multiply(X));
    }
    /** Negation of C11 */
    public Apfloat NEG_C11() {
        return C11().multiply(NEG_N1);
    }

    /**
     * C12 = phi * sqrt(x² + x + 1 + phi) / (2 * x)
     */
    public Apfloat C12() {
        return PHI.multiply(SQRT_xsqr_plus_x_plus_1_plus_phi).divide(N2.multiply(X));
    }

    /** Negation of C12 */
    public Apfloat NEG_C12() {
        return C12().multiply(NEG_N1);
    }
    /**
     * C13 = phi * sqrt(x² + 2*x*phi + 2) / (2 * x)
     */
    public Apfloat C13() {
        return PHI.multiply(SQRT_xsqr_plus_2xphi_plus_2).divide(N2.multiply(X));
    }

    /** Negation of C13 */
    public Apfloat NEG_C13() {
        return C13().multiply(NEG_N1);
    }

    /**
     * C14 = sqrt(x²*(1 + 2*phi) - phi) / 2
     */
    public Apfloat C14() {
        return SQRT_xsqr_mul_1plus2PHI_minPHI.divide(N2);
    }

    /** Negation of C14 */
    public Apfloat NEG_C14() {
        return C14().multiply(NEG_N1);
    }

    /**
     * C15 = phi * sqrt(x² + x) / 2
     */
    public Apfloat C15() {
        return PHI.multiply(SQRT_xsqr_plus_x).divide(N2);
    }

    /** Negation of C15 */
    public Apfloat NEG_C15() {
        return C15().multiply(NEG_N1);
    }

    /**
     * C16 = phi³ * sqrt(x * (x + phi) + 1) / (2 * x²)
     */
    public Apfloat C16() {
        return PHICUB.multiply(SQRT_x_mul_xplusPHI_plus1).divide(N2.multiply(XSQR));
    }
    /** Negation of C16 */
    public Apfloat NEG_C16() {
        return C16().multiply(NEG_N1);
    }

    /**
     * C17 = sqrt(x²*(617 + 842*phi) + x*(919 + 1589*phi) + (627 + 784*phi)) / 62
     */
    public Apfloat C17() {
        return ApfloatMath.sqrt(
                XSQR.multiply(N617.add(N842.multiply(PHI)))
                        .add(X.multiply(N919.add(N1589.multiply(PHI))))
                        .add(N627.add(N784.multiply(PHI)))
        ).divide(N62);
    }

    /** Negation of C17 */
    public Apfloat NEG_C17() {
        return C17().multiply(NEG_N1);
    }

    /**
     * C18 = phi² * sqrt(x * (x + phi) + 1) / (2 * x)
     */
    public Apfloat C18() {
        return PHISQR.multiply(SQRT_x_mul_xplusPHI_plus1).divide(N2.multiply(X));
    }

    /** Negation of C18 */
    public Apfloat NEG_C18() {
        return C18().multiply(NEG_N1);
    }

    /**
     * C19 = phi * sqrt(x * (x + phi) + 1) / 2
     */
    public Apfloat C19() {
        return PHI.multiply(SQRT_x_mul_xplusPHI_plus1).divide(N2);
    }
    /** Negation of C19 */
    public Apfloat NEG_C19() {
        return C19().multiply(NEG_N1);
    }

    public HexecontahedronPentagonal(
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

    /** fills scaled vertices */
    protected void scaledVert(ArrayList<Triad<Apfloat>> basis) {
        for (Triad<Apfloat> vB_i: basis) {
            vertices.add(VectorMath.mult(VectorMath.normalize(vB_i), super.getRadius().toString()));
        }
    }

    /** fills faces_pnt and face_norms_pnt vertices */
    protected void pntFaces(ArrayList<Pentad<Tuple<Apfloat>>> faces) {
        faces_pnt = faces;
        for (Pentad<Tuple<Apfloat>> pnt_i: faces) {
            face_norms_pnt.add(VectorMath.normalPent(
                    (Triad<Apfloat>) pnt_i.fetch(0),
                    (Triad<Apfloat>) pnt_i.fetch(1),
                    (Triad<Apfloat>) pnt_i.fetch(2),
                    (Triad<Apfloat>) pnt_i.fetch(3),
                    (Triad<Apfloat>) pnt_i.fetch(4),
                    true
                    )
            );
        }
    }

    /**
     * Tests whether the given Cartesian point lies inside or on the boundary
     *
     * <p>The test is performed by checking dot products against all
     * face normals (squares and triangles).</p>
     *
     * @param point_cart the Cartesian coordinate point
     * @return {@code true} if inside or on surface, {@code false} if outside
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {
        for (int i = 0; i < faces_pnt.size(); i++) {
            // === Check pentagonal faces ===
            Pentad<Tuple<Apfloat>> face = faces_pnt.get(i);
            Triad<Apfloat> vertA = (Triad<Apfloat>) face.fetch(0);

            // given point p and vertA, calculate vector from centroid -> p:
            // m = p - vertA = (p_x - centroid_x, p_x - centroid_y, p_x -centroid_z)
            Triad<Apfloat> m = VectorMath.subs(point_cart, vertA);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_pnt.get(i);
            Apfloat d = VectorMath.dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(N0) > 0) {
                return false;
            }
        }
        return true;
    }
}
