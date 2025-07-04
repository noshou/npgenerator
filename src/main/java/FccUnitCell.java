import com.oson.tuple.*;
import org.jetbrains.annotations.*;
import org.apfloat.*;

/**
 * Represents a face-centered cubic (FCC) unit cell with a four-atom basis.
 * <p>
 * The FCC unit cell is defined with equal edge lengths (a = b = c) and orthogonal angles (90°),
 * and it places atoms at the canonical fractional coordinates:
 * <ul>
 *   <li>(0, 0, 0)</li>
 *   <li>(0.5, 0.5, 0)</li>
 *   <li>(0.5, 0, 0.5)</li>
 *   <li>(0, 0.5, 0.5)</li>
 * </ul>
 * These positions correspond to the typical FCC lattice symmetry with space group "F m -3 m".
 * </p>
 *
 * <p>
 * This class is suitable for constructing periodic crystals or finite FCC-based nanoparticles.
 * </p>
 */
public class FccUnitCell extends BravaisUnitCell {

    private static final Apfloat HALF = new Apfloat("0.5");

    /**
     * Constructs a fully parameterized FCC unit cell with a four-atom basis.
     * Each atom may have distinct identity and properties (element, radius, charge).
     *
     * @param lattice_constant the unit cell edge length in Å (shared by a, b, c); must not be null
     * @param precision         number of decimal digits for Apfloat calculations
     * @param atom_a_name       element name for atom at (0, 0, 0); must not be null
     * @param atom_a_charge     formal charge of atom A
     * @param atom_a_radius     radius of atom A (Å); must not be null
     * @param atom_b_name       element name for atom at (0.5, 0.5, 0); must not be null
     * @param atom_b_charge     formal charge of atom B
     * @param atom_b_radius     radius of atom B (Å); must not be null
     * @param atom_c_name       element name for atom at (0.5, 0, 0.5); must not be null
     * @param atom_c_charge     formal charge of atom C
     * @param atom_c_radius     radius of atom C (Å); must not be null
     * @param atom_d_name       element name for atom at (0, 0.5, 0.5); must not be null
     * @param atom_d_charge     formal charge of atom D
     * @param atom_d_radius     radius of atom D (Å); must not be null
     */
    public FccUnitCell(
            @NotNull Apfloat lattice_constant,
            int precision,
            @NotNull String atom_a_name,
            int atom_a_charge,
            @NotNull String atom_a_radius,
            @NotNull String atom_b_name,
            int atom_b_charge,
            @NotNull String atom_b_radius,
            @NotNull String atom_c_name,
            int atom_c_charge,
            @NotNull String atom_c_radius,
            @NotNull String atom_d_name,
            int atom_d_charge,
            @NotNull String atom_d_radius
    ) {
        super(
                lattice_constant, lattice_constant, lattice_constant,    // a = b = c
                new Apfloat("90"),
                new Apfloat("90"),
                new Apfloat("90"),                                  // α = β = γ = 90°
                "F m -3 m",                                              // FCC space group
                buildBasis(
                        precision,
                        atom_a_name, atom_a_charge, atom_a_radius,
                        atom_b_name, atom_b_charge, atom_b_radius,
                        atom_c_name, atom_c_charge, atom_c_radius,
                        atom_d_name, atom_d_charge, atom_d_radius
                ),
                LatticeType.FCC,
                precision
        );
    }

    /**
     * Builds a four-atom basis located at canonical FCC fractional positions.
     * Each atom is initialized with its identity, radius, formal charge, and high-precision settings.
     *
     * @param precision      precision for Apfloat volume/radius math
     * @param atom_a_name    element name for atom A at (0, 0, 0); must not be null
     * @param atom_a_charge  formal charge of atom A
     * @param atom_a_radius  atomic radius of atom A; must not be null
     * @param atom_b_name    element name for atom B at (0.5, 0.5, 0); must not be null
     * @param atom_b_charge  formal charge of atom B
     * @param atom_b_radius  atomic radius of atom B; must not be null
     * @param atom_c_name    element name for atom C at (0.5, 0, 0.5); must not be null
     * @param atom_c_charge  formal charge of atom C
     * @param atom_c_radius  atomic radius of atom C; must not be null
     * @param atom_d_name    element name for atom D at (0, 0.5, 0.5); must not be null
     * @param atom_d_charge  formal charge of atom D
     * @param atom_d_radius  atomic radius of atom D; must not be null
     * @return a {@code Polyad<Atom>} containing the four initialized atoms
     */
    @Contract("_, _, _, _, _, _, _, _, _, _, _, _, _ -> new")
    private static @NotNull Polyad<@NotNull Atom> buildBasis(
            int precision,
            @NotNull String atom_a_name,
            int atom_a_charge,
            @NotNull String atom_a_radius,
            @NotNull String atom_b_name,
            int atom_b_charge,
            @NotNull String atom_b_radius,
            @NotNull String atom_c_name,
            int atom_c_charge,
            @NotNull String atom_c_radius,
            @NotNull String atom_d_name,
            int atom_d_charge,
            @NotNull String atom_d_radius
    ) {
        Atom a = new Atom(atom_a_name, atom_a_radius, new Triad<>("0", "0", "0"), atom_a_charge, precision);
        Atom b = new Atom(atom_b_name, atom_b_radius, new Triad<>(HALF.toString(), HALF.toString(), "0"), atom_b_charge, precision);
        Atom c = new Atom(atom_c_name, atom_c_radius, new Triad<>(HALF.toString(), "0", HALF.toString()), atom_c_charge, precision);
        Atom d = new Atom(atom_d_name, atom_d_radius, new Triad<>("0", HALF.toString(), HALF.toString()), atom_d_charge, precision);
        return new Polyad<>(new Atom[]{a, b, c, d});
    }

    /**
     * Returns the atom located at a given fractional coordinate within the FCC unit cell.
     * <p>
     * This method matches the input against the four canonical positions.
     * If no atom exists at the given coordinate, {@code null} is returned.
     * </p>
     *
     * @param frac_x fractional x-coordinate (Apfloat, ≥ 0, < 1)
     * @param frac_y fractional y-coordinate (Apfloat, ≥ 0, < 1)
     * @param frac_z fractional z-coordinate (Apfloat, ≥ 0, < 1)
     * @return the {@link Atom} located at this lattice position, or {@code null} if empty
     */
    @Override
    @Contract(pure = true)
    public @Nullable Atom getLatticePoint(
            @NotNull Apfloat frac_x,
            @NotNull Apfloat frac_y,
            @NotNull Apfloat frac_z
    ) {
        // Normalize fractional coords modulo 1 (assuming inputs ≥ 0)
        Apfloat x = frac_x.mod(Apfloat.ONE);
        Apfloat y = frac_y.mod(Apfloat.ONE);
        Apfloat z = frac_z.mod(Apfloat.ONE);

        if (x.compareTo(Apfloat.ZERO) == 0
                && y.compareTo(Apfloat.ZERO) == 0
                && z.compareTo(Apfloat.ZERO) == 0) {
            return super.getAtom(0);
        } else if (x.compareTo(HALF) == 0
                && y.compareTo(HALF) == 0
                && z.compareTo(Apfloat.ZERO) == 0) {
            return super.getAtom(1);
        } else if (x.compareTo(HALF) == 0
                && y.compareTo(Apfloat.ZERO) == 0
                && z.compareTo(HALF) == 0) {
            return super.getAtom(2);
        } else if (x.compareTo(Apfloat.ZERO) == 0
                && y.compareTo(HALF) == 0
                && z.compareTo(HALF) == 0) {
            return super.getAtom(3);
        } else {
            return null;
        }
    }

//    // unit test for debug purposes
//    // based on gold nanoparticles
//    public static void main(String[] args) {
//        Apfloat lattice_constant = new Apfloat("")
//    }
}
