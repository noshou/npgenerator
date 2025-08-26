package io.github.noshou.npg.atom;
import io.github.noshou.tuple.*;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;

/**
 * Represents a single atom with associated identity, fractional position,
 * radius, volume, and charge.
 * <p>
 * The radius and derived values like volume are stored as {@code String} to
 * support arbitrary-precision
 * Apfloat arithmetic. Atom.Atom data is intended to be converted into cartesian
 * coordinates and written into
 * mmCIF files or other crystallographic formats.
 * <p>
 * Additional fields like index and cartesian are mutable and assist
 * with serialization or spatial placement.
 * 
 */
public class Atom implements Comparable<Atom> {

    /** Chemical element symbol (e.g., "Au") */
    private final @NotNull String element;

    /** Atomic radius in Ångströms (as Apfloat-compatible string) */
    private final @NotNull String radius;

    /** Precomputed atomic volume from radius (as Apfloat-compatible string) */
    private final @NotNull String volume;

    /** Formal charge, stored as string (e.g., "1+", "1-", "0") */
    private final @NotNull String formal_charge;

    /** Fractional coordinates of the atom in the unit cell */
    private final @NotNull Triad<String> fractional_position;

    /** Precision (digits) used for all Apfloat computations */
    private final int precision;

    /** Unique index for this atom in output (default: -1) */
    private int index = -1;

    /** Cartesian coordinates of the atom (default: Nullad) */
    private @NotNull Tuple<String> cartesian = new Nullad<>();

    /** Fractional coordinates of the atom (default: Nullad) */
    private @NotNull Tuple<String> fractional = new Nullad<>();

    /**
     * Constructs an {@code Atom.Atom} with element name, radius, charge, fractional position.
     *
     * @param element     the atomic symbol (e.g., "Na"); must not be null
     * @param radius      atomic radius in Ångströms as a string (Apfloat-compatible)
     * @param frac_pos    fractional coordinates (x, y, z) as strings in a triad
     * @param charge      formal charge of the atom (positive, negative, or 0)
     * @param precision   decimal precision to use for volume calculations
     */
    public Atom(
            @NotNull String element,
            @NotNull String radius,
            @NotNull Triad<String> frac_pos,
            int charge,
            int precision
    ) {
        if (charge > 0) {
            this.formal_charge = String.format("+%d", charge);
        } else if (charge < 0) {
            this.formal_charge = String.format("-%d", Math.abs(charge));
        } else {
            this.formal_charge = "0";
        }

        // element must conform to first letter being upercase, second lower case
        if (element.length() == 1) {
            element = element.trim();
            this.element = Character.toString(element.charAt(0)).toUpperCase();
        } else if (element.length() == 2){
            element = element.trim();
            this.element = Character.toString(element.charAt(0)).toUpperCase()
                            + Character.toString(element.charAt(1)).toLowerCase();
        } else {
            throw new IllegalArgumentException("Invalid element symbol!");
        }
        this.radius = radius;
        this.fractional_position = frac_pos;
        this.precision = precision;

        // Volume = (4/3) * π * (r)^3
        Apfloat PI = ApfloatMath.pi(precision);
        Apfloat THREE = new Apfloat("3", precision);
        Apfloat FOUR = new Apfloat("4", precision);
        Apfloat r = new Apfloat(this.radius, precision);
        Apfloat v_atm_1 = ApfloatMath.pow(r, 3);
        Apfloat v_atm_2 = FOUR.multiply(PI).divide(THREE);
        this.volume = (v_atm_1.multiply(v_atm_2)).toString();
    }

    /**
     * Sets the Cartesian coordinates of the atom.
     * <p>
     * This method is <b>private</b> because the cartesian should not be known at the time of
     * {@code Atom.Atom} construction—it is computed later when converting fractional positions
     * to Cartesian coordinates, typically by a lattice or grid-building class.
     * @param the converted (x, y, z) position in real space relative to the nanoparticle
     */
    @Contract(mutates = "this")
    private void setCartesian(@NotNull Triad<String> cartesian) {
        this.cartesian = cartesian;
    }

    /**
     * Sets the fractional coordinates of the atom.
     * <p>
     * This method is <b>private</b> because the fractional coordinate should not be known at the time of
     * {@code Atom.Atom} construction—it is computed later.
     * @param the converted (x, y, z) position in fractional space relative to the crystal
     */
    @Contract(mutates = "this")
    private void setFractional(@NotNull Triad<String> cartesian) {
        this.cartesian = cartesian;
    }


    /**
     * Gets the Cartesian coordinates of the atom, if set.
     * Default is a {@code Nullad} if unset.
     * @return a {@code Tuple<String>} of coordinates (may be empty/null-like)
     */
    @Contract(pure = true)
    public @NotNull Tuple<String> getCartesian() {
        return this.cartesian;
    }

    /**
     * Gets the fractional coordinates of the atom, if set.
     * Default is a {@code Nullad} if unset.
     * @return a {@code Tuple<String>} of coordinates (may be empty/null-like)
     */
    @Contract(pure = true)
    public @NotNull Tuple<String> getFractional() {
        return this.fractional;
    }

    /**
     * Sets a unique index identifier for this atom.
     * <p>
     * This method is <b>private</b> because the index is not intrinsic to the atom's definition,
     * but rather assigned externally when atoms are laid out or serialized (e.g., for mmCIF export).
     * It ensures the {@code Atom.Atom} constructor remains minimal and strictly about physical properties.
     * @param idx unique index (typically 1-indexed)
     */
    @Contract(mutates = "this")
    private void setIndex(int idx) {
        this.index = idx;
    }

    /**
     * Returns the atomic volume computed from the radius in Å³.
     * @return the volume as a string (Apfloat-compatible)
     */
    @Contract(pure = true)
    public @NotNull String getVolume() {
        return this.volume;
    }

    /**
     * Gets the formal charge (e.g., "-2", "+1", "0").
     * @return the formal charge as a string
     */
    @Contract(pure = true)
    public @NotNull String getFormalCharge() {
        return this.formal_charge;
    }

    /**
     * Returns formal charge as integer
     * @return formal charge as integer
     */
    @Contract(pure = true)
    public int getFormalChargeInt() {
        if (this.formal_charge.startsWith("-")) {
            return -1 * Integer.parseInt(this.formal_charge.replace("-", ""));
        }
        else if (this.formal_charge.startsWith("+")) {
            return Integer.parseInt(this.formal_charge.replace("+", ""));
        }
        else {
            return 0;
        }
    }

    /**
     * Returns the atomic element symbol.
     * @return a non-null chemical symbol
     */
    @Contract(pure = true)
    public @NotNull String getElement() {
        return this.element;
    }

    /**
     * Returns the atomic radius (string form for Apfloat compatibility).
     * @return radius in Å as a string
     */
    @Contract(pure = true)
    public @NotNull String getRadius() {
        return this.radius;
    }

    /**
     * Gets the unique index of the atom.
     * Defaults to -1 if not explicitly set.
     * @return the integer index
     */
    @Contract(pure = true)
    public int getIndex() {
        return this.index;
    }

    /**
     * Returns the precision (in decimal digits) used for Apfloat math.
     * @return numeric precision
     */
    @Contract(pure = true)
    public int getPrecision() {
        return this.precision;
    }

    /**
     * Compares two atoms based on atomic radius using Apfloat math.
     * @param that the atom to compare to
     * @return -1, 0, or 1 depending on comparison
     */
    @Override
    @Contract(value = "_ -> param1", pure = true)
    public int compareTo(@NotNull Atom that) {
        int maxPrecision = Math.max(this.getPrecision(), that.getPrecision());
        Apfloat this_radius = new Apfloat(this.getRadius(), maxPrecision);
        Apfloat that_radius = new Apfloat(that.getRadius(), maxPrecision);
        return this_radius.compareTo(that_radius);
    }

    /**
     * Returns a human-readable string representation of the atom.
     * <p>
     * This includes:
     * <ul>
     *     <li>Element symbol</li>
     *     <li>Radius in Ångströms</li>
     *     <li>Computed atomic volume in Å³</li>
     *     <li>Formal charge</li>
     *     <li>Fractional lattice position</li>
     *     <li>Precision used for numerical calculations</li>
     * </ul>
     * Useful for debugging or displaying atom metadata in logs or user interfaces.
     * @return a formatted multiline string describing the atom’s key properties
     */
    @Override
    @Contract(pure = true)
    public @NotNull String toString() {
        return  String.format("Element:            \t%s\n", this.getElement()) +
                String.format("Radius:             \t%s Å\n", this.getRadius()) +
                String.format("Volume:             \t%s Å³\n", this.getVolume()) +
                String.format("Charge:             \t%s\n", this.getFormalCharge()) +
                String.format("Fractional Position:\t%s\n", this.getFractionalPosition()) +
                String.format("Digits of Precision:\t%d", this.getPrecision());
    }

    /**
     * Returns the atom's fractional coordinates in the unit cell.
     * @return a {@code Triad<String>} containing x, y, z as strings
     */
    @Contract(pure = true)
    public @NotNull Triad<String> getFractionalPosition() {
        return fractional_position;
    }

    /**
     * Assigns lattice metadata to this atom: a unique index and fractional and cartesian coordinates.
     * <p>
     * <strong>Warning:</strong> This method directly mutates the calling object by updating its
     * index and real-space (Cartesian) position. It should be used with caution in contexts
     * where atomic identity is preserved but placement within the crystal or nanoparticle is staged.
     * <p><strong>Typical usage pattern:</strong>
     * <ol>
     *   <li>Create a unit cell (defines element types and fractional positions)</li>
     *   <li>Instantiate a lattice grid (defines spatial positions)</li>
     *   <li>Use {@code latticePoint()} to assign index and coordinates</li>
     * </ol>
     * @param idx         the globally unique atom index (typically 1-based)
     * @param cartesian  the Cartesian coordinates of the atom on the lattice
     * @param fractional the fractional coordinates of the atom on the lattice (must be between 0 and 1)
     */
    @Contract(mutates = "this")
    public void latticePoint(
            int idx, @NotNull Triad<String> cartesian, @NotNull Triad<String> fractional) {
        this.setIndex(idx);
        this.setCartesian(cartesian);
        this.setFracCoords(fractional);
    }
}
