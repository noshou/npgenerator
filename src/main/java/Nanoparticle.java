import com.google.common.primitives.UnsignedInteger;
import org.apfloat.*;

/**
 * Abstract base class representing a nanoparticle composed of atoms.
 * <p>
 * This class provides core functionality for modeling nanoparticles, including:
 * <ul>
 *   <li>Atomic properties (radius, volume)</li>
 *   <li>Nanoparticle volume and atom count</li>
 *   <li>High-precision calculations using Apfloat</li>
 *   <li>Temperature-dependent properties</li>
 * </ul>
 *
 * Subclasses implement specific geometric models to calculate nanoparticle properties.
 * All floating-point calculations use arbitrary precision arithmetic.
 *
 * <h2>Core Parameters</h2>
 * <ul>
 *   <li>{@code r_atm}: Atomic radius (Å) - converted from angstroms during construction</li>
 *   <li>{@code v_atm}: Atomic volume (nm³)</li>
 *   <li>{@code n_atm}: Estimated number of atoms in nanoparticle</li>
 *   <li>{@code v_npt}: Nanoparticle volume (nm³)</li>
 *   <li>{@code temp}: Temperature in Kelvin</li>
 * </ul>
 */
public abstract class Nanoparticle {

    /**
     * Radius of a single atom in Å.
     */
    private Apfloat r_atm;

    /**
     * Volume of a single atom in cubic nanometers (nm³).
     */
    private Apfloat v_atm;

    /**
     * Estimated number of atoms in the nanoparticle (unitless count).
     */
    private Apfloat n_atm;

    /**
     * Number of significant digits used in all Apfloat computations.
     */
    private int precision;

    /***
     * name of atom
     */
    private String name;

    /**
     * volume of nanoparticle
     */
    private Apfloat v_npt;

    /**
     * additional parameters to add to toString method
     */
    private String[] additional_string_params;

    /**
     * number of digits to display when toString called
     */
    private int print_digits;


    /**
     * equilibrium temperature of nanoparticle
     */
    private int temp;

    /**
     * sets volume of nanoparticle
     * @param volume Volume of nanoparticle
     */
    protected void vNpt(Apfloat volume) {
        this.v_npt = volume;
    }

    /**
     * sets volume of atom
     * @param volume Volume of atom
     */
    protected void vAtm(Apfloat volume) {
        this.v_atm = volume;
    }

    /**
     * sets number of atoms in nanoparticle
     * @param number Number of nanoparticles in atom
     */
    protected void nAtm(Apfloat number) {
        this.n_atm = number;
    }

    /**
     * additional parameters to display when
     * Nanoparticle.toString() method called
     * @param params list of additional parameters
     */
    protected void addParams(String[] params) {
        int length = params.length;
        additional_string_params = new String [length];
        for (int i = 0; i < length; i++) {
            additional_string_params[i] = params[i];
        }
    }

    /**
     * @return Volume of nanoparticle
     */
    protected Apfloat volumeNptApfloat() {
        return this.v_npt;
    }

    /**
     * @return Volume of atom
     */
    protected Apfloat volumeAtmApfloat() {
        return this.v_atm;
    }

    /**
     * @return Radius of atom (in Å)
     */
    protected Apfloat radiusAtmApfloat() {
        return this.r_atm;
    }

    /**
     * @return Number of atoms in nanoparticle
     */
    protected Apfloat numberAtmApFloat() {
        return this.n_atm;
    }

    /**
     * Returns the configured decimal precision for all Apfloat calculations.
     *
     * @return the number of significant digits used in calculations
     */
    public int fetchPrecision() {
        return this.precision;
    }

    /**
     * @return Number of digits to display outputted values
     */
    public int fetchDisplayDigits() {
        return this.print_digits;
    }

    /**
     * Returns the volume of a single atom with specified significant digits.
     *
     * @return string representation of atomic volume (nm³)
     */
    public String atomVolume() {
        Apfloat fmt_vAtm = new Apfloat(
                this.v_atm.toString(),
                this.fetchDisplayDigits()
        );
        return fmt_vAtm.toString();
    }

    /**
     * Returns the estimated number of atoms with specified significant digits.
     *
     * @return string representation of estimated atom count
     */
    public String numberOfAtoms() {
        Apfloat fmt_nAtm = new Apfloat(
                this.n_atm.toString(),
                this.fetchDisplayDigits()
        );
        return fmt_nAtm.toString();
    }

    /**
     * Returns the volume of the nanoparticle with specified significant digits.
     *
     * @return string representation of nanoparticle volume (nm³)
     */
    public String volumeOfNanoparticle() {
        Apfloat fmt_vNpt = new Apfloat(
                this.v_npt.toString(),
                this.fetchDisplayDigits()
        );
        return fmt_vNpt.toString();
    }

    /**
     * Returns the radius of a single atom with specified significant digits.
     *
     * @return string representation of the atomic radius (nanometers)
     */
    public String atomRadius() {

        Apfloat fmt_rAtm = new Apfloat(
                this.radiusAtmApfloat().toString(),
                this.print_digits
        );
        return fmt_rAtm.toString();
    }

    /**
     * @return Name of nanoparticle
     */
    public String fetchName() {
        return this.name;
    }

    /**
     * @return Equilibrium temperature of nanoparticle
     */
    public int fetchTemp() { return this.temp; }

    /**
     * Constructs a Nanoparticle base instance.
     *
     * @param r_atm        Atomic radius in angstroms (Å)
     * @param precision    Number of significant digits for calculations
     * @param name         Atom name (mmCIF compatible)
     * @param temp         Temperature in Kelvin (K)
     * @param print_digs   Number of digits for string representation
     */
    public Nanoparticle(
            String r_atm,
            UnsignedInteger precision,
            String name,
            int temp,
            int print_digs
    ) {
        // store name of atoms making up nanoparticle
        this.name = name;

        this.temp = temp;

        // store precision as int
        this.precision = precision.intValue();

        // radius in Å 
        this.r_atm = new Apfloat(r_atm, this.precision);

        // store digits to display when printing
        this.print_digits = print_digs;
    }

    /**
     * Returns a formatted string summary of all key computed values with
     * default precision of 5 digits.
     *
     * @return summary string of nanoparticle and atomic parameters
     */
    public String toString() {
        StringBuilder np = new StringBuilder();
        np.append("\n");
        np.append(String.format("\nAtom name:             \t%s", this.fetchName()));
        np.append(String.format("\nVolume of Atom:        \t%s nm³", this.atomVolume()));
        np.append(String.format("\nRadius of Atom:        \t%s Å", this.atomRadius()));
        np.append(String.format("\nTemperature (Kelvin):  \t%d K", this.fetchTemp()));
        np.append(String.format("\nCalculation  Precision:\t%s digits", this.fetchPrecision()));
        np.append(String.format("\nAtoms  in Nanoparticle:\t%s", this.numberOfAtoms()));
        np.append(String.format("\nVolume of Nanoparticle:\t%s nm³", this.volumeOfNanoparticle()));
        for (String param: this.additional_string_params) {
            np.append(param);
        }
        np.append("\n");
        return np.toString();
    }

}
