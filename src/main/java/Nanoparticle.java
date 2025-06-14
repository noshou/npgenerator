import com.google.common.primitives.UnsignedInteger;
import org.apfloat.*;

public abstract class Nanoparticle {

    /**
     * Radius of a single atom in nanometers (nm), converted internally from angstroms.
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

    private Apfloat v_npt;

    /**
     * temperature of atom
     */
    private int temp;

    protected void vNpt(Apfloat volume) {
        this.v_npt = volume;
    }

    protected void vAtm(Apfloat volume) {
        this.v_atm = volume;
    }

    protected void nAtm(Apfloat number) {
        this.n_atm = number;
    }

    /**
     * additional parameters to add to toString method
     */
    private String[] additional_string_params;

    private int print_digits;


    protected void addParams(String[] params) {
        int length = params.length;
        additional_string_params = new String [length];
        for (int i = 0; i < length; i++) {
            additional_string_params[i] = params[i];
        }
    }

    protected Apfloat volumeNptApfloat() {
        return this.v_npt;
    }

    protected Apfloat volumeAtmApfloat() {
        return this.v_atm;
    }

    protected Apfloat radiusAtmApfloat() {
        return this.r_atm;
    }

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

    public String fetchName() {
        return this.name;
    }

    public int fetchTemp() { return this.temp; }

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
        np.append(String.format("\nRadius of Atom:        \t%s nm", this.atomRadius()));
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
