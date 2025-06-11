import org.apfloat.*;
import java.util.*;
import com.google.common.primitives.UnsignedInteger;

/**
 * Calculates optimally filled sphere of atoms, assuming an icosahedral geometry.
 * <p>
 * This model treats the nanoparticle as a regular icosahedron composed of
 * closely packed atoms, each modeled as a sphere. The computation is based
 * on the radius of the nanoparticle and the radius of the atom composing the nanoparticle.
 * </p>
 *
 * <h2>Mathematical Formulation</h2>
 * <p>
 * Let:
 * <ul>
 *   <li>{@code r_npt} = radius of the nanoparticle (nanometers, nm)</li>
 *   <li>{@code r_atm} = radius of a single atom (angstroms, converted to nanometers)</li>
 *   <li>{@code n_atm} = estimated number of atoms in the nanoparticle (unitless count)</li>
 *   <li>{@code l_ico} = edge length of the icosahedron (nanometers)</li>
 *   <li>{@code v_ico} = volume of the icosahedron (cubic nanometers, nm³)</li>
 *   <li>{@code v_atm} = volume of a single atom (cubic nanometers, nm³)</li>
 * </ul>
 *
 * <p>From icosahedral geometry, we derive the following:</p>
 *
 * <pre>{@code
 * r_npt = l_ico * (sqrt(10 + 2 * sqrt(5)) / 4)
 * => l_ico = (4 / sqrt(10 + 2 * sqrt(5))) * r_npt
 *
 * v_ico = ((15 + 5 * sqrt(5)) / 12) * l_ico^3
 * => v_ico = ((15 + 5 * sqrt(5)) / 12) * ((4 / sqrt(10 + 2 * sqrt(5))) * r_npt)^3
 *
 * v_atm = (4 / 3) * π * r_atm^3
 *
 * n_atm = v_ico / v_atm
 * => n_atm = ((15 + 5 * sqrt(5)) / (16 * π)) * (r_npt / r_atm)^3 * (64 / (10 + 2 * sqrt(5))^(3/2))
 * }</pre>
 *
 * <p>
 * This yields a closed-form expression to estimate the atom count based on known physical constants.
 * </p>
 *
 * <p><b>Reference:</b> The model and formulae are adapted from the study:
 * <a href="https://link-springer-com.cyber.usask.ca/article/10.1007/s11051-016-3587-7">
 * Gieseke et al., "SAXS-based structure analysis of ligand-free gold nanoparticles in aqueous solution",
 * Journal of Nanoparticle Research (2016)
 * </a>
 * </p>
 *
 * @author [Your Name]
 */
public class SpherePacked {

    private Apfloat r_npt;
    private Apfloat r_atm;
    private Apfloat v_atm;
    private Apfloat n_atm;
    private Apfloat v_ico;
    private Apfloat l_ico;
    private int precision;

    /**
     * Constructs a {@code GoldSphereCalibration} instance to model an icosahedral
     * gold nanoparticle and estimate the number of gold atoms it contains.
     *
     * @param r_npt        The radius of the nanoparticle in nanometers (nm). Must be a string representation of a decimal number.
     * @param r_atm       The atomic radius of a gold atom in angstroms (Å). Must be a string representation of a decimal number.
     *                   This value will be internally converted to nanometers (1 Å = 0.1 nm).
     * @param precision  The number of significant digits to retain in all Apfloat computations.
     *                   Higher values ensure more accurate mathematical operations but are computationally heavier.
     */
    public SpherePacked(String r_npt, String r_atm, UnsignedInteger precision) {
        // get precision
        this.precision = precision.intValue();

        // constants
        Apfloat PI = ApfloatMath.pi(this.precision);
        Apfloat TWO = new Apfloat("2", this.precision);
        Apfloat THREE = new Apfloat("3", this.precision);
        Apfloat FOUR = new Apfloat("4", this.precision);
        Apfloat FIVE = new Apfloat("5", this.precision);
        Apfloat TEN = new Apfloat("10", this.precision);
        Apfloat TWELVE = new Apfloat("12", this.precision);
        Apfloat FIFTEEN = new Apfloat("15", this.precision);
        Apfloat SIXTEEN = new Apfloat("16", this.precision);
        Apfloat SIXTY_FOUR = new Apfloat("64", this.precision);

        this.r_atm = new Apfloat(r_atm, this.precision).divide(TEN);
        this.r_npt = new Apfloat(r_npt, this.precision);

        // calculate the length of the icosahedral shell, l_ico, given a radius r_npt
        // l_ico = (4*r_npt)/sqrt(10+2*sqrt(5))
        Apfloat l_ico_1 = FOUR.multiply(this.r_npt);
        Apfloat l_ico_2 = TEN.add(TWO.multiply(ApfloatMath.sqrt(FIVE)));
        this.l_ico = l_ico_1.divide(l_ico_2);

        // calculate volume of icosahedral shell
        // v_ico = (lIco^3)*((15+5*sqrt(5))/12)
        Apfloat v_ico_1 = ApfloatMath.pow(this.l_ico, 3);
        Apfloat v_ico_2 = FIFTEEN.add(FIVE.multiply(ApfloatMath.sqrt(FIVE))).divide(TWELVE);
        this.v_ico = v_ico_1.multiply(v_ico_2);

        // calculate volume of gold atom (given radius r_atm), v_atm
        // v_atm = (r_atm^3)*((4*pi)/3)
        Apfloat v_atm_1 = ApfloatMath.pow(this.r_atm, 3);
        Apfloat v_atm_2 = (FOUR).multiply(PI).divide(THREE);
        this.v_atm = v_atm_1.multiply(v_atm_2);

        // calculate number of gold atoms, n_atm
        // n_atm = ((15+5*sqrt(5)/(16*pi))*((r_npt/r_atm)^3)*(64/(10+2*sqrt(5)^(3/2)))
        Apfloat n_atm_1 = FIFTEEN.add(FIVE.multiply(ApfloatMath.sqrt(FIVE))).divide(SIXTEEN.multiply(PI));
        Apfloat n_atm_2 = ApfloatMath.pow(this.r_npt.divide(this.r_atm), 3);
        Apfloat n_atm_3 = SIXTY_FOUR.divide(ApfloatMath.pow(TEN.add(TWO.multiply(ApfloatMath.sqrt(FIVE))), THREE.divide(TWO)));
        this.n_atm = n_atm_1.multiply(n_atm_2).multiply(n_atm_3);
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
     * Returns the radius of the nanoparticle using the full configured precision.
     *
     * @return string representation of the nanoparticle radius (in nanometers)
     */
    public String nanoParticleRadius() {
        return nanoParticleRadius(this.fetchPrecision());
    }

    /**
     * Returns the radius of the nanoparticle with specified number of significant digits.
     *
     * @param digits desired number of significant digits (ignored if non-positive or equal to current precision)
     * @return string representation of the nanoparticle radius (in nanometers)
     */
    public String nanoParticleRadius(int digits) {
        if (digits == this.fetchPrecision() || digits <= 0) {
            return this.r_npt.toString();
        }
        Apfloat fmt_rNpt = new Apfloat(this.r_npt.toString(), digits);
        return fmt_rNpt.toString();
    }

    /**
     * Returns the radius of a single atom using the full configured precision.
     *
     * @return string representation of the atomic radius (in nanometers)
     */
    public String atomRadius() {
        return atomRadius(this.fetchPrecision());
    }

    /**
     * Returns the radius of a single atom with specified number of significant digits.
     *
     * @param digits desired number of significant digits (ignored if non-positive or equal to current precision)
     * @return string representation of the atomic radius (in nanometers)
     */
    public String atomRadius(int digits) {
        if (digits == this.fetchPrecision() || digits <= 0) {
            return this.r_atm.toString();
        }
        Apfloat fmt_rAtm = new Apfloat(this.r_atm.toString(), digits);
        return fmt_rAtm.toString();
    }

    /**
     * Returns the volume (nm³) of a single atom using the full configured precision.
     *
     * @return string representation of the atomic volume (in nm³)
     */
    public String atomVolume() {
        return atomVolume(this.fetchPrecision());
    }

    /**
     * Returns the volume (nm³) of a single atom with specified number of significant digits.
     *
     * @param digits desired number of significant digits (ignored if non-positive or equal to current precision)
     * @return string representation of the atomic volume (in nm³)
     */
    public String atomVolume(int digits) {
        if (digits == this.fetchPrecision() || digits <= 0) {
            return this.v_atm.toString();
        }
        Apfloat fmt_vAtm = new Apfloat(this.v_atm.toString(), digits);
        return fmt_vAtm.toString();
    }

    /**
     * Returns the number of atoms estimated in the nanoparticle using the full configured precision.
     *
     * @return string representation of atom count
     */
    public String atomsInNanoParticle() {
        return atomsInNanoParticle(this.fetchPrecision());
    }

    /**
     * Returns the number of atoms estimated in the nanoparticle with specified number of significant digits.
     *
     * @param digits desired number of significant digits (ignored if non-positive or equal to current precision)
     * @return string representation of atom count
     */
    public String atomsInNanoParticle(int digits) {
        if (digits == this.fetchPrecision() || digits <= 0) {
            return this.n_atm.toString();
        }
        Apfloat fmt_nAtm = new Apfloat(this.n_atm.toString(), digits);
        return fmt_nAtm.toString();
    }

    /**
     * Returns the computed volume (nm³) of the nanoparticle using the full configured precision.
     *
     * @return string representation of nanoparticle volume (in nm³)
     */
    public String volumeOfNanoparticle() {
        return volumeOfNanoparticle(this.fetchPrecision());
    }

    /**
     * Returns the computed volume (nm³) of the nanoparticle with specified number of significant digits.
     *
     * @param digits desired number of significant digits (ignored if non-positive or equal to current precision)
     * @return string representation of nanoparticle volume (in nm³)
     */
    public String volumeOfNanoparticle(int digits) {
        if (digits == this.fetchPrecision() || digits <= 0) {
            return this.v_ico.toString();
        }
        Apfloat fmt_vIco = new Apfloat(this.v_ico.toString(), digits);
        return fmt_vIco.toString();
    }

    /**
     * Returns the edge length (nm) of the icosahedron using the full configured precision.
     *
     * @return string representation of the edge length (in nanometers)
     */
    public String lengthOfEdge() {
        return lengthOfEdge(this.fetchPrecision());
    }

    /**
     * Returns the edge length (nm) of the icosahedron with specified number of significant digits.
     *
     * @param digits desired number of significant digits (ignored if non-positive or equal to current precision)
     * @return string representation of the edge length (in nanometers)
     */
    public String lengthOfEdge(int digits) {
        if (digits == this.fetchPrecision() || digits <= 0) {
            return this.l_ico.toString();
        }
        Apfloat fmt_lIco = new Apfloat(this.l_ico.toString(), digits);
        return fmt_lIco.toString();
    }

    public String toString() {
        StringBuilder np = new StringBuilder();
        np.append("\n");
        np.append(String.format("\nCalculation  Precision:\t%s digits", fetchPrecision()));
        np.append(String.format("\nAtoms  in Nanoparticle:\t%s", this.atomsInNanoParticle(5)));
        np.append(String.format("\nVolume of Nanoparticle:\t%snm³", this.volumeOfNanoparticle(5)));
        np.append(String.format("\nRadius of Nanoparticle:\t%snm", this.nanoParticleRadius(5)));
        np.append(String.format("\nLength of Nanoparticle:\t%snm", this.lengthOfEdge(5)));
        np.append(String.format("\nVolume of Atom:        \t%snm³", this.atomVolume(5)));
        np.append(String.format("\nRadius of Atom:        \t%snm", this.atomRadius(5)));
        np.append("\n");
        return np.toString();
    }

    public static void main(String[] args) {

        // get user inputs
        Scanner inp = new Scanner(System.in);
        System.out.println("Input radius of nanoparticle (nm):");
        String r_npt = inp.nextLine();
        System.out.println("Input radius of atom (Å):");
        String r_atm = inp.nextLine();
        System.out.println("Enter digits of precision:");
        String _prec = inp.nextLine();
        UnsignedInteger precision = UnsignedInteger.valueOf(_prec);
        System.out.println(precision);
        // output
        SpherePacked result = new SpherePacked(r_npt, r_atm, precision);
        System.out.println(result);
    }
}
