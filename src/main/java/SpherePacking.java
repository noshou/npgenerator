import org.apfloat.*;
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
public class SpherePacking {

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
    public SpherePacking(String r_npt, String r_atm, UnsignedInteger precision) {

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

        this.precision = precision.intValue();
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

    public String nanoParticleRadius() {
        return this.r_npt.toString();
    }

    public String atomRadius() {
        return this.r_atm.toString();
    }

    public String atomVolume() {
        return this.v_atm.toString();
    }

    public String atomsInNanoParticle() {
        return this.n_atm.toString();
    }

    public String volumeOfNanoparticle() {
        return this.v_ico.toString();
    }

    public String lengthOfEdge() {
        return this.l_ico.toString();
    }

    public int fetchPrecision() {
        return this.precision;
    }

    public String toString() {
        StringBuilder np = new StringBuilder();
        np.append("\n");
        np.append(String.format("Calculation  Precision:\t%s digits", fetchPrecision()));
        np.append(String.format("Atoms  in Nanoparticle:\t%s", this.atomsInNanoParticle()));
        np.append(String.format("Volume of Nanoparticle:\t%snm³", this.volumeOfNanoparticle()));
        np.append(String.format("Radius of Nanoparticle:\t%snm", this.nanoParticleRadius()));
        np.append(String.format("Length of Nanoparticle:\t%snm", this.lengthOfEdge()));
        np.append(String.format("Volume of Atom:        \t%snm³", this.atomVolume()));
        np.append(String.format("Radius of Atom:        \t%snm", this.atomRadius()));
        np.append("\n");
        return np.toString();
    }
}
