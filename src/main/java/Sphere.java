import org.apfloat.*;
import java.util.*;
import com.google.common.primitives.UnsignedInteger;

/**
 * Models a nanoparticle as an icosahedral sphere and calculates atomic packing.
 * <p>
 * This implementation treats the nanoparticle as a regular icosahedron composed of
 * closely packed atomic spheres. Calculations are based on:
 * <ul>
 *   <li>Nanoparticle radius</li>
 *   <li>Atomic radius</li>
 *   <li>Icosahedral geometry relationships</li>
 * </ul>
 *
 * <h2>Key Formulas</h2>
 * <pre>
 * r_npt = l_ico × (√(10 + 2√5)) / 4
 * l_ico = (4 × r_npt) / √(10 + 2√5)
 * v_npt = (15 + 5√5)/12 × l_ico³
 * v_atm = (4/3)π(r_atm)³
 * n_atm = v_npt / v_atm
 * </pre>
 *
 * <h2>References</h2>
 * Implements the geometric model from:
 * Gieseke et al., "SAXS-based structure analysis of ligand-free gold nanoparticles
 * in aqueous solution", Journal of Nanoparticle Research (2016)
 * DOI: 10.1007/s11051-016-3587-7
 */
public class Sphere extends Nanoparticle {

    /**
     * Volume of the icosahedral nanoparticle in cubic nanometers (nm³).
     */
    private Apfloat v_npt;

    /**
     * Edge length of the icosahedron in nanometers (nm).
     */
    private Apfloat l_ico;

    /**
     * Radius of the nanoparticle in nanometers (nm).
     */
    private Apfloat r_npt;

    /**
     * Returns the radius of the nanoparticle as a string with full configured precision.
     *
     * @return string representation of the nanoparticle radius (Å)
     */
    public String nanoParticleRadius() {
        Apfloat fmt_rIco = new Apfloat (
                this.r_npt.toString(),
                super.fetchDisplayDigits()
        );
        return fmt_rIco.toString();
    }

    /**
     * Returns the edge length of the icosahedron with specified significant digits.
     *
     * @return string representation of icosahedron edge length (nm)
     */
    public String lengthOfEdge() {
        Apfloat fmt_lIco = new Apfloat(
                this.l_ico.toString(),
                super.fetchDisplayDigits()
        );
        return fmt_lIco.toString();
    }


    /**
     * Constructs a sphere-modeled nanoparticle and calculates properties.
     *
     * @param r_npt        Nanoparticle radius in nanometers (nm)
     * @param r_atm        Atomic radius in angstroms (Å)
     * @param precision    Calculation precision (significant digits)
     * @param name         Atom name (mmCIF compatible)
     * @param temp         Equilibrium temperature in Kelvin (K)
     * @param print_digits Display precision for string outputs
     */
    public Sphere(
            String r_npt,
            String r_atm,
            UnsignedInteger precision,
            String name,
            int temp,
            int print_digits
    ) {
        // construct nanoparticle instance
        super(
            r_atm,
            precision,
            name,
            temp,
            print_digits
        );

        // Define constants with precision
        Apfloat PI = ApfloatMath.pi(super.fetchPrecision());
        Apfloat TWO = new Apfloat("2", super.fetchPrecision());
        Apfloat THREE = new Apfloat("3", super.fetchPrecision());
        Apfloat FOUR = new Apfloat("4", super.fetchPrecision());
        Apfloat FIVE = new Apfloat("5", super.fetchPrecision());
        Apfloat TEN = new Apfloat("10", super.fetchPrecision());
        Apfloat TWELVE = new Apfloat("12", super.fetchPrecision());
        Apfloat FIFTEEN = new Apfloat("15", super.fetchPrecision());
        Apfloat SIXTEEN = new Apfloat("16", super.fetchPrecision());
        Apfloat SIXTY_FOUR = new Apfloat("64", super.fetchPrecision());

        // Set nanoparticle radius in nm
        this.r_npt = new Apfloat(r_npt, super.fetchPrecision());

        // Calculate icosahedron edge length l_ico from nanoparticle radius r_npt:
        // l_ico = (4 * r_npt) / sqrt(10 + 2 * sqrt(5))
        Apfloat l_ico_1 = FOUR.multiply(this.r_npt);
        Apfloat l_ico_2 = TEN.add(TWO.multiply(ApfloatMath.sqrt(FIVE)));
        this.l_ico = l_ico_1.divide(ApfloatMath.sqrt(l_ico_2));

        // Calculate volume of icosahedron:
        // v_npt = ((15 + 5 * sqrt(5)) / 12) * l_ico^3
        Apfloat v_npt_1 = ApfloatMath.pow(this.l_ico, 3);
        Apfloat v_npt_2 = FIFTEEN.add(FIVE.multiply(ApfloatMath.sqrt(FIVE))).divide(TWELVE);
        super.vNpt(v_npt_1.multiply(v_npt_2));

        // Calculate volume of a single atom:
        // v_atm = (4/3) * π * r_atm^3
        Apfloat v_atm_1 = ApfloatMath.pow(super.radiusAtmApfloat().divide(TEN), 3);
        Apfloat v_atm_2 = FOUR.multiply(PI).divide(THREE);
        super.vAtm(v_atm_1.multiply(v_atm_2));

        // Calculate the number of atoms n_atm using derived formula:
        // n_atm = ⌊((15+5*sqrt(5)) / (16 * π)) * (r_npt / r_atm)^3 * (64 / (10 + 2 * sqrt(5))^(3/2))⌋
        Apfloat n_atm_1 = FIFTEEN.add(FIVE.multiply(ApfloatMath.sqrt(FIVE))).divide(SIXTEEN.multiply(PI));
        Apfloat n_atm_2 = ApfloatMath.pow(this.r_npt.divide(super.radiusAtmApfloat().divide(TEN)), 3);
        Apfloat n_atm_3 = SIXTY_FOUR.divide(ApfloatMath.pow(TEN.add(TWO.multiply(ApfloatMath.sqrt(FIVE))), THREE.divide(TWO)));
        super.nAtm(ApfloatMath.floor(n_atm_1.multiply(n_atm_2).multiply(n_atm_3)));

        // add additional string parameters to display output
        String [] ico_params = new String[] {
            String.format("\nRadius of Nanoparticle:\t%s Å",
                    this.nanoParticleRadius()
            ),
            String.format("\nLength of Nanoparticle:\t%s nm",
                    this.lengthOfEdge()
            )
        };

        super.addParams(ico_params);

        // print status
        System.out.println("\nDone preliminary calculations...");
        System.out.println(super.toString());
    }

    /**
     * Entry point for user interaction.
     * Prompts for nanoparticle radius, atomic radius, and desired precision,
     * then prints out all calculated values.
     *
     * @param args command line arguments (unused)
     */
    public static void main(String[] args) {

        // Prompt user inputs
        Scanner inp = new Scanner(System.in);
        System.out.println("Input temperature of nanoparticle:");
        int t_npt = inp.nextInt();
        System.out.println("Input radius of nanoparticle (nm):");
        String r_npt = inp.nextLine();
        System.out.println("Input name   of atom    :");
        String a_name = inp.nextLine();
        System.out.println("Input radius of atom (Å):");
        String r_atm = inp.nextLine();
        System.out.println("Enter digits of precision:");
        String _prec = inp.nextLine();
        UnsignedInteger precision = UnsignedInteger.valueOf(_prec);

        // Perform calculation and output results
        Sphere result = new Sphere(r_npt, r_atm, precision, a_name, t_npt, 5);
        System.out.println(result);
    }
}
