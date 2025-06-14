import org.apfloat.*;
import java.util.*;
import com.google.common.primitives.UnsignedInteger;
import org.rcsb.cif.model.binary.BinaryFile;

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
 *   <li>{@code v_npt} = volume of the icosahedron (cubic nanometers, nm³)</li>
 *   <li>{@code v_atm} = volume of a single atom (cubic nanometers, nm³)</li>
 * </ul>
 *
 * <p>From icosahedral geometry, we derive the following:</p>
 *
 * <pre>{@code
 * r_npt = l_ico * (sqrt(10 + 2 * sqrt(5)) / 4)
 * => l_ico = (4 / sqrt(10 + 2 * sqrt(5))) * r_npt
 *
 * v_npt = ((15 + 5 * sqrt(5)) / 12) * l_ico^3
 * => v_npt = ((15 + 5 * sqrt(5)) / 12) * ((4 / sqrt(10 + 2 * sqrt(5))) * r_npt)^3
 *
 * v_atm = (4 / 3) * π * r_atm^3
 *
 * n_atm = v_npt / v_atm
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
     * @return string representation of the nanoparticle radius (nanometers)
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
     * Constructs a {@code SpherePacked} instance modeling a nanoparticle and estimates
     * the number of atoms it contains.
     *
     * @param r_npt         Radius of the nanoparticle in nanometers (nm), passed as string.
     * @param r_atm         Radius of a single atom in angstroms (Å), passed as string.
     *                      This value is internally converted to nanometers.
     * @param precision     Number of significant digits for Apfloat calculations.
     * @param name          Name of atom (must be mmcif compatible)
     * @param temp          Temperature of atom (K)
     * @param print_digits  Number of digits to print
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

        // print status
        System.out.println("\nDone preliminary calculations...");
        System.out.println(super.toString());

//        // calculate atoms, process mmCIF file
//        Atom atom;
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
