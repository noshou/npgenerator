package Lattice;

import Atom.Atom;
import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;

/**
 * Abstract representation of a crystallographic unit cell.
 * <p>
 * Encapsulates geometric parameters (unit cell lengths and angles),
 * symmetry space group, basis atoms (each with fractional positions),
 * and lattice type. Concrete subclasses define specific logic for retrieving
 * atomic data at lattice coordinates.
 * </p>
 * <p>
 * This base class supports extension for Bravais systems such as FCC, BCC,
 * HCP, etc., and can be used to construct periodic or finite crystalline structures.
 * </p>
 */
public abstract class BravaisUnitCell extends  UnitCell {

    /** Unit cell lengths along the a, b, and c axes (in Ångströms). */
    protected final @NotNull Apfloat a, b, c;

    /** Interaxial angles α (between b and c), β (between a and c),
     * γ (between a and b) in degrees. */
    protected final @NotNull Apfloat alpha, beta, gamma;

    /**
     * Constructs a general crystallographic unit cell with required geometric
     * and atomic properties.
     *
     * @param a            unit cell edge length a (in Å); must not be null
     * @param b            unit cell edge length b (in Å); must not be null
     * @param c            unit cell edge length c (in Å); must not be null
     * @param alpha        interaxial angle α (in degrees); must not be null
     * @param beta         interaxial angle β (in degrees); must not be null
     * @param gamma        interaxial angle γ (in degrees); must not be null
     * @param space_group  Hermann–Mauguin space group label; must not be null
     * @param basis        list of basis atoms in fractional coordinates; must not be null
     * @param lattice_type type of the lattice; must not be null
     * @param precision    digits of precision
     */
    public BravaisUnitCell(
            @NotNull Apfloat a,
            @NotNull Apfloat b,
            @NotNull Apfloat c,
            @NotNull Apfloat alpha,
            @NotNull Apfloat beta,
            @NotNull Apfloat gamma,
            @NotNull String space_group,
            @NotNull Polyad<Atom> basis,
            @NotNull LatticeType lattice_type,
            int precision
    ) {
        super(precision, lattice_type, basis, space_group);
        this.a = a;
        this.b = b;
        this.c = c;
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
    }
    /** @return interaxial angle α (degrees, between b and c) */
    @Contract(pure = true)
    private @NotNull Dyad<String> getAlpha() {
        return new Dyad<>("alpha", this.alpha.toString());
    }
    /** @return interaxial angle β (degrees, between a and c) */
    @Contract(pure = true)
    private @NotNull Dyad<String> getBeta() {
        return new Dyad<>("beta", this.beta.toString());
    }
    /** @return interaxial angle γ (degrees, between a and b) */
    @Contract(pure = true)
    private @NotNull Dyad<String> getGamma() {
        return new Dyad<>("gamma", this.gamma.toString());
    }

    /**
     * Retrieves all interaxial angles (<code>a</code>, <code>β</code>,
     * and <code> γ</code>) as a typed {@link Polyad}
     * of {@link Tuple} instances. Each tuple contains a label and its associated value.
     *
     * <p>
     * This method uses a raw array cast to generic type,
     * but is safe under internal constraints.
     * </p>
     *
     * @return a non-null {@link Polyad} of {@link Tuple}&lt;String&gt;
     * representing the interaxial unit cell angles in degrees.
     */
    @Override
    @SuppressWarnings("unchecked")
    public @NotNull Polyad<Tuple<String>> getCellAngles() {
        return new Polyad<Tuple<String>>(
                new Tuple[]{
                        this.getAlpha(),
                        this.getBeta(),
                        this.getGamma()
                }
        );
    }

    /**
     * Retrieves the unit cell edge length <code>a</code> in Ångströms.
     *
     * @return a {@link Dyad} containing the label "a"
     * and the corresponding edge length value
     */
    @Contract(pure = true)
    private @NotNull Dyad<String> getA() {
        return new Dyad<>("a", this.a.toString());
    }

    /**
     * Retrieves the unit cell edge length <code>b</code> in Ångströms.
     *
     * @return a {@link Dyad} containing the label "b" and
     * the corresponding edge length value
     */
    @Contract(pure = true)
    private @NotNull Dyad<String> getB() {
        return new Dyad<>("b", this.b.toString());
    }

    /**
     * Retrieves the unit cell edge length <code>c</code> in Ångströms.
     *
     * @return a {@link Dyad} containing the label "c" and
     * the corresponding edge length value
     */
    @Contract(pure = true)
    private @NotNull Dyad<String> getC() {
        return new Dyad<>("c", this.c.toString());
    }

    /**
     * Retrieves all unit cell lengths (<code>a</code>, <code>b</code>,
     * and <code>c</code>) as a typed {@link Polyad}
     * of {@link Tuple} instances. Each tuple contains a label and its associated value.
     *
     * <p>
     * This method uses a raw array cast to generic type, but is safe under
     * internal constraints.
     * </p>
     *
     * @return a non-null {@link Polyad} of {@link Tuple}&lt;String&gt;
     * representing the unit cell lengths
     */
    @Override
    @SuppressWarnings("unchecked")
    public @NotNull Polyad<Tuple<String>> getCellLengths() {
        return new Polyad<Tuple<String>>(
                new Tuple[]{
                        this.getA(),
                        this.getB(),
                        this.getC()
                }
        );
    }
}
