/**
 * REFERENCE: https://lampz.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
 */
import java.math.*;
import java.util.*;
import com.oson.tuple.*;
import org.apfloat.*;
import com.google.common.primitives.UnsignedInteger;

/**
 * Represents a form factor calculated from coefficient parameters.
 * <p>
 * The form factor is calculated using the formula:
 * <pre>
 * f(|G|) = Σ (i=1 to 4) a_i * exp(-b_i * (q / 4π)^2) + c
 * </pre>
 * where {@code a} and {@code b} are {@link Tetrad}s of coefficients,
 * and {@code c} is a {@link Monad} constant coefficient.
 * </p>
 */
public class FormFactor {

    /**
     * Name identifier for this {@link FormFactor} instance.
     */
    private String name;

    /**
     * The 'a' coefficients for the form factor calculation,
     * stored as a {@link Tetrad} of {@link BigDecimal}.
     */
    private Tetrad<Apfloat> a;

    /**
     * The 'b' coefficients for the form factor calculation,
     * stored as a {@link Tetrad} of {@link BigDecimal}.
     */
    private Tetrad<Apfloat> b;

    /**
     * The constant 'c' coefficient for the form factor,
     * stored as a {@link Monad} of {@link BigDecimal}.
     */
    private Monad<Apfloat> c;

    private int percision;

    /**
     * Constructs a {@link FormFactor} with specified coefficients.
     *
     * @param name A name identifier for this {@link FormFactor}.
     * @param a {@link Tetrad} of {@link String} representing the 'a' coefficients.
     * @param b {@link Tetrad} of {@link String} representing the 'b' coefficients.
     * @param c {@link Monad} of {@link String} representing the constant 'c' coefficient.
     */
    public FormFactor(
            String name,
            Tetrad<String> a,
            Tetrad<String> b,
            Monad<String> c,
            UnsignedInteger percision
    ) {
        this.name = name;
        this.percision = percision.intValue();

        // create coeff tetrads
        this.a = new Tetrad<Apfloat>(
            new Apfloat(a.fetch(0), percision.intValue()),
            new Apfloat(a.fetch(1), percision.intValue()),
            new Apfloat(a.fetch(2), percision.intValue()),
            new Apfloat(a.fetch(3), percision.intValue())
        );
        this.b = new Tetrad<Apfloat>(
            new Apfloat(b.fetch(0), percision.intValue()),
            new Apfloat(b.fetch(1), percision.intValue()),
            new Apfloat(b.fetch(2), percision.intValue()),
            new Apfloat(b.fetch(3), percision.intValue())
        );
        this.c = new Monad<Apfloat>(new Apfloat(c.fetch(0), percision.intValue()));
    }

    public int getPercision() {
        return this.percision;
    }
    /**
     * Computes the form factor value for a given scattering vector magnitude {@code q}.
     *
     * @param q The scattering vector magnitude in Å⁻¹. Must be between 0 and 25.
     * @return The computed form factor value.
     * @throws IllegalArgumentException if {@code q} is outside the valid range [0, 25].
     */
    public String getFormFactor(String q) {
        Apfloat q = new Apfloat(q, this.getPercision());
        Apfloat twenty_five = new Apfloat("25", this.getPercision());

        if (q.compareTo(Apfloat.ZERO) < 0 || q.compareTo(twenty_five) > 0) {
            throw new IllegalArgumentException("q must be between 0 and 25 Å⁻¹");
        }

        Apfloat pi = ApfloatMath.pi(new Apfloat(1, 1000)); // get π to 1000 digits
        Apfloat q_fourpi_sqr = q.divide(
                new Apfloat("4", 1000).multiply(pi)
        ).pow(2);
        Apfloat form_factor = Apfloat("0", this.getPercision());

        for (int i = 0; i < 4; i++) {
            Apfloat a_i = this.a.fetch(i);
            Apfloat b_i = this.b.fetch(i);
            Apfloat f_i = a_i.multiply(ApfloatMath.exp(b_i.negate().multiply(q_fourpi_sqr)));
            form_factor = form_factor.add(f_i);
        }
        form_factor = form_factor.add(this.c.fetch(0));

        return form_factor.toString();
    }
}

