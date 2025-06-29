import com.oson.tuple.*;
import org.apfloat.*;
import java.io.IOException;


public class Sphere extends Shape {

    private NpMmcifBuilder file;
    public Sphere(
            String radius,
            LatticeType lattice_type,
            int precision,
            Polyad<Atom> basis,
            String lattice_constant,
            String file_name,
            String structure_name,
            String structure_index
    ) {
        super(
                radius,
                lattice_type,
                precision,
                basis,
                lattice_constant,
                file_name,
                structure_name,
                structure_index
        );
    }

    @Override
    public void build() {

        // get file instance
        try {
            this.file = NpMmcifBuilder.getInstance(this.file_name);
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }

        // initialize writing shape
        try {
            this.file.initShape(super.getThis());  // pass a Shape instance `s` here, not `Shape s`
        } catch (IOException e1) {
            try {
                this.file.abort();
            } catch (IOException abortException) {
                e1.addSuppressed(abortException);
            }
            throw new RuntimeException(e1);
        }

        // write atoms
        int index = 0;
        Apfloat r = new Apfloat(this.radius, this.precision);
        Apfloat a = new Apfloat(this.lattice_constant, this.precision);
        Triad<String> curr = coordinates.getPosition();

        while (curr != null) {
            Apfloat x_frac = new Apfloat(curr.fetch(0), this.precision);
            Apfloat y_frac = new Apfloat(curr.fetch(1), this.precision);
            Apfloat z_frac = new Apfloat(curr.fetch(2), this.precision);
            Apfloat x_cart = x_frac.multiply(a);
            Apfloat y_cart = y_frac.multiply(a);
            Apfloat z_cart = z_frac.multiply(a);

            // noinspection SuspiciousNameCombination
            Apfloat circle = ApfloatMath.pow(x_cart, 2)
                    .add(ApfloatMath.pow(y_cart, 2))
                    .add(ApfloatMath.pow(z_cart, 2));

            // x_cart^2 + y^2 + z^2 ≤ r^2
            if (circle.compareTo(ApfloatMath.pow(r, 2)) <= 0) {
                Atom curr_atom = unit_cell.getLatticePoint(
                        x_frac.doubleValue(),
                        y_frac.doubleValue(),
                        z_frac.doubleValue()
                );
                if (curr_atom != null) {
                    curr_atom.latticePoint(
                            index,
                            new Triad<>(
                                    x_cart.toString(),
                                    y_cart.toString(),
                                    z_cart.toString()
                            )
                    );
                    try {
                        this.file.addAtom(curr_atom);
                    } catch (IOException e2) {
                        try {
                            this.file.abort();
                        } catch (IOException abortException) {
                            e2.addSuppressed(abortException);
                        }
                        throw new RuntimeException(e2);
                    }
                }
            }
        index++;
        curr = coordinates.getPosition();
        }
    }
}