package FileWriter;

import org.apfloat.Apfloat;
import org.jetbrains.annotations.*;
import java.io.IOException;

/**
 * A writer for debugging atomic coordinates, exporting them in CSV format.
 * <p>
 * The output contains both fractional and Cartesian coordinates, along with
 * an occupancy flag. This is primarily used for diagnostics or visualization.
 */
public class CoordsDebugWriter extends FileWriter {

    /**
     * Constructs a {@code FileWriter.CoordsDebugWriter} that writes to a CSV file.
     *
     * @param file_name The base name of the output file (without extension).
     * @throws IOException If the underlying file cannot be created or opened for writing.
     */
    public CoordsDebugWriter(@NotNull String file_name) throws IOException {
        super(file_name, ".csv");
    }

    /**
     * Initializes the debug coordinate writer with a fixed CSV header.
     * <p>
     * This method must be called with a {@code null} initializer. Any non-null
     * argument will result in an exception.
     *
     * @param initializer Must be {@code null}.
     * @throws IOException              If writing the header fails.
     * @throws IllegalArgumentException If {@code initializer} is not {@code null}.
     */
    @Override
    @Contract("!null -> fail")  // if param is not null, throw exception
    public void init(@Nullable Object initializer) throws IOException {
        if (initializer != null) {
            throw new IllegalArgumentException("initializer must be null!");
        }
        writer.write("x_frac,y_frac,z_frac,x_cart,y_cart,z_cart,is_occupied\n");
    }

    /**
     * Appends a line of coordinate data to the debug CSV file.
     *
     * @param x_frac      The fractional x-coordinate.
     * @param y_frac      The fractional y-coordinate.
     * @param z_frac      The fractional z-coordinate.
     * @param x_cart      The Cartesian x-coordinate.
     * @param y_cart      The Cartesian y-coordinate.
     * @param z_cart      The Cartesian z-coordinate.
     * @param is_occupied {@code true} if the site is occupied, {@code false} otherwise.
     * @throws IOException If writing to the file fails.
     */
    @Contract(pure = false)
    public void addCoordinate(
            @NotNull Apfloat x_frac,
            @NotNull Apfloat y_frac,
            @NotNull Apfloat z_frac,
            @NotNull Apfloat x_cart,
            @NotNull Apfloat y_cart,
            @NotNull Apfloat z_cart,
            boolean is_occupied
    ) throws IOException {
        String line = String.format(
                "%s,%s,%s,%s,%s,%s,%b\n",
                x_frac.toString(),
                y_frac.toString(),
                z_frac.toString(),
                x_cart.toString(),
                y_cart.toString(),
                z_cart.toString(),
                is_occupied
        );
        writer.write(line);
    }
}
