import com.oson.tuple.*;
import org.apfloat.*;

/**
 * Represents a 3D axis-aligned bounding box that encloses a sphere centered at a given point
 * with a specified radius. The class computes the eight corner vertices and six face centroids
 * of the bounding box using high-precision arithmetic via Apfloat.
 *
 * <p>This class assumes the bounding box is a cube aligned with the coordinate axes, where
 * each face is offset from the center by a given radius in the ±x, ±y, ±z directions.
 *
 * <p>The inputs are provided as strings and wrapped in generic tuple types (Triad, Hexad, Octad).
 */
public class BoundingBox {

    /** The eight corner vertices of the bounding box. */
    private Octad<Tuple<String>> vertices;

    /** The six face-centered points (centroids) of the bounding box. */
    private Hexad<Tuple<String>> centroids;

    /** The center of the bounding box (which matches the center of the enclosed sphere). */
    private Triad<String> center;

    /**
     * Constructs a bounding box around a sphere defined by its center and radius.
     *
     * @param centroid A Triad of strings representing the x, y, and z coordinates of the sphere's center.
     * @param radius A string representation of the sphere's radius.
     * @param precision The precision used for all Apfloat calculations (number of decimal digits).
     */
    public BoundingBox(
            Triad<String> centroid,
            String radius,
            int precision) {

        this.center = new Triad<>(
                centroid.fetch(0),
                centroid.fetch(1),
                centroid.fetch(2)
        );

        // Convert input strings to Apfloat values
        Apfloat x = new Apfloat(this.center.fetch(0), precision);
        Apfloat y = new Apfloat(this.center.fetch(1), precision);
        Apfloat z = new Apfloat(this.center.fetch(2), precision);
        Apfloat r = new Apfloat(radius, precision);

        // Compute the six face centroids by offsetting one axis at a time
        Triad<String> cx1 = new Triad<>(x.add(r).toString(), y.toString(), z.toString());
        Triad<String> cx2 = new Triad<>(x.subtract(r).toString(), y.toString(), z.toString());
        Triad<String> cy1 = new Triad<>(x.toString(), y.add(r).toString(), z.toString());
        Triad<String> cy2 = new Triad<>(x.toString(), y.subtract(r).toString(), z.toString());
        Triad<String> cz1 = new Triad<>(x.toString(), y.toString(), z.add(r).toString());
        Triad<String> cz2 = new Triad<>(x.toString(), y.toString(), z.subtract(r).toString());
        this.centroids = new Hexad<>(cx1, cx2, cy1, cy2, cz1, cz2);

        // Compute all eight bounding box vertices (corners) from combinations of ±r offsets
        Triad<String> v_x1_y1_z1 = new Triad<>(
                x.add(r).toString(),
                y.add(r).toString(),
                z.add(r).toString()
        );
        Triad<String> v_x1_y1_z2 = new Triad<>(
                x.add(r).toString(),
                y.add(r).toString(),
                z.subtract(r).toString()
        );
        Triad<String> v_x1_y2_z2 = new Triad<>(
                x.add(r).toString(),
                y.subtract(r).toString(),
                z.subtract(r).toString()
        );
        Triad<String> v_x2_y2_z2 = new Triad<> (
                x.subtract(r).toString(),
                y.subtract(r).toString(),
                z.subtract(r).toString()
        );
        Triad<String> v_x2_y2_z1 = new Triad<> (
                x.subtract(r).toString(),
                y.subtract(r).toString(),
                z.add(r).toString()
        );
        Triad<String> v_x1_y2_z1 = new Triad<> (
                x.add(r).toString(),
                y.subtract(r).toString(),
                z.add(r).toString()
        );
        Triad<String> v_x2_y1_z2 = new Triad<> (
                x.subtract(r).toString(),
                y.add(r).toString(),
                z.subtract(r).toString()
        );
        Triad<String> v_x2_y1_z1 = new Triad<> (
                x.subtract(r).toString(),
                y.add(r).toString(),
                z.add(r).toString()
        );
        this.vertices = new Octad<Tuple<String>> (
                v_x1_y1_z1,
                v_x1_y1_z2,
                v_x1_y2_z2,
                v_x2_y2_z2,
                v_x2_y2_z1,
                v_x1_y2_z1,
                v_x2_y1_z2,
                v_x2_y1_z1
        );
    }

    /**
     * Returns the eight vertices of the bounding box.
     *
     * @return An {@code Octad} containing eight {@code Triad<String>} instances representing the corners.
     */
    public Octad<Tuple<String>> fetchVertices() {
        return this.vertices;
    }

    /**
     * Returns the six face centroids of the bounding box.
     *
     * @return A {@code Hexad} of {@code Triad<String>} instances representing the face centers.
     */
    public Hexad<Tuple<String>> fetchCentroids() {
        return this.centroids;
    }

    /**
     * Returns the center of the bounding box (same as the center of the inner sphere).
     *
     * @return A {@code Triad<String>} representing the center coordinates.
     */
    public Triad<String> fetchCenter() {
        return this.center;
    }

    /**
     * Returns a string representation of the bounding box, including its vertices,
     * face centroids, and center point.
     *
     * @return A formatted string listing the key components of the bounding box.
     */
    @Override
    public String toString() {
        StringBuilder bbox = new StringBuilder();
        bbox.append("Vertices: \t").append(this.fetchVertices());
        bbox.append("\nCentroids:\t").append(this.fetchCentroids());
        bbox.append("\nCenter:   \t").append(this.fetchCenter());
        return bbox.toString();
    }
}
