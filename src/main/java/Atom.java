import com.oson.tuple.*;

/**
 * Represents a single atom with its identity, position, bounding box, and related properties.
 * <p>
 * This class assumes all inputs are valid and properly formatted.
 * The radius is stored as a String to accommodate Apfloat objects handled externally.
 * </p>
 */
class Atom {

    /** The chemical element symbol of the atom (e.g., "Au") */
    private final String element;

    /** The atomic radius as a string (Angstroms) */
    private final String radius;

    /** The unique index or identifier of this atom */
    private final int index;

    /** The Cartesian coordinates of the atom as a triad (x, y, z) */
    private final Triad<Double> coordinates;

    /** The B-factor (temperature factor) associated with the atom */
    private final double b_factor;

    /**
     * The bounding box of the atom represented by eight corner points,
     * stored as an Octad of tuples where each tuple is a (x, y, z) coordinate.
     */
    private final Octad<Tuple<Double>> bounding_box;

    private final String formal_charge;

    private final String b_iso;

    /**
     * Constructs an Atom instance with specified properties.
     *
     * @param name the name of the atom (e.g., "Au1")
     * @param element the element symbol (e.g., "Au")
     * @param radius the radius of the atom as a String for Apfloat compatibility
     * @param index unique index of the atom
     * @param x_coord x-coordinate of the atom in Cartesian space
     * @param y_coord y-coordinate of the atom in Cartesian space
     * @param z_coord z-coordinate of the atom in Cartesian space
     * @param b_factor the B-factor (temperature factor) of the atom
     * @param bbox_vert_1 the first vertex (corner) of the bounding box as a Triad of Doubles (x, y, z)
     * @param bbox_vert_2 the second vertex (corner) of the bounding box as a Triad of Doubles (x, y, z)
     * @param bbox_vert_3 the third vertex (corner) of the bounding box as a Triad of Doubles (x, y, z)
     * @param bbox_vert_4 the fourth vertex (corner) of the bounding box as a Triad of Doubles (x, y, z)
     * @param bbox_vert_5 the fifth vertex (corner) of the bounding box as a Triad of Doubles (x, y, z)
     * @param bbox_vert_6 the sixth vertex (corner) of the bounding box as a Triad of Doubles (x, y, z)
     * @param bbox_vert_7 the seventh vertex (corner) of the bounding box as a Triad of Doubles (x, y, z)
     * @param bbox_vert_8 the eighth vertex (corner) of the bounding box as a Triad of Doubles (x, y, z)
     */
    public Atom(
            String name,
            String element,
            String radius,
            int formal_charge,
            int index,
            double x_coord,
            double y_coord,
            double z_coord,
            double b_factor,
            Triad<Double> bbox_vert_1,
            Triad<Double> bbox_vert_2,
            Triad<Double> bbox_vert_3,
            Triad<Double> bbox_vert_4,
            Triad<Double> bbox_vert_5,
            Triad<Double> bbox_vert_6,
            Triad<Double> bbox_vert_7,
            Triad<Double> bbox_vert_8
    ) {
        this.element = element;
        this.radius = radius;

        if (formal_charge > 0) {
            this.formal_charge = new String(String.format("%d+", formal_charge));
        } else if (formal_charge < 0) {
            this.formal_charge = new String(String.format("%d-", formal_charge));
        } else {
            this.formal_charge = "0";
        }

        this.index = index;
        this.coordinates = new Triad<>(x_coord, y_coord, z_coord);
        this.b_factor = b_factor;
        this.bounding_box = new Octad<>(
                bbox_vert_1,
                bbox_vert_2,
                bbox_vert_3,
                bbox_vert_4,
                bbox_vert_5,
                bbox_vert_6,
                bbox_vert_7,
                bbox_vert_8
        );

        // calculate b-factor
        this.b_iso = "-1";
    }


    public String fetchFormalCharge() {
        return this.formal_charge;
    }

    /**
     * Returns the element symbol of the atom.
     *
     * @return the element symbol
     */
    public String fetchElement() {
        return this.element;
    }

    /**
     * Returns the atomic radius as a String.
     *
     * @return the radius string (Apfloat compatible)
     */
    public String fetchRadius() {
        return this.radius;
    }

    /**
     * Returns the unique index of the atom
     *
     * @return the atom index
     */
    public int fetchIndex() {
        return this.index;
    }

    /**
     * Returns the Cartesian coordinates of the atom.
     *
     * @return a {@code Triad<Double>} representing (x, y, z)
     */
    public Triad<Double> fetchCoordinates() {
        return this.coordinates;
    }

    /**
     * Returns the B-factor (temperature factor) of the atom.
     *
     * @return the B-factor value in Å²
     */
    public String bIso() {
        return this.b_iso;
    }

    /**
     * Returns the bounding box of the atom as a dyad of tuples representing the two corners.
     *
     * @return the bounding box as a {@code Dyad<Tuple<Double>>}
     */
    public Octad<Tuple<Double>> boundingBox() {
        return this.bounding_box;
    }
}
