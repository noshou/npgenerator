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

    /**
     * The bounding box of the atom represented by eight corner points,
     * six face centroids, and the center point
     */
    private final BoundingBox bbox;

    private final String formal_charge;

    /**
     * Constructs an Atom instance with specified properties.
     *
     * @param name      the name of the atom (e.g., "Au1")
     * @param element   the element symbol (e.g., "Au")
     * @param radius    the radius of the atom as a String for Apfloat compatibility
     * @param index     unique index of the atom
     * @param charge    formal charge of atom
     * @param centroid  the center point of the atom (x,y,z) which defines the bounding box
     */
    public Atom(
            String name,
            String element,
            String radius,
            int formal_charge,
            int index,
            double charge,
            Triad<Tuple<String>> centroid
    ) {
        this.element = element;
        this.radius = radius;

        if (charge > 0) {
            this.formal_charge = new String(String.format("%d+", formal_charge));
        } else if (formal_charge < 0) {
            this.formal_charge = new String(String.format("%d-", formal_charge));
        } else {
            this.formal_charge = "0";
        }

        this.index = index;

        this.bbox = new BoundingBox(centroid, radius);
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
     * Returns the bounding box (internal use only)
     *
     * @return the bounding box
     */
    private BoundingBox boundingBox() {
        return this.bbox;
    }

    /**
     * Returns the Cartesian coordinates of the atom (center of bounding box)
     *
     * @return a {@code Triad<Double>} representing (x, y, z)
     */
    public Triad<Tuple<String>> fetchCoordinates() {
        return this.boundingBox().fetchCenter();
    }

    /*
    * Returns verticies of bounding box
    */
    public Octad<Tuple<String>> fetchBboxVertices() {
        return this.boundingBox().fetchVertices();
    }

    public Hexad<Tuple<String>> fetchBboxFaceCentroids() {
        return this.boundingBox().fetchCentroids();
    }
}
