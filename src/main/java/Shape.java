import com.oson.tuple.Polyad;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * Abstract base class representing a geometric nanostructure based on a lattice arrangement of atoms.
 * <p>
 * Encapsulates the unit cell, coordinate system, and physical parameters such as radius and lattice constant
 * required to construct a shape within a specified lattice system.
 * </p>
 *
 * <p><b>Contract:</b> This class is immutable after construction. All fields are non-null unless otherwise noted.
 * Subclasses must implement {@link #build()} to complete the generation of atomic data.</p>
 */
public abstract class Shape {

    /**
     * The unit cell defining the structure's basis and symmetry. Non-null.
     */
    protected @NotNull UnitCell unit_cell;

    /**
     * Coordinate system implementation for generating atomic positions. Non-null.
     */
    protected @NotNull AtomicCoordinates coordinates;

    /**
     * The output file name prefix (without extension) used when writing the final structure. Non-null.
     */
    protected final @NotNull String file_name;

    /**
     * The user-defined name of the structure. Non-null.
     */
    protected final @NotNull String structure_name;

    /**
     * An identifier for the structure instance, typically unique across builds. Non-null.
     */
    protected final @NotNull String structure_index;

    /**
     * The numeric precision (number of significant digits) used for all high-precision computations.
     */
    protected final int precision;

    /**
     * The lattice constant (edge length of the unit cell) as a string to preserve Apfloat precision. Non-null.
     */
    protected final @NotNull Apfloat lattice_constant;

    /**
     * The atomic radius of the constituent atom in angstroms, represented as a string for Apfloat precision. Non-null.
     */
    protected final @NotNull Apfloat radius_angstroms;

    /**
     * Constructs a new shape instance, resolving units and initializing the lattice.
     *
     * @param radius           the radius of the atom as a string (interpreted using {@code radius_type}), non-null
     * @param radius_type      the unit of the radius (must be "pm", "A", "Å", or "nm"), non-null
     * @param lattice_type     the lattice type (currently only {@code LatticeType.FCC} is supported), non-null
     * @param precision        the numeric precision for Apfloat operations
     * @param basis            the atom basis for the unit cell (must contain exactly four atoms for FCC), non-null
     * @param lattice_constant the lattice constant (edge length of unit cell) as a string, non-null
     * @param file_name        the base name for any file output operations, non-null
     * @param structure_name   a user-defined name for the structure, non-null
     * @param structure_index  a unique structure ID used for tracking, non-null
     * @throws IllegalArgumentException if the radius unit or lattice type is not supported
     */
    public Shape(
            @NotNull String radius,
            @NotNull String radius_type,
            @NotNull LatticeType lattice_type,
            int precision,
            @NotNull Polyad<Atom> basis,
            @NotNull String lattice_constant,
            @NotNull String file_name,
            @NotNull String structure_name,
            @NotNull String structure_index
    ) {
        Apfloat TEN = new Apfloat("10", precision);
        Apfloat ONE_HUNDRED = new Apfloat("100", precision);
        Apfloat r = new Apfloat(radius, precision);

        if (
                radius_type.equalsIgnoreCase("pm")
                        || radius_type.equalsIgnoreCase("pico-meters")
                        || radius_type.equalsIgnoreCase("pico meters")
        ) {
            this.radius_angstroms = r.divide(ONE_HUNDRED);
        } else if (
                radius_type.equalsIgnoreCase("A")
                        || radius_type.equalsIgnoreCase("Å")
                        || radius_type.equalsIgnoreCase("Angstrom")
        ) {
            this.radius_angstroms = r;
        } else if (
                radius_type.equalsIgnoreCase("nm")
                        || radius_type.equalsIgnoreCase("nanometer")
        ) {
            this.radius_angstroms = r.multiply(TEN);
        } else {
            throw new IllegalArgumentException("Illegal radius parameter!");
        }

        this.structure_index = structure_index;
        this.file_name = file_name;
        this.structure_name = structure_name;
        this.precision = precision;
        this.lattice_constant = new Apfloat(lattice_constant, this.precision);

        if (lattice_type == LatticeType.FCC) {
            this.coordinates = new FccCoordinates(this.radius_angstroms);
            this.unit_cell = new FccUnitCell(
                    this.lattice_constant,
                    this.precision,
                    basis.fetch(0).getElement(),
                    basis.fetch(0).getFormalChargeInt(),
                    basis.fetch(0).getRadius(),
                    basis.fetch(1).getElement(),
                    basis.fetch(1).getFormalChargeInt(),
                    basis.fetch(1).getRadius(),
                    basis.fetch(2).getElement(),
                    basis.fetch(2).getFormalChargeInt(),
                    basis.fetch(2).getRadius(),
                    basis.fetch(3).getElement(),
                    basis.fetch(3).getFormalChargeInt(),
                    basis.fetch(3).getRadius()
            );
        } else {
            throw new IllegalArgumentException("Illegal lattice type!");
        }
    }

    /**
     * Returns the unit cell used to define atomic basis and lattice geometry.
     *
     * @return the unit cell instance, never null
     */
    @Contract(pure = true)
    public @NotNull UnitCell getUnitCell() {
        return this.unit_cell;
    }

    /**
     * Returns the coordinate system used for placing atoms.
     *
     * @return the coordinate generation implementation, never null
     */
    @Contract(pure = true)
    public @NotNull AtomicCoordinates getCoordinates() {
        return this.coordinates;
    }

    /**
     * Returns a reference to this instance.
     * <p>
     * May be overridden in subclasses that refine the return type.
     * </p>
     *
     * @return this shape instance, never null
     */
    @Contract(pure = true)
    public @NotNull Shape getThis() {
        return this;
    }

    /**
     * Returns the name of the structure.
     *
     * @return the structure name, never null
     */
    @Contract(pure = true)
    public @NotNull String getStructureName() {
        return this.structure_name;
    }

    /**
     * Returns the index of the structure, typically a unique string.
     *
     * @return the structure index, never null
     */
    @Contract(pure = true)
    public @NotNull String getStructureIndex() {
        return this.structure_index;
    }

    /**
     * Returns the lattice constant used for generating the structure.
     *
     * @return the lattice constant as a string, never null
     */
    @Contract(pure = true)
    public @NotNull String getLatticeConstant() {
        return this.lattice_constant.toString();
    }

    /**
     * Returns the file name (without extension) for structure output.
     *
     * @return the file name base, never null
     */
    @Contract(pure = true)
    public @NotNull String getFileName() {
        return this.file_name;
    }

    /**
     * Returns the atomic radius in angstroms as a string.
     *
     * @return the radius in angstroms, never null
     */
    @Contract(pure = true)
    public @NotNull String getRadius() {
        return this.radius_angstroms.toString();
    }

    /**
     * Abstract method that constructs or generates the shape using the lattice and coordinate definitions.
     * <p>
     * Implementing classes must define the full build procedure including atom generation and file output.
     * </p>
     */
    public abstract void build();
}
