import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;
import java.io.IOException;

/**
 * Abstract base class representing a geometric nanostructure based on a lattice arrangement of atoms.
 * <p>
 * Encapsulates the unit cell, coordinate system, and physical parameters such as radius and lattice constant
 * required to construct a shape within a specified lattice system.
 * </p>
 *
 * <p><b>Contract:</b> This class is immutable after construction. All fields are non-null unless otherwise noted.
 * Subclasses must implement {@link #inBounds()} to complete the generation of atomic data.</p>
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
    @Contract("-> this")
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
    public @NotNull Apfloat getLatticeConstant() {
        return this.lattice_constant;
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
    public @NotNull Apfloat getRadius() {
        return this.radius_angstroms;
    }

    /**
     * Determines whether the Cartesian point (x, y, z) lies within the spatial bounds of the structure.
     * <p>
     * Implementations should define the exact geometric inclusion criteria (e.g., within a sphere or cube).
     * </p>
     *
     * @param x_cart the x-coordinate in Cartesian space, non-null
     * @param y_cart the y-coordinate in Cartesian space, non-null
     * @param z_cart the z-coordinate in Cartesian space, non-null
     * @return true if the point is within the bounds of the shape, false otherwise
     */
    @Contract(pure = true)
    protected abstract boolean inBounds(
            @NotNull Apfloat x_cart,
            @NotNull Apfloat y_cart,
            @NotNull Apfloat z_cart
    );

    /**
     * Builds the atomic structure and writes it to a CIF file.
     * <p>
     * Coordinates are iterated and filtered through {@link #inBounds(Apfloat, Apfloat, Apfloat)}.
     * Each valid coordinate is transformed into a lattice atom and recorded in the output file.
     * </p>
     *
     * <p><b>Contract:</b> This method must be called only once per instance. If writing fails at any point,
     * the temporary output is aborted.</p>
     *
     * @throws RuntimeException if an I/O error occurs during file writing or abortion
     */
    @Contract("-> fail")  // method may throw at runtime
    public void build() {

        // get file instance, initialize shape
        // RADIUS IS IN NANOMETERS !!!
        NpMmcifBuilder file;
        try {
            file = new NpMmcifBuilder(this.file_name);
            file.initShape(this.getThis());
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }

        // write atoms
        int index = 0;
        Triad<Apfloat> curr = this.coordinates.getPosition();

        // loop through coordinates, check if point is in bounds
        while (curr != null) {
            Apfloat x_frac = curr.fetch(0);
            Apfloat y_frac = curr.fetch(1);
            Apfloat z_frac = curr.fetch(2);
            Apfloat x_cart = x_frac.multiply(this.lattice_constant);
            Apfloat y_cart = y_frac.multiply(this.lattice_constant);
            Apfloat z_cart = z_frac.multiply(this.lattice_constant);

            // check if point is in the unit lattice
            if (inBounds(x_cart, y_cart, z_cart)) {

                // if atom is null -> not in unit cell
                // else, write atom to file
                Atom curr_atom = this.getUnitCell().getLatticePoint(
                        x_frac,
                        y_frac,
                        z_frac
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
                        file.addAtom(curr_atom);
                    } catch (IOException e2) {
                        try {
                            file.abort();
                        } catch (IOException abortException) {
                            e2.addSuppressed(abortException);
                        }
                        throw new RuntimeException(e2);
                    }
                }
            }
            index++;
            curr = this.getCoordinates().getPosition();
        }

        // write files
        try {
            file.writeFile();
        } catch (IOException e2) {
            try {
                file.abort();
            } catch (IOException abortException) {
                e2.addSuppressed(abortException);
            }
            throw new RuntimeException(e2);
        }
    }

    /**
     * Builds the atomic structure and writes it to a CIF file, with optional debug coordinate logging.
     * <p>
     * This variant behaves identically to {@link #build()} but also emits a debug trace of included and
     * excluded atoms to a secondary file if {@code debug} is {@code true}.
     * </p>
     *
     * <p><b>Contract:</b> This method must be called only once per instance. If writing or logging fails,
     * all temporary files are aborted.</p>
     *
     * @param debug whether to emit coordinate debug information to an auxiliary file
     * @throws RuntimeException if an I/O error occurs during writing, logging, or abortion
     */
    @Contract("_ -> fail")  // method may throw at runtime
    public void build(boolean debug) {

        // initialize debug log (if indicated)
        CoordsDebugWriter dlog = null;
        if (debug) {
            try {
                dlog = new CoordsDebugWriter("build_debug_"+this.file_name);
                dlog.initLog();
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        // get file instance, initialize shape
        // RADIUS IS IN NANOMETERS !!!
        NpMmcifBuilder file;
        try {
            file = new NpMmcifBuilder(this.file_name);
            file.initShape(this.getThis());
        }
        catch (IOException e) {
            if (dlog != null) {
                try {
                    dlog.abort();
                } catch (IOException ex2) {
                    e.addSuppressed(ex2);  // Optional: add extra context
                }
            }
            throw new RuntimeException(e);
        }

        // write atoms
        int index = 0;
        Triad<Apfloat> curr = this.coordinates.getPosition();

        // loop through coordinates, check if point is in bounds
        while (curr != null) {
            Apfloat x_frac = curr.fetch(0);
            Apfloat y_frac = curr.fetch(1);
            Apfloat z_frac = curr.fetch(2);
            Apfloat x_cart = x_frac.multiply(this.lattice_constant);
            Apfloat y_cart = y_frac.multiply(this.lattice_constant);
            Apfloat z_cart = z_frac.multiply(this.lattice_constant);

            // check if point is in the unit lattice
            if (inBounds(x_cart, y_cart, z_cart)) {

                // if atom is null -> not in unit cell
                // else, write atom to file
                Atom curr_atom = this.getUnitCell().getLatticePoint(
                        x_frac,
                        y_frac,
                        z_frac
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
                    index++;
                    try {
                        file.addAtom(curr_atom);
                        if (dlog != null) {
                            dlog.addCoordinate(
                                    x_frac,
                                    y_frac,
                                    z_frac,
                                    x_cart,
                                    y_cart,
                                    z_cart,
                                    true);
                        }
                    } catch (IOException e2) {
                        try {
                            file.abort();
                            if (dlog != null) {
                                dlog.abort();
                            }
                        } catch (IOException abortException) {
                            e2.addSuppressed(abortException);
                        }
                        throw new RuntimeException(e2);
                    }
                } else if (dlog != null) {
                    try {
                        dlog.addCoordinate(
                                x_frac,
                                y_frac,
                                z_frac,
                                x_cart,
                                y_cart,
                                z_cart,
                                false);
                    } catch (IOException e2) {
                        try {
                            file.abort();
                            dlog.abort();
                        } catch (IOException abortException) {
                            e2.addSuppressed(abortException);
                        }
                        throw new RuntimeException(e2);
                    }
                }
            }
            curr = this.getCoordinates().getPosition();
        }

        // write files
        try {
            file.writeFile();
            if (dlog != null) {
                dlog.writeFile();
            }
        } catch (IOException e2) {
            try {
                file.abort();
                if (dlog != null) {
                    dlog.abort();
                }
            } catch (IOException abortException) {
                e2.addSuppressed(abortException);
            }
            throw new RuntimeException(e2);
        }
    }

}
