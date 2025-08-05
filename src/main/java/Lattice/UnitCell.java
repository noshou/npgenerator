package Lattice;

import com.oson.tuple.*;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.*;
import Atom.*;

public abstract class UnitCell {
    /** Digits of precision */
    private final int precision;

    /** Enum representing the lattice type (FCC, BCC, etc.). */
    protected final @NotNull LatticeType lattice_type;

    /** Polyad of atoms forming the asymmetric basis for this unit cell. */
    protected final @NotNull Polyad<Atom> basis;

    /** Hermann–Mauguin space group string (e.g., "P 1", "F m -3 m"). */
    protected final @NotNull String space_group;

    public UnitCell(
            int precision,
            @NotNull LatticeType lattice_type,
            @NotNull Polyad<Atom> basis,
            @NotNull String space_group
    ) {
        this.precision = precision;
        this.lattice_type = lattice_type;
        this.basis = basis;
        this.space_group = space_group;
    }

    /** @return digits of precision */
    @Contract(pure = true)
    protected int getPrecision() {
        return this.precision;
    }
    /**
     * Returns the lattice type (FCC, BCC, etc.) associated with this unit cell.
     *
     * @return non-null enum indicating lattice symmetry
     */
    @Contract(pure = true)
    public @NotNull LatticeType latticeType() {
        return this.lattice_type;
    }

    /** @return Hermann–Mauguin space group string */
    @Contract(pure = true)
    public @NotNull String getSpaceGroup() { return space_group; }

    /** @return Atom.Atom based on index */
    @Contract(pure = true)
    protected @NotNull Atom getAtom(int idx) {
        return this.basis.fetch(idx);
    }

    /** @return unit cell edge lengths (Ångströms) */
    @Contract(pure = true)
    public @NotNull abstract Polyad<Tuple<String>> getCellLengths();

    /** @return interaxial angles  */
    @Contract(pure = true)
    public @NotNull abstract Polyad<Tuple<String>> getCellAngles();

    /**
     * Returns the basis atom to be placed at a specific lattice point.
     * <p>
     * This is implemented by child classes which define the logic for
     * positioning atoms at fractional coordinates based on lattice type
     * and geometry.
     * </p>
     *
     * @param frac_x fractional x-coordinate in unit cell (e.g., "0.5"); must not be null
     * @param frac_y fractional y-coordinate in unit cell (e.g., "0"); must not be null
     * @param frac_z fractional z-coordinate in unit cell (e.g., "0.5"); must not be null
     * @return the corresponding basis {@code Atom.Atom} with position and properties, or {@code null} if no atom exists at this site
     */
    @Contract(pure = true)
    public abstract @Nullable Atom getLatticePoint(Apfloat frac_x, Apfloat frac_y, Apfloat frac_z);

}
