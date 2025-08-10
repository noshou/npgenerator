package io.github.noshou.npg.lattice;

import io.github.noshou.tuple.Triad;
import org.apfloat.Apfloat;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * {@code Lattice.FccCoordinates} lazily generates a grid of fractional coordinates for a nanoparticle
 * modeled using a face-centered cubic (FCC) lattice.
 * <p> The coordinate grid spans a cube of lattice points centered at the origin, with 0.5 step size,
 * and bounded by the specified nanoparticle radius.
 * <p> Coordinates are generated on demand using internal counters, which simulate 3 nested loops over x, y, and z.
 * This avoids precomputing and storing large numbers of positions in memory.
 */
public class FccCoordinates implements AtomicCoordinates {

    /** Whether all positions have been iterated. */
    private boolean is_finished = false;

    /** Nanoparticle bounding radius (rounded up), in lattice units. */
    private final @NotNull Apfloat radius;

    /** Constant step of 0.5 units used for FCC grid spacing. */
    private final @NotNull Apfloat face_const;

    /** Constant multiplier -1.0, used to initialize to negative radius. */
    private final @NotNull Apfloat negate;

    /** Current x position in iteration. */
    private @NotNull Apfloat x;

    /** Current y position in iteration. */
    private @NotNull Apfloat y;

    /** Current z position in iteration. */
    private @NotNull Apfloat z;

    /**
     * Increments the current x position by 0.5.
     */
    @Contract(mutates = "this")
    private void incrX() {
        this.x = this.x.add(face_const);
    }

    /**
     * Increments the current y position by 0.5.
     */
    @Contract(mutates = "this")
    private void incrY() {
        this.y = this.y.add(face_const);
    }

    /**
     * Increments the current z position by 0.5.
     */
    @Contract(mutates = "this")
    private void incrZ() {
        this.z = this.z.add(face_const);
    }

    /**
     * Resets x to -radius after completing one full x-sweep.
     */
    @Contract(mutates = "this")
    private void resetX() {
        this.x = this.radius.multiply(negate);
    }

    /**
     * Resets y to -radius after completing one full y-sweep.
     */
    @Contract(mutates = "this")
    private void resetY() {
        this.y = this.radius.multiply(negate);
    }

    /**
     * Constructs a lazily-evaluated FCC grid of positions based on nanoparticle radius.
     * @param radius    the target nanoparticle radius (in same units as lattice constant), must not be null
     */
    @Contract(pure = false)
    public FccCoordinates(
            @NotNull Apfloat radius
    ) {
        this.radius = radius.ceil();
        this.negate = new Apfloat("-1", this.radius.precision());
        this.face_const = new Apfloat("0.5", this.radius.precision());
        this.x = this.radius.multiply(negate);
        this.y = this.radius.multiply(negate);
        this.z = this.radius.multiply(negate);
    }

    /**
     * Returns the next 3D fractional coordinate in the FCC grid.
     * <p>
     * The iteration spans a cube from {@code -radius} to {@code +radius}, using a step of 0.5,
     * and proceeds in nested order: x sweeps fastest, then y, then z.
     * <p>
     * The returned coordinate is a {@link Triad} of string representations suitable for Apfloat-based usage.
     * Once all coordinates have been emitted, returns {@code null}.
     * @return the next FCC grid coordinate as a {@code Triad<String>}, or {@code null} if finished
     */
    @Override
    @Contract(pure = false)
    public @Nullable Triad<@NotNull Apfloat> getPosition() {
        if (this.is_finished) {
            return null;
        }

        // Check if we're at the final point
        if (
                x.compareTo(this.radius) == 0
                && y.compareTo(this.radius) == 0
                && z.compareTo(this.radius) == 0
        ) {
            this.is_finished = true;
        }


        // Generate current position before advancing
        Triad<@NotNull Apfloat> output = new Triad<>(
                this.x,
                this.y,
                this.z
        );

        // Advance to the next position
        incrX();
        if (this.x.compareTo(this.radius) > 0) {
            resetX();
            incrY();
            if (this.y.compareTo(this.radius) > 0) {
                resetY();
                incrZ();
                if (this.z.compareTo(this.radius) > 0) {
                    this.is_finished = true;
                    return null;
                }
            }

        }
        return output;
    }
}
