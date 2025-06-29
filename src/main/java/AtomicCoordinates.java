import com.oson.tuple.*;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Abstract for lattice-based grids of atomic positions.
 * <p>
 * Provides methods to access all precomputed positions or to consume them one by one.
 * </p>
 */
public abstract class AtomicCoordinates {

    /** A thread-safe queue of 3D fractional positions representing lattice sites. */
    private final ConcurrentLinkedQueue<Triad<String>> grid_positions;

    public AtomicCoordinates(
            ConcurrentLinkedQueue<Triad<String>> grid_positions
    ) {
        this.grid_positions = grid_positions;
    }

    /**
     * Retrieves and removes the next Cartesian position from the grid.
     *
     * @return the next available position as a {@link Triad} of string coordinates, or {@code null} if empty
     */
    public Triad<String> getPosition() {
        return this.grid_positions.poll();
    }

    /**
     * Returns the full queue of all precomputed grid positions.
     *
     * @return a {@link ConcurrentLinkedQueue} of fractional positions as string triads
     */
    public ConcurrentLinkedQueue<Triad<String>> getPositions() {
        return this.grid_positions;
    }

}