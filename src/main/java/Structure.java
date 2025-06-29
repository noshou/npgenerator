import com.oson.tuple.*;
import org.jetbrains.annotations.*;

public abstract class Structure {
    private final Shape shape;
    private final LatticeType lattice_type;

    public Structure(
            Shape shape,
            LatticeType lattice_type
    ) {
        this.shape = shape;
        this.lattice_type = lattice_type;
    }

    // add factory methods here and abstract returns
    public abstract @NotNull UnitCell getBasis();
}
