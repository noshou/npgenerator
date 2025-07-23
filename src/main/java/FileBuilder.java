import org.jetbrains.annotations.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.*;

/**
 * An abstract base class for safely building files with atomic finalization.
 * <p>
 * The {@code FileBuilder} writes all data to a temporary file (with `.tmp` extension) first.
 * Once finalized, the file is atomically renamed to its final name to prevent partial writes
 * in case of errors or interruption. Subclasses must implement the {@link #init(Object)} method.
 */
abstract class FileBuilder {

    /** The full file name including the extension (excluding `.tmp` suffix during writing). */
    private final String file_name;

    /** Writer for the temporary file. Accessible to subclasses for writing content. */
    protected final BufferedWriter writer;

    /** Flag indicating whether {@link #writeFile()} has been called. */
    protected boolean is_finished = false;

    /**
     * Constructs a file builder that writes to a temporary file.
     * The final file will be named {@code file_name + extension}.
     *
     * @param file_name The base name of the file (without extension).
     * @param extension The extension to append (e.g., ".cif", ".csv").
     * @throws IOException If the temporary file cannot be created or opened.
     */
    public FileBuilder(
            @NotNull String file_name,
            @NotNull String extension
    ) throws IOException {
        this.file_name = file_name + extension;
        Path temp_path = Paths.get(this.file_name + ".tmp");
        this.writer = Files.newBufferedWriter(temp_path);
    }

    /**
     * Aborts the build and deletes the temporary file.
     * <p>
     * After calling this method, the builder is considered invalid.
     *
     * @throws IOException If the temporary file cannot be deleted.
     */
    public void abort() throws IOException {
        writer.close();
        Files.deleteIfExists(Paths.get(this.file_name + ".tmp"));
    }

    /**
     * Finalizes and renames the file to its intended name, replacing any existing file.
     * <p>
     * This method performs an atomic move to ensure that the file either fully exists or
     * not at all. It is safe to call this method only once.
     *
     * @throws IOException If the temporary file cannot be closed or moved.
     */
    public void writeFile() throws IOException {
        if (!is_finished) {
            writer.close();
            Files.move(
                    Paths.get(this.file_name + ".tmp"),
                    Paths.get(this.file_name),
                    StandardCopyOption.ATOMIC_MOVE,
                    StandardCopyOption.REPLACE_EXISTING
            );
            is_finished = true;
        }
    }

    /**
     * Initializes the builder's output by writing metadata or headers.
     * <p>
     * Subclasses should define the semantics and expected type of {@code initializer}.
     *
     * @param initializer               A context-specific object (e.g., a Shape or null) used to generate headers.
     * @throws IOException              If writing fails.
     * @throws IllegalArgumentException If the initializer is invalid for the implementation.
     */
    public abstract void init(@Nullable Object initializer) throws IOException;
}
