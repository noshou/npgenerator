import org.apfloat.*;
import java.io.*;
import java.nio.file.*;

class CoordsDebugWriter {
    private static CoordsDebugWriter instance = null;
    private final String file_name;
    private final BufferedWriter writer;
    private boolean is_finished;
    private CoordsDebugWriter(String file_name) throws IOException {
        // create tmp file
        this.file_name = file_name + ".csv"; // full final name
        Path temp_path = Paths.get(this.file_name + ".tmp"); // → test.cif.tmp
        this.writer = Files.newBufferedWriter(temp_path);
    }
        public static CoordsDebugWriter getInstance(String file_name) throws IOException {
        if (instance == null) {
            instance = new CoordsDebugWriter(file_name);
        }
        return instance;
    }

    public void initLog() throws IOException {
        writer.write(
                "x_frac,y_frac,z_frac,x_cart,y_cart,z_cart,is_occupied\n"
        );
    }


    public void addCoordinate(
            Apfloat x_frac,
            Apfloat y_frac,
            Apfloat z_frac,
            Apfloat x_cart,
            Apfloat y_cart,
            Apfloat z_cart,
            boolean is_occupied
    ) throws IOException {
        String line = String.format(
                "%s,%s,%s,%s,%s,%s,%b\n",
                x_frac.toString(),
                y_frac.toString(),
                z_frac.toString(),
                x_cart.toString(),
                y_cart.toString(),
                z_cart.toString(),
                is_occupied
        );
        writer.write(line);
    }

    /**
      * Aborts the build and deletes the temporary file.
      *
      * @throws IOException If file deletion fails
      */
    public void abort() throws IOException {
        writer.close();
        Files.deleteIfExists(Paths.get(this.file_name + ".tmp"));
    }

    /**
      * Finalizes and writes the file to disk.
      *
      * @throws IOException If file write fails
      */
    public void writeFile() throws IOException {
        if (!is_finished) {
            writer.close();

            Files.move(
                    Paths.get(this.file_name + ".tmp"),   // from
                    Paths.get(this.file_name),            // to
                    StandardCopyOption.ATOMIC_MOVE,
                    StandardCopyOption.REPLACE_EXISTING   // optional: overwrites existing file
            );

            is_finished = true;
        }
    }
}