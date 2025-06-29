import com.oson.tuple.*;
import java.io.*;
import java.nio.file.*;

//
class NpMmcifBuilder {
    private static NpMmcifBuilder instance = null;
    private final String file_name;
    private final BufferedWriter writer;
    private boolean is_finished;
    private NpMmcifBuilder(String file_name) throws IOException {
        // create tmp file
        this.file_name = file_name + ".mmcif";
        Path temp_path = Paths.get(file_name + ".tmp");
        this.writer = Files.newBufferedWriter(temp_path);
    }
        public static NpMmcifBuilder getInstance(String file_name) throws IOException {
        if (instance == null) {
            instance = new NpMmcifBuilder(file_name);
        }
        return instance;
    }
    public void initShape(Shape s) throws IOException {
        writer.write("data_" + s.getStructureIndex() + "\n");
        writer.write("_entry.id " + s.getStructureIndex() + "\n\n");
        writer.write("_cell.entry_idx " + s.getStructureIndex() + "\n");

        // write unit cell side lengths and names
        Polyad<Tuple<String>> cell_lengths = s.getUnitCell().getCellLengths();
        for (int i = 0; i < cell_lengths.fetchSize(); i++) {
            Dyad<String> length = (Dyad<String>)cell_lengths.fetch(i);
            String name = length.fetch(0);
            String size = length.fetch(1);
            writer.write("_cell.length_" + name + " " + size + "\n");
        }

        // write unit cell interaxial angles and names
        Polyad<Tuple<String>> cell_angles = s.getUnitCell().getCellAngles();
        for (int i = 0; i < cell_angles.fetchSize(); i++) {
            Dyad<String> angle = (Dyad<String>)cell_angles.fetch(i);
            String name = angle.fetch(0);
            String size = angle.fetch(1);
            writer.write("_cell.angle_" + name + " " + size + "\n");
        }

        // write symmetry information
        writer.write("_symmetry.entry_id " + s.getStructureIndex() + "\n");
        writer.write("_symmetry.space_group_name_" + s.getUnitCell().getSpaceGroup()+"\n");

        // write loop information; structure now ready for addAtom()
        writer.write("loop_\n");
        writer.write("_atom_site.group_PDB\n");
        writer.write("_atom_site.id\n");
        writer.write("_atom_site.type_symbol\n");
        writer.write("_atom_site.label_atom_id\n");
        writer.write("_atom_site.Cartn_x\n");
        writer.write("_atom_site.Cartn_y\n");
        writer.write("_atom_site.Cartn_z\n");
        writer.write("_atom_site.pdbx_formal_charge\n");
        writer.write("_atom_site.occupancy\n");
        writer.write("_atom_site.auth_asym_id\n");
        writer.write("_atom_site.custom_radius_Ångströms\n");
    }

    public void addAtom(Atom a) throws IOException {
        if (this.is_finished) {
            throw new IllegalStateException("Builder has already been finalized!");
        }
        String line = String.format(
                "HETATM %d %s %s %s %s %s %s 1 A %s\n",
                a.getIndex(),
                a.getElement(),
                String.format("%s%d", a.getElement(), a.getIndex()),
                ((Triad<String>) a.getCentroid()).fetch(0),
                ((Triad<String>) a.getCentroid()).fetch(1),
                ((Triad<String>) a.getCentroid()).fetch(2),
                a.getFormalCharge(),
                a.getRadius()
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
        Files.deleteIfExists(Paths.get(file_name + ".tmp"));
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
                    Paths.get(file_name + ".tmp"),
                    Paths.get(file_name),
                    StandardCopyOption.ATOMIC_MOVE
            );
            is_finished = true;
        }
    }
}