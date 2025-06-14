import java.io.*;
import java.nio.file.*;

class NpMmcifBuilder {
    private final String file_name;
    private final BufferedWriter writer;
    private int current_idx = 0;
    private boolean isFinished = false;

    public NpMmcifBuilder(
            String structure_name,
            String data_id,
            String atom_radius_angstroms
    ) throws IOException {

        // Create temporary file
        this.file_name = data_id + ".mmcif";
        Path tempPath = Paths.get(file_name + ".tmp");
        this.writer = Files.newBufferedWriter(tempPath);

        // Write header with proper formatting
        writer.write("data_" + data_id + "\n");
        writer.write("_entry.id " + data_id + "\n\n");
        writer.write("_cell.entry_id " + data_id + "\n");
        writer.write("_cell.length_a " + atom_radius_angstroms + "\n");
        writer.write("_cell.length_b " + atom_radius_angstroms + "\n");
        writer.write("_cell.length_c " + atom_radius_angstroms + "\n");
        writer.write("_cell.angle_alpha 90\n");
        writer.write("_cell.angle_beta 90\n");
        writer.write("_cell.angle_gamma 90\n\n");
        writer.write("_symmetry.entry_id " + data_id + "\n");
        writer.write("_symmetry.space_group_name_H-M 'P 1'\n\n");
        writer.write("loop_\n");
        writer.write("_entity.id\n");
        writer.write("_entity.type\n");
        writer.write("_entity.details\n");
        writer.write(String.format("1 non-polymer man '%s'\n\n", structure_name));
        writer.write("loop_\n");
        writer.write("_atom_site.group_PDB\n");
        writer.write("_atom_site.id\n");
        writer.write("_atom_site.type_symbol\n");
        writer.write("_atom_site.label_atom_id\n");
        writer.write("_atom_site.Cartn_x\n");
        writer.write("_atom_site.Cartn_y\n");
        writer.write("_atom_site.Cartn_z\n");
        writer.write("_atom_site.B_iso_or_equiv\n");
        writer.write("_atom_site.pdbx_formal_charge\n");
        writer.write("_atom_site.occupancy\n");
        writer.write("_atom_site.auth_asym_id\n");
        writer.write("_atom_site.radius_Å\n");
        writer.write("_atom_site.bounding_box\n");
        writer.write("_atom_site.temperature_K");
    }

    public void addAtom(Atom atom) throws IOException {
        if (isFinished) {
            throw new IllegalStateException("Builder has already been finalized");
        }
        if (atom.fetchIndex() <= current_idx) {
            abort();
            throw new IllegalAccessError("Atom index must be unique and 1-indexed!");
        }

        // Format coordinates with proper decimal handling
        String line = new String(String.format(
                "HETATM %d %s %lf %lf %lf %s 1 A %s %s %s %s",
                atom.fetchIndex(),
                String.format(
                        "%sd",
                        atom.fetchElement().toUpperCase(),
                        atom.fetchIndex()
                ),
                atom.fetchCoordinates().fetch(0),
                atom.fetchCoordinates().fetch(1),
                atom.fetchCoordinates().fetch(2),
                atom.bIso(),
                atom.fetchFormalCharge(),
                atom.fetchRadius(),
                atom.boundingBox().toString(),
                atom.temp()
            )
        );
        writer.write(line);
        current_idx++;
    }

    public void writeFile() throws IOException {
        if (!isFinished) {
            writer.close();
            Files.move(
                    Paths.get(file_name + ".tmp"),
                    Paths.get(file_name),
                    StandardCopyOption.ATOMIC_MOVE
            );
            isFinished = true;
        }
    }

    public void abort() throws IOException {
        writer.close();
        Files.deleteIfExists(Paths.get(file_name + ".tmp"));
        isFinished = true;
    }
}


