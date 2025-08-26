package io.github.noshou.npg.npwriter;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.shapes.*;
import io.github.noshou.tuple.*;
import org.jetbrains.annotations.*;
import java.io.IOException;

/**
 * A builder class for constructing mmCIF files representing atomistic nanoparticle shapes.
 * <p>
 * This builder writes to a temporary file (with `.cif` extension) and ensures the output is mmCIF-compliant.
 * Structural metadata such as unit cell dimensions, angles, and symmetry are included, along with atom records.
 * <p>
 * Once finalized, the file cannot be modified.
 */
public class MmCifWriter extends FileWriter {

    /**
     * Constructs a builder that writes to the specified file with a `.cif` extension.
     * @param file_name The base name (without extension) of the mmCIF file.
     * @throws IOException If file creation or access fails.
     */
    public MmCifWriter(@NotNull String file_name) throws IOException {
        super(file_name, ".cif");
    }

    /**
     * Initializes the mmCIF structure with unit cell parameters, symmetry group, and atom loop headers.
     * <p>
     * This method requires the initializer to be a {@link Shape} instance. Any other type will result in an exception.
     * @param initializer The {@link Shape} object from which to extract structural metadata.
     * @throws IOException              If writing to the file fails.
     * @throws IllegalArgumentException If {@code initializer} is null or not of type {@link Shape}.
     */
    @Override
    @Contract("null -> fail")
    public void init(@Nullable Object initializer) throws IOException {
        if (!(initializer instanceof Shape s)) {
            throw new IllegalArgumentException("initializer must be of type Shapes.Shape!");
        }

        writer.write("data_" + s.getStructureIndex() + "\n");
        writer.write("_entry.id " + s.getStructureIndex() + "\n\n");
        writer.write("_cell.entry_idx " + s.getStructureIndex() + "\n");

        // Unit cell lengths
        Polyad<Tuple<String>> cell_lengths = s.getUnitCell().getCellLengths();
        for (int i = 0; i < cell_lengths.fetchSize(); i++) {
            Dyad<String> length = (Dyad<String>) cell_lengths.fetch(i);
            writer.write(
                    "_cell.length_"
                        + length.fetch(0)
                        + " " + length.fetch(1)
                        + "\n"
            );
        }

        // Unit cell angles
        Polyad<Tuple<String>> cell_angles = s.getUnitCell().getCellAngles();
        for (int i = 0; i < cell_angles.fetchSize(); i++) {
            Dyad<String> angle = (Dyad<String>) cell_angles.fetch(i);
            writer.write(
                    "_cell.angle_"
                        + angle.fetch(0)
                        + " "
                        + angle.fetch(1)
                        + "\n"
            );
        }

        // Symmetry
        writer.write(
                    "_symmetry.entry_id "
                        + s.getStructureIndex()
                        + "\n"
        );
        writer.write(
                    "_symmetry.space_group_name_H-M \""
                        + s.getUnitCell().getSpaceGroup()
                        + "\"\n\n"
        );

        // Atom.Atom loop headers
        writer.write("loop_\n");
        writer.write("_atom_site.group_PDB\n");
        writer.write("_atom_site.id\n");
        writer.write("_atom_site.type_symbol\n");
        writer.write("_atom_site.label_atom_id\n");
        writer.write("_atom_site.label_alt_id\n");
        writer.write("_atom_site.label_comp_id\n");
        writer.write("_atom_site.label_asym_id\n");
        writer.write("_atom_site.label_entity_id\n");
        writer.write("_atom_site.label_seq_id\n");
        writer.write("_atom_site.pdbx_PDB_ins_code\n");
        writer.write("_atom_site.Cartn_x\n");
        writer.write("_atom_site.Cartn_y\n");
        writer.write("_atom_site.Cartn_z\n");
        writer.write("_atom_site.occupancy\n");
        writer.write("_atom_site.B_iso_or_equiv\n");
        writer.write("_atom_site.pdbx_formal_charge\n");
        writer.write("_atom_site.auth_seq_id\n");
        writer.write("_atom_site.auth_comp_id\n");
        writer.write("_atom_site.auth_asym_id\n");
        writer.write("_atom_site.auth_atom_id\n");
        writer.write("_atom_site.pdbx_PDB_model_num\n");
    }

    /**
     * Adds an atom to the mmCIF structure.
     * <p>
     * The atom's coordinates, element type, charge, and related metadata are written in a
     * format matching the mmCIF atom site loop. Unused fields are filled with placeholders
     * as per mmCIF convention.
     * @param a The {@link Atom} to add. Must not be null and must have a non-empty element symbol.
     * @throws IOException              If writing to the file fails.
     * @throws IllegalArgumentException If the atom element is null or empty.
     * @throws IllegalStateException    If the builder has already been finalized and closed.
     */
    public void addAtom(@NotNull Atom a) throws IOException {
        if (is_finished) {
            throw new IllegalStateException(
                    "Builder has already been finalized!"
            );
        }
        String element = a.getElement();
        if (element.isEmpty()) {
            throw new IllegalArgumentException(
                    "Atom.Atom element cannot be null or empty"
            );
        }

        // Normalize element symbol
        element = element.substring(0, 1).toUpperCase() +
                (element.length() > 1 ? element.substring(1).toLowerCase() : "");

        Triad<String> coords = (Triad<String>) a.getCentroid();
        String charge = a.getFormalCharge();

        // Format atom entry as per mmCIF loop order
        String[] tokens = new String[] {
                "HETATM",                           // group_PDB
                String.valueOf(a.getIndex()),       // id
                element,                            // type_symbol
                element + a.getIndex(),             // label_atom_id
                ".",                                // label_alt_id
                element,                            // label_comp_id
                "A",                                // label_asym_id
                "1",                                // label_entity_id
                String.valueOf(a.getIndex()),       // label_seq_id
                ".",                                // pdbx_PDB_ins_code
                coords.fetch(0),                // Cartn_x
                coords.fetch(1),                // Cartn_y
                coords.fetch(2),                // Cartn_z
                "1.00",                             // occupancy
                "1.00",                             // B_iso_or_equiv
                charge,                             // formal_charge
                String.valueOf(a.getIndex()),       // auth_seq_id
                element,                            // auth_comp_id
                "A",                                // auth_asym_id
                element + a.getIndex(),             // auth_atom_id
                "1"                                 // model_num
        };

        writer.write(String.join(" ", tokens) + "\n");
    }
}
