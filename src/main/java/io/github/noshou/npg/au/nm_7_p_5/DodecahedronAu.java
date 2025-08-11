package io.github.noshou.npg.au.nm_7_p_5;

import io.github.noshou.npg.atom.Atom;
import io.github.noshou.npg.lattice.LatticeType;
import io.github.noshou.npg.nputil.FormatExecTime;
import io.github.noshou.npg.shapes.*;
import io.github.noshou.npg.shapes.platonic.Dodecahedron;
import io.github.noshou.tuple.*;
import org.jetbrains.annotations.*;
/**
 * Constructs and benchmarks a gold (Au) nanoparticle using
 * a face-centered cubic (FCC) lattice.
 * <p>
 * This class generates a default FCC atomic basis for gold and builds a
 * spherical nanoparticle using the {@link Dodecahedron} class,
 * configured with user-defined radius, precision, and lattice parameters.
 */
public class DodecahedronAu {

    /**
     * Returns a default FCC basis for gold (Au) containing four atoms
     * located at the standard fractional positions within the FCC unit cell.
     * <p>
     * Each atom is initialized with:
     * <ul>
     *   <li>Element: "Au"</li>
     *   <li>Atomic radius: 1.44 Å (as a string)</li>
     *   <li>Fractional coordinates: FCC standard positions</li>
     *   <li>B-factor: 0</li>
     *   <li>Unique index: 100</li>
     * </ul>
     * @return a non-null {@link Polyad} containing the FCC atomic basis for gold
     * @implNote All atoms share the same radius and index,
     * and are positioned for symmetry.
     */
    @Contract(" -> new")
    private static @NotNull Polyad<Atom> getBasis() {
        Atom _1 = new Atom(
                "Au",
                "1.44",
                new Triad<>("0", "0", "0"),
                0,
                100
        );
        Atom _2 = new Atom(
                "Au",
                "1.44",
                new Triad<>("0.5", "0.5", "0"),
                0,
                100
        );
        Atom _3 = new Atom(
                "Au",
                "1.44",
                new Triad<>("0.5", "0", "0.5"),
                0,
                100
        );
        Atom _4 = new Atom(
                "Au",
                "1.44",
                new Triad<>("0", "0.5", "0.5"),
                0,
                100
        );
        Atom[] atoms = new Atom[]{_1, _2, _3, _4};
        return new Polyad<>(atoms);
    }

    /**
     * The main entry point for building and benchmarking a
     * gold FCC spherical nanoparticle.
     * <p>
     * It creates a {@link Dodecahedron} using a specified radius and lattice constant,
     * builds the structure, and prints the execution time in a human-readable format.
     * @param args the command-line arguments (not used)
     */
    public static void main(String @NotNull [] args) {
        long start_time_init = System.nanoTime();

        // Configuration parameters
        String radius = "7.5"; // in nanometers
        String units = "nm";
        LatticeType lattice_type = LatticeType.FCC;
        int precision = 500;
        Polyad<Atom> atoms = getBasis();
        String lattice_constant = "4.33"; // Ångstroms
        String id = "DodecahedronAu";

        System.out.println("\nNanoparticle Radius:  \t" + radius + " " + units);

        long end_time_init = System.nanoTime();
        System.out.println(
                "Init time:   \t" + FormatExecTime.formatDuration(
                        (end_time_init - start_time_init)
                ) + "\n"
        );
        // Create and build the nanoparticle
        Shape auDodecahedron = new Dodecahedron(
                radius,
                units,
                lattice_type,
                precision,
                atoms,
                lattice_constant,
                id, // label
                id, // CIF ID
                id  // output file stem
        );

        auDodecahedron.build();

        long end_time_build = System.nanoTime();

        System.out.println(
                "Build time:  \t" + FormatExecTime.formatDuration(
                        (end_time_build - end_time_init)
                ) + "\n"
        );
        System.out.println(
                "Exec time:   \t" + FormatExecTime.formatDuration(
                        (end_time_build - start_time_init)
                ) + "\n"
        );
    }
}
