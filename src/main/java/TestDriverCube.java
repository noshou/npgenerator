import com.oson.tuple.*;

public class TestDriverCube {
    public static void main(String[] args) {

        // create basis set of Au gold atoms
        Atom _1 = new Atom(
                "Au",
                "1.44",
                new Triad<String>("0","0","0"),
                0,
                100
        );
        Atom _2 = new Atom(
                "Au",
                "1.44",
                new Triad<String>("0.5","0.5","0"),
                0,
                100
        );
        Atom _3 = new Atom(
                "Au",
                "1.44",
                new Triad<String>("0.5","0","0.5"),
                0,
                100
        );
        Atom _4 = new Atom(
                "Au",
                "1.44",
                new Triad<>("0","0.5","0.5"),
                0,
                100
        );
        Atom[] atoms = new Atom[]{_1,_2,_3,_4};
        Polyad<Atom> basis = new Polyad<>(atoms);

        // measure time it takes to build
        long start_time = System.nanoTime();
        Shape test = new Cube(
                "1.5",
                "nanometer",
                LatticeType.FCC,
                100,
                basis,
                "4.0786",
                "test",
                "test_sphere",
                "test_me_now"
        );
        test.build(true);
        long end_time =  System.nanoTime();
        long measured = (end_time - start_time) / (long) 6e10;
        System.out.print("Execution time (minutes): " + measured);
    }
}
