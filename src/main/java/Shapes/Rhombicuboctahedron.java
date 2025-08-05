package Shapes;

import Utilities.VectorMath;
import com.oson.tuple.*;
import org.apfloat.*;
import org.jetbrains.annotations.*;
import Lattice.*;
import Atom.*;

public class Rhombicuboctahedron extends Shape {

    // the 8 triangular faces
    private final Octad<Tuple<Tuple<Apfloat>>> faces_tri;
    private final Octad<Tuple<Apfloat>> face_norms_tri;

    // the 18 square faces
    private final Octakaidecad<Tuple<Tuple<Apfloat>>> faces_sqr;
    private final Octakaidecad<Tuple<Apfloat>> face_norms_sqr;

    // Apfloat constants
    private final Apfloat NEG = new Apfloat("-1", super.precision);
    private final Apfloat ZERO = new Apfloat("0", super.precision);
    private final Apfloat ONE = new Apfloat("1", super.precision);
    private final Apfloat TWO = new Apfloat("2", super.precision);
    private final Apfloat THREE = new Apfloat("3", super.precision);
    private final Apfloat FOUR = new Apfloat("4", super.precision);
    private final Apfloat ONE_PLUS_SQRT2 = ONE.add(ApfloatMath.sqrt(TWO));

    public Rhombicuboctahedron(
            @NotNull String radius,
            @NotNull String radius_type,
            @NotNull LatticeType lattice_type,
            int precision,
            @NotNull Polyad<Atom> basis,
            @NotNull String lattice_constant,
            @NotNull String file_name,
            @NotNull String structure_name,
            @NotNull String structure_index
    ) {
        super(
                radius,
                radius_type,
                lattice_type,
                precision,
                basis,
                lattice_constant,
                file_name,
                structure_name,
                structure_index
        );

        // ==== CANONICAL BASIS VERTICES ====

        // (±1,±1,±(1+√2))
        Triad<Apfloat> vB0 = new Triad<>(ONE, ONE, ONE_PLUS_SQRT2);                     // (+1,+1,+(1+√2))
        Triad<Apfloat> vB1 = new Triad<>(ONE, ONE, ONE_PLUS_SQRT2.multiply(NEG));       // (+1,+1,-(1+√2))
        Triad<Apfloat> vB2 = new Triad<>(ONE, NEG, ONE_PLUS_SQRT2);                     // (+1,-1,+(1+√2))
        Triad<Apfloat> vB3 = new Triad<>(NEG, ONE, ONE_PLUS_SQRT2);                     // (-1,+1,+(1+√2))
        Triad<Apfloat> vB4 = new Triad<>(NEG, NEG, ONE_PLUS_SQRT2);                     // (-1,-1,+(1+√2))
        Triad<Apfloat> vB5 = new Triad<>(NEG, NEG, ONE_PLUS_SQRT2.multiply(NEG));       // (-1,-1,-(1+√2))
        Triad<Apfloat> vB6 = new Triad<>(NEG, ONE, ONE_PLUS_SQRT2.multiply(NEG));       // (-1,+1,-(1+√2))
        Triad<Apfloat> vB7 = new Triad<>(ONE, NEG, ONE_PLUS_SQRT2.multiply(NEG));       // (+1,-1,-(1+√2))

        // (±1,±(1+√2),±1)
        Triad<Apfloat> vB8 = new Triad<>(ONE, ONE_PLUS_SQRT2, ONE);                     // (+1,+(1+√2),+1)
        Triad<Apfloat> vB9 = new Triad<>(ONE, ONE_PLUS_SQRT2, NEG);                     // (+1,+(1+√2),-1)
        Triad<Apfloat> vB10= new Triad<>(ONE, ONE_PLUS_SQRT2.multiply(NEG), ONE);       // (+1,-(1+√2),+1)
        Triad<Apfloat> vB11= new Triad<>(NEG, ONE_PLUS_SQRT2, ONE);                     // (-1,+(1+√2),+1)
        Triad<Apfloat> vB12= new Triad<>(NEG, ONE_PLUS_SQRT2.multiply(NEG), ONE);       // (-1,-(1+√2),+1)
        Triad<Apfloat> vB13= new Triad<>(NEG, ONE_PLUS_SQRT2.multiply(NEG), NEG);       // (-1,-(1+√2),-1)
        Triad<Apfloat> vB14= new Triad<>(NEG, ONE_PLUS_SQRT2, NEG);                     // (-1,+(1+√2),-1)
        Triad<Apfloat> vB15= new Triad<>(ONE, ONE_PLUS_SQRT2.multiply(NEG), NEG);       // (+1,-(1+√2),-1)

        // (±(1+√2),±1,±1)
        Triad<Apfloat> vB16= new Triad<>(ONE_PLUS_SQRT2, ONE, ONE);                     // (+(1+√2),+1,+1)
        Triad<Apfloat> vB17= new Triad<>(ONE_PLUS_SQRT2, ONE, NEG);                     // (+(1+√2),+1,-1)
        Triad<Apfloat> vB18= new Triad<>(ONE_PLUS_SQRT2, NEG, ONE);                     // (+(1+√2),-1,+1)
        Triad<Apfloat> vB19= new Triad<>(ONE_PLUS_SQRT2.multiply(NEG), ONE, ONE);       // (-(1+√2),+1,+1)
        Triad<Apfloat> vB20= new Triad<>(ONE_PLUS_SQRT2.multiply(NEG), NEG, ONE);       // (-(1+√2),-1,+1)
        Triad<Apfloat> vB21= new Triad<>(ONE_PLUS_SQRT2.multiply(NEG), NEG, NEG);       // (-(1+√2),-1,-1)
        Triad<Apfloat> vB22= new Triad<>(ONE_PLUS_SQRT2.multiply(NEG), ONE, NEG);       // (-(1+√2),+1,-1)
        Triad<Apfloat> vB23= new Triad<>(ONE_PLUS_SQRT2, NEG, NEG);                     // (+(1+√2),-1,-1)

        // ==== SCALING FACTOR ====
        // dist from origin = sqrt(x²+y²+z²) = sqrt(1²+1²+(1+√2)²) = sqrt(2+(1+√2)²)
        Apfloat dist = ApfloatMath.sqrt(ONE.add(ONE).add(ApfloatMath.pow(ONE_PLUS_SQRT2, 2)));
        Apfloat r = super.getRadius();
        String r_dist = r.divide(dist).toString();

        // ==== SCALED CARTESIAN VERTICES ====

        // (±(r/sqrt(2+(1+√2)²),±(r/sqrt(2+(1+√2)²),±((r+r*√2)/(sqrt(2+(1+√2)²))
        Triad<Apfloat> vC0 = VectorMath.mult(vB0, r_dist);  // (+(r/sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²),+((r+r*√2)/(sqrt(2+(1+√2)²))
        Triad<Apfloat> vC1 = VectorMath.mult(vB1, r_dist);  // (+(r/sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²),-((r+r*√2)/(sqrt(2+(1+√2)²))
        Triad<Apfloat> vC2 = VectorMath.mult(vB2, r_dist);  // (+(r/sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²),+((r+r*√2)/(sqrt(2+(1+√2)²))
        Triad<Apfloat> vC3 = VectorMath.mult(vB3, r_dist);  // (-(r/sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²),+((r+r*√2)/(sqrt(2+(1+√2)²))
        Triad<Apfloat> vC4 = VectorMath.mult(vB4, r_dist);  // (-(r/sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²),+((r+r*√2)/(sqrt(2+(1+√2)²))
        Triad<Apfloat> vC5 = VectorMath.mult(vB5, r_dist);  // (-(r/sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²),-((r+r*√2)/(sqrt(2+(1+√2)²))
        Triad<Apfloat> vC6 = VectorMath.mult(vB6, r_dist);  // (-(r/sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²),-((r+r*√2)/(sqrt(2+(1+√2)²))
        Triad<Apfloat> vC7 = VectorMath.mult(vB7, r_dist);  // (+(r/sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²),-((r+r*√2)/(sqrt(2+(1+√2)²))

        // (±(r/sqrt(2+(1+√2)²),±((r+r*√2)/(sqrt(2+(1+√2)²),±(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC8 = VectorMath.mult(vB8, r_dist);  // (+(r/sqrt(2+(1+√2)²),+((r+r*√2)/(sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC9 = VectorMath.mult(vB9, r_dist);  // (+(r/sqrt(2+(1+√2)²),+((r+r*√2)/(sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC10= VectorMath.mult(vB10,r_dist);  // (+(r/sqrt(2+(1+√2)²),-((r+r*√2)/(sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC11= VectorMath.mult(vB11,r_dist);  // (-(r/sqrt(2+(1+√2)²),+((r+r*√2)/(sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC12= VectorMath.mult(vB12,r_dist);  // (-(r/sqrt(2+(1+√2)²),-((r+r*√2)/(sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC13= VectorMath.mult(vB13,r_dist);  // (-(r/sqrt(2+(1+√2)²),-((r+r*√2)/(sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC14= VectorMath.mult(vB14,r_dist);  // (-(r/sqrt(2+(1+√2)²),+((r+r*√2)/(sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC15= VectorMath.mult(vB15,r_dist);  // (+(r/sqrt(2+(1+√2)²),-((r+r*√2)/(sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²))

        // (±((r+r*√2)/(sqrt(2+(1+√2)²),±(r/sqrt(2+(1+√2)²),±(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC16= VectorMath.mult(vB16,r_dist);  // (+((r+r*√2)/(sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC17= VectorMath.mult(vB17,r_dist);  // (+((r+r*√2)/(sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC18= VectorMath.mult(vB18,r_dist);  // (+((r+r*√2)/(sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC19= VectorMath.mult(vB19,r_dist);  // (-((r+r*√2)/(sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC20= VectorMath.mult(vB20,r_dist);  // (-((r+r*√2)/(sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC21= VectorMath.mult(vB21,r_dist);  // (-((r+r*√2)/(sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC22= VectorMath.mult(vB22,r_dist);  // (-((r+r*√2)/(sqrt(2+(1+√2)²),+(r/sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²))
        Triad<Apfloat> vC23= VectorMath.mult(vB23,r_dist);  // (+((r+r*√2)/(sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²),-(r/sqrt(2+(1+√2)²))


        // ==== AXIS ALIGNED QUADRILATERAL FACE DEFINITIONS ====

        // z = +(1+√2)
        // Vertices: vB0(+1,+1,+(1+√2)), vB2(+1,-1,+(1+√2)), vB4(-1,-1,+(1+√2)), vB3(-1,+1,+(1+√2))
        Tetrad<Tuple<Apfloat>> posZ = new Tetrad<>(vC0,vC2,vC4,vC3);    // (±x,±y,+(1+√2))
        Triad<Apfloat> posZ_norm = VectorMath.normalQuad(vC0,vC2,vC4,vC3,true);

        // z = -(1+√2)
        // Vertices: vB1(+1,+1,-(1+√2)), vB7(+1,-1,-(1+√2)), vB5(-1,-1,-(1+√2)), vB6(-1,+1,-(1+√2))
        Tetrad<Tuple<Apfloat>> negZ = new Tetrad<>(vC1,vC7,vC5,vC6);    // (±x,±y,-(1+√2))
        Triad<Apfloat> negZ_norm = VectorMath.normalQuad(vC1,vC7,vC5,vC6,true);

        // y = +(1+√2)
        // Vertices: vB8(+1,+(1+√2),+1), vB11(-1,+(1+√2),+1), vB14(-1,+(1+√2),-1), vB9(+1,+(1+√2),-1)
        Tetrad<Tuple<Apfloat>> posY = new Tetrad<>(vC8,vC11,vC14,vC9);  // (±x,+(1+√2),±z)
        Triad<Apfloat> posY_norm = VectorMath.normalQuad(vC8,vC11,vC14,vC9,true);

        // y = -(1+√2)
        // Vertices: vB10(+1,-(1+√2),+1), vB15(+1,-(1+√2),-1), vB13(-1,-(1+√2),-1), vB12(-1,-(1+√2),+1)
        Tetrad<Tuple<Apfloat>> negY = new Tetrad<>(vC10,vC15,vC13,vC12); // (±x,-(1+√2),±z)
        Triad<Apfloat> negY_norm = VectorMath.normalQuad(vC10,vC15,vC13,vC12,true);

        // x = +(1+√2)
        // Vertices: vB16(+(1+√2),+1,+1), vB17(+(1+√2),+1,-1), vB23(+(1+√2),-1,-1), vB18(+(1+√2),-1,+1)
        Tetrad<Tuple<Apfloat>> posX = new Tetrad<>(vC16,vC17,vC23,vC18); // (+(1+√2),±y,±z)
        Triad<Apfloat> posX_norm = VectorMath.normalQuad(vC16,vC17,vC23,vC18,true);

        // x = -(1+√2)
        // Vertices: vB19(-(1+√2),+1,+1), vB22(-(1+√2),+1,-1), vB21(-(1+√2),-1,-1), vB20(-(1+√2),-1,+1)
        Tetrad<Tuple<Apfloat>> negX = new Tetrad<>(vC19,vC22,vC21,vC20); // (-(1+√2),±y,±z)
        Triad<Apfloat> negX_norm = VectorMath.normalQuad(vC19,vC22,vC21,vC20,true);

        // ==== DIAGONAL QUADRILATERAL FACE DEFINITIONS ====
        // Faces bridging between coordinate planes at 45° angles

        // posX_posY: vertices in positive x and y quadrant
        // vB0(+1,+1,+(1+√2)), vB8(+1,+(1+√2),+1), vB16(+(1+√2),+1,+1), vB3(-1,+1,+(1+√2))
        Tetrad<Tuple<Apfloat>> posXposY = new Tetrad<>(vC0,vC8,vC16,vC3);
        Triad<Apfloat> posXposY_norm = VectorMath.normalQuad(vC0,vC8,vC16,vC3,true);

        // posXnegY: vertices in positive x, negative y quadrant
        // vB1(+1,+1,-(1+√2)), vB9(+1,+(1+√2),-1), vB17(+(1+√2),+1,-1), vB6(-1,+1,-(1+√2))
        Tetrad<Tuple<Apfloat>> posXnegY = new Tetrad<>(vC1,vC9,vC17,vC6);
        Triad<Apfloat> posXnegY_norm = VectorMath.normalQuad(vC1,vC9,vC17,vC6,true);

        // negXposY: vertices in negative x, positive y quadrant
        // vB2(+1,-1,+(1+√2)), vB10(+1,-(1+√2),+1), vB18(+(1+√2),-1,+1), vB4(-1,-1,+(1+√2))
        Tetrad<Tuple<Apfloat>> negXposY = new Tetrad<>(vC2,vC10,vC18,vC4);
        Triad<Apfloat> negXposY_norm = VectorMath.normalQuad(vC2,vC10,vC18,vC4,true);

        // negXnegY: vertices in negative x, negative y quadrant
        // vB7(+1,-1,-(1+√2)), vB15(+1,-(1+√2),-1), vB23(+(1+√2),-1,-1), vB5(-1,-1,-(1+√2))
        Tetrad<Tuple<Apfloat>> negXnegY = new Tetrad<>(vC7,vC15,vC23,vC5);
        Triad<Apfloat> negXnegY_norm = VectorMath.normalQuad(vC7,vC15,vC23,vC5,true);

        // posYposZ: vertices in positive y, positive z quadrant
        // vB0(+1,+1,+(1+√2)), vB16(+(1+√2),+1,+1), vB17(+(1+√2),+1,-1), vB1(+1,+1,-(1+√2))
        Tetrad<Tuple<Apfloat>> posYposZ = new Tetrad<>(vC0,vC16,vC17,vC1);
        Triad<Apfloat> posYposZ_norm = VectorMath.normalQuad(vC0,vC16,vC17,vC1,true);

        // posYnegZ: vertices in positive y, negative z quadrant
        // vB2(+1,-1,+(1+√2)), vB18(+(1+√2),-1,+1), vB23(+(1+√2),-1,-1), vB7(+1,-1,-(1+√2))
        Tetrad<Tuple<Apfloat>> posYnegZ = new Tetrad<>(vC2,vC18,vC23,vC7);
        Triad<Apfloat> posYnegZ_norm = VectorMath.normalQuad(vC2,vC18,vC23,vC7,true);

        // negYposZ: vertices in negative y, positive z quadrant
        // vB3(-1,+1,+(1+√2)), vB19(-(1+√2),+1,+1), vB22(-(1+√2),+1,-1), vB6(-1,+1,-(1+√2))
        Tetrad<Tuple<Apfloat>> negYposZ = new Tetrad<>(vC3,vC19,vC22,vC6);
        Triad<Apfloat> negYposZ_norm = VectorMath.normalQuad(vC3,vC19,vC22,vC6,true);

        // negYnegZ: vertices in negative y, negative z quadrant
        // vB4(-1,-1,+(1+√2)), vB20(-(1+√2),-1,+1), vB21(-(1+√2),-1,-1), vB5(-1,-1,-(1+√2))
        Tetrad<Tuple<Apfloat>> negYnegZ = new Tetrad<>(vC4,vC20,vC21,vC5);
        Triad<Apfloat> negYnegZ_norm = VectorMath.normalQuad(vC4,vC20,vC21,vC5,true);

        // posZposX: vertices in positive z, positive x quadrant
        // vB8(+1,+(1+√2),+1), vB16(+(1+√2),+1,+1), vB18(+(1+√2),-1,+1), vB10(+1,-(1+√2),+1)
        Tetrad<Tuple<Apfloat>> posZposX = new Tetrad<>(vC8,vC16,vC18,vC10);
        Triad<Apfloat> posZposX_norm = VectorMath.normalQuad(vC8,vC16,vC18,vC10,true);

        // posZnegX: vertices in positive z, negative x quadrant
        // vB9(+1,+(1+√2),-1), vB17(+(1+√2),+1,-1), vB23(+(1+√2),-1,-1), vB15(+1,-(1+√2),-1)
        Tetrad<Tuple<Apfloat>> posZnegX = new Tetrad<>(vC9,vC17,vC23,vC15);
        Triad<Apfloat> posZnegX_norm = VectorMath.normalQuad(vC9,vC17,vC23,vC15,true);

        // negZposX: vertices in negative z, positive x quadrant
        // vB11(-1,+(1+√2),+1), vB19(-(1+√2),+1,+1), vB20(-(1+√2),-1,+1), vB12(-1,-(1+√2),+1)
        Tetrad<Tuple<Apfloat>> negZposX = new Tetrad<>(vC11,vC19,vC20,vC12);
        Triad<Apfloat> negZposX_norm = VectorMath.normalQuad(vC11,vC19,vC20,vC12,true);

        // negZnegX: vertices in negative z, negative x quadrant
        // vB14(-1,+(1+√2),-1), vB22(-(1+√2),+1,-1), vB21(-(1+√2),-1,-1), vB13(-1,-(1+√2),-1)
        Tetrad<Tuple<Apfloat>> negZnegX = new Tetrad<>(vC14,vC22,vC21,vC13);
        Triad<Apfloat> negZnegX_norm = VectorMath.normalQuad(vC14,vC22,vC21,vC13,true);

        faces_sqr = new Octakaidecad<>(
                posZ,
                negZ,
                posY,
                negY,
                posX,
                negX,
                posXposY,
                posXnegY,
                negXposY,
                negXnegY,
                posYposZ,
                posYnegZ,
                negYposZ,
                negYnegZ,
                posZposX,
                posZnegX,
                negZposX,
                negZnegX
        );

        face_norms_sqr = new Octakaidecad<>(
                posZ_norm,
                negZ_norm,
                posY_norm,
                negY_norm,
                posX_norm,
                negX_norm,
                posXposY_norm,
                posXnegY_norm,
                negXposY_norm,
                negXnegY_norm,
                posYposZ_norm,
                posYnegZ_norm,
                negYposZ_norm,
                negYnegZ_norm,
                posZposX_norm,
                posZnegX_norm,
                negZposX_norm,
                negZnegX_norm

        );

        // ==== TRIANGULAR FACE DEFINITIONS ====
        // Each triangular face connects one vertex from each coordinate type
        // Named by the coordinate signs of the three vertices in CCW order

        // posXposYposZ: vertices with positive coordinates
        // vB0(+1,+1,+(1+√2)), vB8(+1,+(1+√2),+1), vB16(+(1+√2),+1,+1)
        Triad<Tuple<Apfloat>> posXposYposZ = new Triad<>(vC0,vC8,vC16);
        Triad<Apfloat> posXposYposZ_norm = VectorMath.normalTriple(vC0,vC8,vC16,true);

        // posXposYnegZ: vertices with positive x,y and negative z
        // vB1(+1,+1,-(1+√2)), vB9(+1,+(1+√2),-1), vB17(+(1+√2),+1,-1)
        Triad<Tuple<Apfloat>> posXposYnegZ = new Triad<>(vC1,vC9,vC17);
        Triad<Apfloat> posXposYnegZ_norm = VectorMath.normalTriple(vC1,vC9,vC17,true);

        // posXnegYposZ: vertices with positive x,z and negative y
        // vB2(+1,-1,+(1+√2)), vB10(+1,-(1+√2),+1), vB18(+(1+√2),-1,+1)
        Triad<Tuple<Apfloat>> posXnegYposZ = new Triad<>(vC2,vC10,vC18);
        Triad<Apfloat> posXnegYposZ_norm = VectorMath.normalTriple(vC2,vC10,vC18,true);

        // posXnegYnegZ: vertices with positive x and negative y,z
        // vB7(+1,-1,-(1+√2)), vB15(+1,-(1+√2),-1), vB23(+(1+√2),-1,-1)
        Triad<Tuple<Apfloat>> posXnegYnegZ = new Triad<>(vC7,vC15,vC23);
        Triad<Apfloat> posXnegYnegZ_norm = VectorMath.normalTriple(vC7,vC15,vC23,true);

        // negXposYposZ: vertices with negative x and positive y,z
        // vB3(-1,+1,+(1+√2)), vB11(-1,+(1+√2),+1), vB19(-(1+√2),+1,+1)
        Triad<Tuple<Apfloat>> negXposYposZ = new Triad<>(vC3,vC11,vC19);
        Triad<Apfloat> negXposYposZ_norm = VectorMath.normalTriple(vC3,vC11,vC19,true);

        // negXposYnegZ: vertices with negative x,z and positive y
        // vB6(-1,+1,-(1+√2)), vB14(-1,+(1+√2),-1), vB22(-(1+√2),+1,-1)
        Triad<Tuple<Apfloat>> negXposYnegZ = new Triad<>(vC6,vC14,vC22);
        Triad<Apfloat> negXposYnegZ_norm = VectorMath.normalTriple(vC6,vC14,vC22,true);

        // negXnegYposZ: vertices with negative x,y and positive z
        // vB4(-1,-1,+(1+√2)), vB12(-1,-(1+√2),+1), vB20(-(1+√2),-1,+1)
        Triad<Tuple<Apfloat>> negXnegYposZ = new Triad<>(vC4,vC12,vC20);
        Triad<Apfloat> negXnegYposZ_norm = VectorMath.normalTriple(vC4,vC12,vC20,true);

        // negXnegYnegZ: vertices with all negative coordinates
        // vB5(-1,-1,-(1+√2)), vB13(-1,-(1+√2),-1), vB21(-(1+√2),-1,-1)
        Triad<Tuple<Apfloat>> negXnegYnegZ = new Triad<>(vC5,vC13,vC21);
        Triad<Apfloat> negXnegYnegZ_norm = VectorMath.normalTriple(vC5,vC13,vC21,true);

        faces_tri = new Octad<>(
                posXposYposZ,
                posXposYnegZ,
                posXnegYposZ,
                posXnegYnegZ,
                negXposYposZ,
                negXposYnegZ,
                negXnegYposZ,
                negXnegYnegZ
        );

        face_norms_tri = new Octad<>(
                posXposYposZ_norm,
                posXposYnegZ_norm,
                posXnegYposZ_norm,
                posXnegYnegZ_norm,
                negXposYposZ_norm,
                negXposYnegZ_norm,
                negXnegYposZ_norm,
                negXnegYnegZ_norm
        );
    }

    /**
     * Tests whether the given Cartesian point lies inside or on the boundary
     * of the Rhombicuboctahedron.
     *
     * <p>The test is performed by checking dot products against all
     * face normals (squares and triangles).</p>
     *
     * @param point_cart the Cartesian coordinate point
     * @return {@code true} if inside or on surface, {@code false} if outside
     */
    @Override
    @Contract(pure = true)
    protected boolean inBounds(@NotNull Triad<Apfloat> point_cart) {

        for (int i = 0; i < faces_sqr.fetchSize(); i++) {

            if (i < faces_tri.fetchSize()) {
                // === Check triangular faces ===
                Triad<Tuple<Apfloat>> face = (Triad<Tuple<Apfloat>>) faces_tri.fetch(i);
                Triad<Apfloat> vert0 = (Triad<Apfloat>) face.fetch(0);
                Triad<Apfloat> vert1 = (Triad<Apfloat>) face.fetch(1);
                Triad<Apfloat> vert2 = (Triad<Apfloat>) face.fetch(2);

                // centroid of triangle
                Triad<Apfloat> centroid = new Triad<>(
                        vert0.fetch(0).add(vert1.fetch(0)).add(vert2.fetch(0)).divide(THREE),
                        vert0.fetch(1).add(vert1.fetch(1)).add(vert2.fetch(1)).divide(THREE),
                        vert0.fetch(2).add(vert1.fetch(2)).add(vert2.fetch(2)).divide(THREE)
                );

                // given point p and vertA, calculate vector from centroid -> p:
                // m = p - centroid = (p_x - centroid_x, p_x - centroid_y, p_x -centroid_z)
                Triad<Apfloat> m = VectorMath.subs(point_cart, centroid);

                // compute dot product of m and face_norm:
                // d = face_norm ⋅ m
                //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
                //   = face_norm_x*(p_x-vertA_x)
                //   + face_norm_y*(p_y-vertA_y)
                //   + face_norm_z*(p_z-vertA_z)
                Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_tri.fetch(i);
                Apfloat d = VectorMath.dot_prod(norm, m);

                // point outside if dot product > 0
                if (d.compareTo(ZERO) > 0) {
                    return false;
                }
            }

            // === Check square faces ===
            Tetrad<Tuple<Apfloat>> face = (Tetrad<Tuple<Apfloat>>) faces_sqr.fetch(i);
            Triad<Apfloat> vert0 = (Triad<Apfloat>) face.fetch(0);
            Triad<Apfloat> vert1 = (Triad<Apfloat>) face.fetch(1);
            Triad<Apfloat> vert2 = (Triad<Apfloat>) face.fetch(2);
            Triad<Apfloat> vert3 = (Triad<Apfloat>) face.fetch(3);

            // centroid of square
            Triad<Apfloat> centroid = new Triad<>(
                    vert0.fetch(0).add(vert1.fetch(0)).add(vert2.fetch(0)).add(vert3.fetch(0)).divide(FOUR),
                    vert0.fetch(1).add(vert1.fetch(1)).add(vert2.fetch(1)).add(vert3.fetch(1)).divide(FOUR),
                    vert0.fetch(2).add(vert1.fetch(2)).add(vert2.fetch(2)).add(vert3.fetch(2)).divide(FOUR)
            );

            // given point p and vertA, calculate vector from centroid -> p:
            // m = p - centroid = (p_x - centroid_x, p_x - centroid_y, p_x -centroid_z)
            Triad<Apfloat> m = VectorMath.subs(point_cart, centroid);

            // compute dot product of m and face_norm:
            // d = face_norm ⋅ m
            //   = face_norm_x*m_x + face_norm_y*m_y + face_norm_z*m_z
            //   = face_norm_x*(p_x-vertA_x)
            //   + face_norm_y*(p_y-vertA_y)
            //   + face_norm_z*(p_z-vertA_z)
            Triad<Apfloat> norm = (Triad<Apfloat>) face_norms_sqr.fetch(i);
            Apfloat d = VectorMath.dot_prod(norm, m);

            // point outside if dot product > 0
            if (d.compareTo(ZERO) > 0) {
                return false;
            }
        }
        return true;
    }
}


