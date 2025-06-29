import com.oson.tuple.Polyad;

public abstract class Shape {
    protected final UnitCell unit_cell;
    protected final AtomicCoordinates coordinates;
    protected final String file_name;
    protected final String structure_name;
    protected final String structure_index;
    protected final int precision;
    protected final String lattice_constant;
    protected final String radius;

    public Shape(
            String radius,
            LatticeType lattice_type,
            int precision,
            Polyad<Atom> basis,
            String lattice_constant,
            String file_name,
            String structure_name,
            String structure_index
    ) {
        this.structure_index = structure_index;
        this.file_name = file_name;
        this.structure_name = structure_name;
        this.precision = precision;
        this.lattice_constant = lattice_constant;
        if (lattice_type == LatticeType.FCC) {
            this.coordinates = new FccCoordinates(
                    lattice_constant,
                    radius,
                    precision
            );
            this.unit_cell = new FccUnitCell(
                    lattice_constant,
                    precision,
                    basis.fetch(0).getElement(),
                    basis.fetch(0).getFormalChargeInt(),
                    basis.fetch(0).getRadius(),
                    basis.fetch(1).getElement(),
                    basis.fetch(1).getFormalChargeInt(),
                    basis.fetch(1).getRadius(),
                    basis.fetch(2).getElement(),
                    basis.fetch(2).getFormalChargeInt(),
                    basis.fetch(2).getRadius(),
                    basis.fetch(3).getElement(),
                    basis.fetch(3).getFormalChargeInt(),
                    basis.fetch(3).getRadius()
            );
            this.radius = radius;
        } else {
            throw new IllegalArgumentException("Illegal lattice type!");
        }
    }

    public UnitCell getUnitCell() {
        return this.unit_cell;
    }
    public AtomicCoordinates getCoordinates() {
        return this.coordinates;
    }
    public Shape getThis() {
        return this;
    }
    public String getStructureName() {
        return this.structure_name;
    }
    public String getStructureIndex() {
        return this.structure_index;
    }
    public String getLatticeConstant() {
        return this.lattice_constant;
    }
    public String getFileName() {
        return this.file_name;
    }
    public String getRadius() {
        return this.radius;
    }
    public abstract void build();
}

