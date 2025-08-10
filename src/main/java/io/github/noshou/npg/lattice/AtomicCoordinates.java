package io.github.noshou.npg.lattice;

import io.github.noshou.tuple.Triad;
import org.apfloat.Apfloat;

public interface AtomicCoordinates {
    Triad<Apfloat> getPosition();
}