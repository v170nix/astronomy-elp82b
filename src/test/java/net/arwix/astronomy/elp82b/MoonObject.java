package net.arwix.astronomy.elp82b;

import net.arwix.astronomy.core.Epoch;
import net.arwix.astronomy.core.coordinates.EclipticCoordinates;
import net.arwix.astronomy.core.vector.Vector;
import org.jetbrains.annotations.NotNull;

public class MoonObject implements EclipticCoordinates<String> {


    @NotNull
    @Override
    public Epoch getEpoch() {
        return Epoch.J2000;
    }

    @Override
    public String getIdObject() {
        return "Moon";
    }

    @NotNull
    @Override
    public Vector getEclipticCoordinates(double T) {
        return MoonElp82b.getElp82bCoordinates(T);
    }
}
