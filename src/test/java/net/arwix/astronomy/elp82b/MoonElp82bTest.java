package net.arwix.astronomy.elp82b;

import com.sun.org.apache.xerces.internal.impl.xpath.regex.Match;
import net.arwix.astronomy.core.AstronomyMatrix;
import net.arwix.astronomy.core.Constants;
import net.arwix.astronomy.core.calendar.CalendarExtensions;
import net.arwix.astronomy.core.vector.*;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.Calendar;
import java.util.Locale;
import java.util.TimeZone;
import java.util.concurrent.TimeUnit;

import static java.lang.Math.sin;
import static org.junit.Assert.*;

public class MoonElp82bTest {
    @Before
    public void setUp() throws Exception {

    }

    @Test
    public void getCoor() throws Exception {

        double au = 149597870.700;
        double t = (2469000.5 - 2451545.0) / 36525.0;
//        System.out.println(t);
//        RectangularVector cR = (RectangularVector) MoonElp82b.getCoordinates(t).getVectorOfType(VectorType.RECTANGULAR);
//        System.out.println(Arrays.toString(cR.toArray()));
        RectangularVector vectorR = (RectangularVector) MoonElp82b.getElp82bCoordinates(t).getVectorOfType(VectorType.RECTANGULAR);
        System.out.println(Arrays.toString(vectorR.toArray()));
        assertEquals("2469000.5 xE2000", -361602.98537, vectorR.x * au , 1e-5);
        assertEquals("2469000.5 yE2000", 44996.9951, vectorR.y * au , 1e-5);
        assertEquals("2469000.5 zE2000", -30696.65316, vectorR.z * au , 1e-5);

        Calendar calendar = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
        calendar.set(2047, Calendar.OCTOBER, 17, 0,0,6);
  //      calendar.set(1989, Calendar.JANUARY, 1, 0,0,0);
        t = CalendarExtensions.getJT(calendar, true);
        System.out.println(t);
        System.out.println(CalendarExtensions.getMJD(calendar) );
        System.out.println("TDB-UT " + CalendarExtensions.getDeltaT(calendar, TimeUnit.SECONDS));
//        System.out.println(CalendarExtensions.getMJD(calendar));
//        vectorR = (RectangularVector) MoonElp82b.getCoor(t).getVectorOfType(VectorType.RECTANGULAR);
//        assertEquals("2469000.5 xE2000", -361602.98481, vectorR.x * au , 1e-5);
//        assertEquals("2469000.5 yE2000", 44996.99625, vectorR.y * au , 1e-5);
//        assertEquals("2469000.5 zE2000", -30696.65152, vectorR.z * au , 1e-5);

    //    double eps = Math.toRadians(23.0 + 26.0 / 60.0 + 21.40883/ 60.0 / 60.0);
    //    Matrix var10000 = Matrix.Companion.transpose(Matrix.Companion.invoke(Matrix.Axis.X, eps));


        Vector vector = MoonElp82b.getElp82bCoordinates(t);
        double dT = vector.normalize() / Constants.C_Light / 36525.0;
        System.out.println(dT * 36525.0 * 24.0 * 60.0 );
        double innerT = t - dT;
        vector = MoonElp82b.getElp82bCoordinates(innerT);
   //    SphericalVector eclVector = (SphericalVector) vector.getVectorOfType(VectorType.SPHERICAL);

   //     System.out.println(eclVector.r * au);
  //      System.out.println(Math.toDegrees(eclVector.theta));
   //     System.out.println(Math.toDegrees(eclVector.phi));

//        eclVector = new SphericalVector( Math.toRadians(70.3037769), Math.toRadians(-3.1776621), eclVector.r);
//        vector = eclVector;
        double k = calendar.get(Calendar.YEAR) - 1988.0;
        double tdb = 56.059 + 0.2822 * k + 0.02223 * k * k;
        System.out.println("ydb" + tdb);



        SphericalVector equ = (SphericalVector)
//                AstronomyMatrix.INSTANCE.createTransformationCoordinates(t, AstronomyMatrix.Coordinates.ECLIPTIC, AstronomyMatrix.Coordinates.EQUATORIAL)
//                        .times(vector).getVectorOfType(VectorType.SPHERICAL);
//
                AstronomyMatrix.INSTANCE.createTransformationCoordinates(0, AstronomyMatrix.Coordinates.ECLIPTIC, AstronomyMatrix.Coordinates.EQUATORIAL)
                        .times(vector).getVectorOfType(VectorType.SPHERICAL);
//

    //    Matrix PN = AstronomyMatrix.INSTANCE.createNutation(innerT).times(AstronomyMatrix.INSTANCE.createPrecession(0, innerT, AstronomyMatrix.Coordinates.EQUATORIAL));



//
  //     equ = (SphericalVector) AstronomyMatrix.INSTANCE.createNutation(innerT).times(equ).getVectorOfType(VectorType.SPHERICAL);

  //      equ = (SphericalVector) AstronomyMatrix.INSTANCE.createPrecession(0, innerT, AstronomyMatrix.Coordinates.EQUATORIAL).times(equ).getVectorOfType(VectorType.SPHERICAL);

   //     equ = (SphericalVector) PN.times(equ).getVectorOfType(VectorType.SPHERICAL);
//        Matrix Ecl2Equ = AstronomyMatrix.INSTANCE.createTransformationCoordinates(t, AstronomyMatrix.Coordinates.ECLIPTIC, AstronomyMatrix.Coordinates.EQUATORIAL);
//        Vector ev = KeplerianOrbit.Planet.EARTH.getOrbitalPlane(t).component2().div(Constants.C_Light);
//        Vector ve = Ecl2Equ.times(ev);
//        equ = (SphericalVector) equ.plus(ve).getVectorOfType(VectorType.SPHERICAL);
//        equ = (SphericalVector) equ.div(equ.normalize()).getVectorOfType(VectorType.SPHERICAL);


        System.out.println(equ.r);
        System.out.println(printLong(equ));
        System.out.println(printLat(equ));

//        equ = (SphericalVector) equ
//                .times(AstronomyMatrix.INSTANCE.createNutation(t))
//                .getVectorOfType(VectorType.SPHERICAL);
//        System.out.println(equ.r);
//        System.out.println(printLong(equ));
//        System.out.println(printLat(equ));
    }

    static private String printLong(Vector p) {
        SphericalVector vector = (SphericalVector) p.getVectorOfType(VectorType.SPHERICAL);

        double hours = (Constants.DEG * vector.phi / 15.0);

        final int hour = (int) hours;
        final double minutes = (hours - hour) * 60.0;
        final int minute = (int) minutes;
        final double seconds = (minutes - minute) * 60.0;

        return String.format(Locale.ENGLISH, "%1$02d:%2$02d:%3$.2f", hour, minute, seconds);
    }

    static private String printLat(Vector p) {
        SphericalVector vector = (SphericalVector) p.getVectorOfType(VectorType.SPHERICAL);

        int g = (int) Math.toDegrees(vector.theta);
        double mm = (Math.toDegrees(vector.theta) - g) * 60.0;
        int m = (int) mm;
        double s = (mm - m) * 60.0;
        return String.format(Locale.ENGLISH, "%1$02d %2$02d %3$.1f", g, Math.abs(m), Math.abs(s));
    }

}