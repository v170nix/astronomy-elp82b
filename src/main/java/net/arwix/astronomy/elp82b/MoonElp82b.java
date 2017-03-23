package net.arwix.astronomy.elp82b;


import net.arwix.astronomy.core.vector.RectangularVector;
import net.arwix.astronomy.core.vector.Vector;
import net.arwix.astronomy.core.vector.VectorType;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;

import static java.lang.Math.*;

abstract public class MoonElp82b {

    private static
    final int[] elp82b_max_lambda_factor = {
            10,
            6,
            8,
            6,
            2,
            8,
            4,
            6,
            5,
            17,
            61,
            65,
            71,
            14,
            12,
            4,
            4,
    };

    private static
    final double[] elp82b_constants = {
            0.000000000000000000e+00,
            0.000000000000000000e+00,
            3.850005584468033630e+05,
            -1.100000000000000039e-04,
            0.000000000000000000e+00,
            2.059999999999999765e-03,
            0.000000000000000000e+00,
            0.000000000000000000e+00,
            0.000000000000000000e+00,
    };

    private static final List<Double> elp82b_coefficients = new ArrayList<Double>(37514 * 2);

    private static final List<Short> elp82b_instructions = new ArrayList<Short>(125000);

    static {
        try {
            // Initializing arrays from Elp82b.data file
            InputStream stream = MoonElp82b.class.getClassLoader().getResourceAsStream("Elp82b.data");
            LineNumberReader reader = new LineNumberReader(new InputStreamReader(stream));
            String line = reader.readLine();
            if (!line.startsWith("elp82b_coefficients")) {
                throw new IOException("File Elp82b.data should contain elp82b_coefficients at line " + reader.getLineNumber());
            }
            while (line != null) {
                line = reader.readLine();
                if (line.contains("}")) {
                    break;
                } else {
                    StringTokenizer st = new StringTokenizer(line, " \t,");
                    elp82b_coefficients.add(Double.parseDouble(st.nextToken()));
                    elp82b_coefficients.add(Double.parseDouble(st.nextToken()));
                }
            }

            line = reader.readLine();
            if (!line.startsWith("elp82b_instructions")) {
                throw new IOException("File Elp82b.data should contain elp82b_instructions at line " + reader.getLineNumber());
            }
            while (line != null) {
                line = reader.readLine();
                if (line.contains("}")) {
                    break;
                } else {
                    StringTokenizer st = new StringTokenizer(line, " \t,");
                    while (st.hasMoreTokens()) {
                        String token = st.nextToken();
                        if (token.startsWith("0x"))
                            elp82b_instructions.add(Short.parseShort(token.substring(2), 16));
                        else if (token.contains("+16*"))
                            elp82b_instructions.add((short) (Short.parseShort(token.substring(0, 1)) + 16 * Short.parseShort(token.substring(5))));
                        else
                            elp82b_instructions.add(Short.parseShort(token));
                    }
                }
            }
      //      saveToFile("elp82b_inst", elp82b_instructions);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }



    /* Delaunay's arguments */
    private static
    final double[] del = {
            (1732559343.73604 - 129597742.2758) * (PI / (180 * 3600)),
            (129597742.2758 - 1161.2283) * (PI / (180 * 3600)),
            (1732559343.73604 - 14643420.2632) * (PI / (180 * 3600)),
            (1732559343.73604 - -6967919.3622) * (PI / (180 * 3600)),
            (-5.8883 - -0.0202) * (PI / (180 * 3600)),
            (-0.0202 - 0.5327) * (PI / (180 * 3600)),
            (-5.8883 - -38.2776) * (PI / (180 * 3600)),
            (-5.8883 - 6.3622) * (PI / (180 * 3600)),
            (0.006604 - 9e-6) * (PI / (180 * 3600)),
            (9e-6 - -1.38e-4) * (PI / (180 * 3600)),
            (0.006604 - -0.045047) * (PI / (180 * 3600)),
            (0.006604 - 0.007625) * (PI / (180 * 3600)),
            (-3.169e-5 - 1.5e-7) * (PI / (180 * 3600)),
            (1.5e-7 - 0.0) * (PI / (180 * 3600)),
            (-3.169e-5 - 2.1301e-4) * (PI / (180 * 3600)),
            (-3.169e-5 - -3.586e-5) * (PI / (180 * 3600))
    };

    /* Precession */
    private static final double[] zeta = {
            (1732559343.73604 + 5029.0966) * (PI / (180 * 3600))//w[3]+preces
    };

    /* Planetary arguments */
    private static final double[] p = {
            538101628.68898 * (PI / (180 * 3600)),
            210664136.43355 * (PI / (180 * 3600)),
            129597742.2758 * (PI / (180 * 3600)),
            68905077.59284 * (PI / (180 * 3600)),
            10925660.42861 * (PI / (180 * 3600)),
            4399609.65932 * (PI / (180 * 3600)),
            1542481.19393 * (PI / (180 * 3600)),
            786550.32074 * (PI / (180 * 3600))
    };

    /* Polynom for incrementing r1 */
    private static final double[] w = {
            218 * 3600 + 18 * 60 + 59.95571,
            1732559343.73604,
            -5.8883,
            0.006604,
            -3.169e-5,
    };


    private static final double a0_div_ath_times_au =
            384747.9806448954 / (384747.9806743165 * 149597870.691);

    /* Polynoms for transformation matrix */
    private static final double p1 = 1.0180391e-5;
    private static final double p2 = 4.7020439e-7;
    private static final double p3 = -5.417367e-10;
    private static final double p4 = -2.507948e-12;
    private static final double p5 = 4.63486e-15;
    private static final double q1 = -1.13469002e-4;
    private static final double q2 = 1.2372674e-7;
    private static final double q3 = 1.265417e-9;
    private static final double q4 = -1.371808e-12;
    private static final double q5 = -3.20334e-15;

    static public Vector getElp82bCoor(double t) {
        RectangularVector vector = (RectangularVector) getCoor(t).getVectorOfType(VectorType.RECTANGULAR);
//        final double rh = vector.z * cos(vector.y);
//        final double x3 = vector.z * sin(vector.y);
//        final double x1 = rh * cos(vector.x);
//        final double x2 = rh * sin(vector.x);

        double pw = t * (p1 + t * (p2 + t * (p3 + t * (p4 + t * p5))));
        double qw = t * (q1 + t * (q2 + t * (q3 + t * (q4 + t * q5))));
        final double pwq = pw * pw;
        final double qwq = qw * qw;
        final double pwqw = 2.0 * pw * qw;
        final double pw2 = 1.0 - 2.0 * pwq;
        final double qw2 = 1.0 - 2.0 * qwq;
        final double ra = 2.0 * sqrt(1.0 - pwq - qwq);
        pw *= ra;
        qw *= ra;
        return new RectangularVector(
                pw2 * vector.x + pwqw * vector.y + pw * vector.z,
                pwqw * vector.x + qw2 * vector.y - qw * vector.z,
                -pw * vector.x + qw * vector.y + (pw2 + qw2 - 1.0) * vector.z);
    }

    static public Vector getCoor(double t) {
    //    final double t = (jd - 2451545.0) / 36525.0;
        double[] lambda = new double[17];
        int i, k;
        for (i = 0; i < 4; i++) {
            lambda[i] = 0.0;
            for (k = 3; k >= 0; k--) {
                lambda[i] += del[k * 4 + i];
                lambda[i] *= t;
            }
            lambda[5 + i] = del[i] * t;
        }
        lambda[4] = zeta[0] * t;
        for (i = 0; i < 8; i++) {
            lambda[9 + i] = p[i] * t;
        }

        List<Double> cosSinLambda = prepareLambdaArray(elp82b_max_lambda_factor, lambda);
        double[] accu = new double[elp82b_constants.length];
        System.arraycopy(elp82b_constants, 0, accu, 0, elp82b_constants.length);
        double[] stack = new double[17 * 2];
        accumulateTerms(elp82b_instructions, elp82b_coefficients, cosSinLambda, accu, stack);

        final double r1 = (accu[0] + w[0] + t * (accu[3] + w[1] + t * (accu[6] + w[2] + t * (w[3] + t * w[4])))) * (PI / (180.0 * 3600.0));
        final double r2 = (accu[1] + t * (accu[4] + t * accu[7])) * (PI / (180.0 * 3600.0));
        final double r3 = (accu[2] + t * (accu[5] + t * accu[8])) * a0_div_ath_times_au;

                final double rh = r3 * cos(r2);
        final double x3 = r3 * sin(r2);
        final double x1 = rh * cos(r1);
        final double x2 = rh * sin(r1);

        return new RectangularVector(x1, x2, x3);

//        final double rh = r3 * cos(r2);
//        final double x3 = r3 * sin(r2);
//        final double x1 = rh * cos(r1);
//        final double x2 = rh * sin(r1);
//
//        double pw = t * (p1 + t * (p2 + t * (p3 + t * (p4 + t * p5))));
//        double qw = t * (q1 + t * (q2 + t * (q3 + t * (q4 + t * q5))));
//        final double pwq = pw * pw;
//        final double qwq = qw * qw;
//        final double pwqw = 2.0 * pw * qw;
//        final double pw2 = 1.0 - 2.0 * pwq;
//        final double qw2 = 1.0 - 2.0 * qwq;
//        final double ra = 2.0 * sqrt(1.0 - pwq - qwq);
//        pw *= ra;
//        qw *= ra;
//
//        // VSOP87 coordinates:
//
//        // VSOP87 coordinates:
//        double x = pw2 * x1 + pwqw * x2 + pw * x3;
//        double y = pwqw * x1 + qw2 * x2 - qw * x3;
//        double z = -pw * x1 + qw * x2 + (pw2 + qw2 - 1.0) * x3;
//
////        System.out.println(x * 149597870.7000);
////        System.out.println(y * 149597870.7000);
////        System.out.println(z * 149597870.7000);
//
//        return new RectangularVector(
//                pw2 * x1 + pwqw * x2 + pw * x3,
//                pwqw * x1 + qw2 * x2 - qw * x3,
//                -pw * x1 + qw * x2 + (pw2 + qw2 - 1.0) * x3);

//        printf("Moon: %f  %22.15f %22.15f %22.15f\n",
//               jd,xyz[0],xyz[1],xyz[2]);
    }

    static List<Double> prepareLambdaArray(int maxLambdaFactor[],
                                           double lambda[]) {
        List<Double> result = new ArrayList<Double>(203 * 4);
        /* initialize result:
           (cos,sin)(1*lambda[0]),(cos,sin)(-1*lambda[0]),(cos,sin)(2*lambda[0]),...
           (cos,sin)(1*lambda[1]),(cos,sin)(-1*lambda[1]),(cos,sin)(2*lambda[1]),...
        */
        int cosSinLambdaP = 0;
        int i;
        for (i = 0; i < 17; i++) {
            int maxFactor = maxLambdaFactor[i];
            double cosLambda = cos(lambda[i]);
            double sinLambda = sin(lambda[i]);
            result.add(cosLambda);
            result.add(sinLambda);
            result.add(cosLambda);
            result.add(-sinLambda);
            int m;
            double cosm0;
            double cosm1;
            double sinm0;
            double sinm1;
            for (m = 2; m <= maxFactor; m++) {
                //cslp += 4;
                /* addition theorem:
                   cos(m*l) = cos(m0*l+m1*l) = cos(m0*l)*cos(m1*l)-sin(m0*l)*sin(m1*l)
                   sin(m*l) = sin(m0*l+m1*l) = cos(m0*l)*sin(m1*l)+sin(m0*l)*cos(m1*l)
                */
                int m0 = ((((m) >> 1) - 1) << 2);
                int m1 = ((((m + 1) >> 1) - 1) << 2);
                cosm0 = result.get(cosSinLambdaP + m0);
                cosm1 = result.get(cosSinLambdaP + m1);
                sinm0 = result.get(cosSinLambdaP + m0 + 1);
                sinm1 = result.get(cosSinLambdaP + m1 + 1);
                cosLambda = cosm0 * cosm1 - sinm0 * sinm1;
                sinLambda = cosm0 * sinm1 + sinm0 * cosm1;
                result.add(cosLambda);
                result.add(sinLambda);
                result.add(cosLambda);
                result.add(-sinLambda);
            }
            cosSinLambdaP += (maxFactor << 2);
        }

        return result;
    }



    static class IndexHolder {
        int spIdx = 0;
        int lambdaIndex = 0;
        int termCount = 0;
    }

    static void accumulateTerms(List<Short> instructions,
                                List<Double> coefficients,
                                List<Double> cosSinLambda,
                                double accu[],
                                double[] sp) {
        /* Accumulates the series given in instructions/coefficients.
           The argument of the series is cos_sin_lambda which has
           been initialized using prepareLambdaArray().
           accu is the output accumulator.
           sp must point to a memory area holding 2*nr_of_lambdas double.
           area must be supplied by the caller, it will be destroyed during the
           calculation.
        */
        Iterator<Short> instructionsIterator = instructions.iterator();
        Iterator<Double> coefficientsIterator = coefficients.iterator();
        IndexHolder idx = new IndexHolder();
        sp[0] = 1.0;
        sp[1] = 0.0;
            for (; ;) {
                idx.termCount = instructionsIterator.next();
                if (idx.termCount < 0xFE) {
                    idx.lambdaIndex = ((idx.termCount & 15) << 8) | (instructionsIterator.next());
                    idx.termCount >>= 4;
                    calculateNewArg(cosSinLambda, accu, sp,
                            instructionsIterator, coefficientsIterator,
                            idx);
                } else {
                    if (idx.termCount == 0xFF) break;
                    /* pop argument from the stack */
                    idx.spIdx -= 2;
                }
            }
    }

    private static void calculateNewArg(List<Double> cosSinLambda,
                                        double[] accu,
                                        double[] sp,
                                        Iterator<Short> instructionsIterator,
                                        Iterator<Double> coefficientsIterator,
                                        IndexHolder idx) {
        /* calculate new argument and push it on the stack */
        double cosLambda = cosSinLambda.get(idx.lambdaIndex << 1);
        double sinLambda = cosSinLambda.get((idx.lambdaIndex << 1) + 1);
        sp[idx.spIdx + 2] = cosLambda * sp[idx.spIdx] - sinLambda * sp[idx.spIdx + 1];
        sp[idx.spIdx + 3] = cosLambda * sp[idx.spIdx + 1] + sinLambda * sp[idx.spIdx];
        idx.spIdx += 2;
        while (--idx.termCount >= 0) {
            Double coef1 = coefficientsIterator.next();
            Double coef2 = coefficientsIterator.next();
            accu[instructionsIterator.next()] += (coef1 * sp[idx.spIdx] + coef2 * sp[idx.spIdx + 1]);
        }
    }


}
