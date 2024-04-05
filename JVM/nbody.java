// Java #5 - https://benchmarksgame-team.pages.debian.net/benchmarksgame/program/nbody-java-5.html
/* The Computer Language Benchmarks Game
   https://salsa.debian.org/benchmarksgame-team/benchmarksgame/

   contributed by Mark C. Lewis
   double[] instead of Object[] by Han Kai
*/

import java.util.random.RandomGenerator;


public final class nbody {
   public static void main(String[] args) {
      int steps = Integer.parseInt(args[0]);
      int numBodies = Integer.parseInt(args[1]);

      NBodySystem bodies = new NBodySystem(numBodies);
      System.out.printf("%e\n", bodies.energy());
      for (int i = 0; i < steps; ++i)
         bodies.advance(0.01);
      System.out.printf("%e\n", bodies.energy());
   }

   final static class NBodySystem {
      private static final int BODY_SIZE = 8;

      private static final int x = 0;
      private static final int y = 1;
      private static final int z = 2;
      private static final int vx = 3;
      private static final int vy = 4;
      private static final int vz = 5;
      private static final int mass = 6;

      private final double[] _bodies;
      private final int numBodies;

      public NBodySystem(int numBodies) {
         this.numBodies = numBodies;
         _bodies = circularOrbits(numBodies);
      }

      public void advance(double dt) {
         final double[] bodies = _bodies;

         for (int i = 0; i < numBodies; ++i) {
            final int offset = BODY_SIZE * i;

            for (int j = i + 1; j < numBodies; ++j) {
               final int ioffset = offset;
               final int joffset = BODY_SIZE * j;

               final double dx = bodies[ioffset + x] - bodies[joffset + x];
               final double dy = bodies[ioffset + y] - bodies[joffset + y];
               final double dz = bodies[ioffset + z] - bodies[joffset + z];

               final double dSquared = dx * dx + dy * dy + dz * dz;
               final double distance = Math.sqrt(dSquared);
               final double mag = dt / (dSquared * distance);

               final double jmass = bodies[joffset + mass];

               bodies[ioffset + vx] -= dx * jmass * mag;
               bodies[ioffset + vy] -= dy * jmass * mag;
               bodies[ioffset + vz] -= dz * jmass * mag;

               final double imass = bodies[ioffset + mass];
               bodies[joffset + vx] += dx * imass * mag;
               bodies[joffset + vy] += dy * imass * mag;
               bodies[joffset + vz] += dz * imass * mag;
            }
         }

         for (int i = 0; i < numBodies; ++i) {
            final int ioffset = BODY_SIZE * i;

            bodies[ioffset + x] += dt * bodies[ioffset + vx];
            bodies[ioffset + y] += dt * bodies[ioffset + vy];
            bodies[ioffset + z] += dt * bodies[ioffset + vz];
         }
      }

      public double energy() {
         final double[] bodies = _bodies;

         double dx, dy, dz, distance;
         double e = 0.0;

         for (int i = 0; i < numBodies; ++i) {
            final int offset = BODY_SIZE * i;

            final double ivx = bodies[offset + vx];
            final double ivy = bodies[offset + vy];
            final double ivz = bodies[offset + vz];
            final double imass = bodies[offset + mass];

            e += 0.5 * imass * (ivx * ivx + ivy * ivy + ivz * ivz);

            for (int j = i + 1; j < numBodies; ++j) {
               final int ioffset = offset;
               final int joffset = BODY_SIZE * j;

               final double ix = bodies[ioffset + x];
               final double iy = bodies[ioffset + y];
               final double iz = bodies[ioffset + z];

               dx = ix - bodies[joffset + x];
               dy = iy - bodies[joffset + y];
               dz = iz - bodies[joffset + z];

               distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
               e -= (imass * bodies[joffset + mass]) / distance;
            }
         }

         return e;
      }

      static void setBody(int index, double[] sys, double x, double y, double z, double vx, double vy, double vz, double radius, double mass) {
         int base = index * BODY_SIZE;
         sys[base + 0] = x;
         sys[base + 1] = y;
         sys[base + 2] = z;
         sys[base + 3] = vx;
         sys[base + 4] = vy;
         sys[base + 5] = vz;
         sys[base + 6] = radius;
         sys[base + 7] = mass;
      }

      
      static double[] circularOrbits(int n) {
         var rand = RandomGenerator.getDefault();
         var system = new double[(n+1) * BODY_SIZE];
         setBody(0, system, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00465047, 1.0 );
         
         for (int i = 0; i < n; ++i) {
            double d = 0.1 + (i * 5.0 / n);
            double v = Math.sqrt(1.0 / d);
            double theta = rand.nextDouble(0.0, 6.28);
            double x = d * Math.cos(theta);
            double y = d * Math.sin(theta);
            double vx = -v * Math.sin(theta);
            double vy = v * Math.cos(theta);
            setBody(i+1, system,
               x, y, 0.0,
               vx, vy, 0.0,
               1e-14,
               1e-7
            );
         }
         return system;
      }
   }

}
