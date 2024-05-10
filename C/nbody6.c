// C #6 - https://benchmarksgame-team.pages.debian.net/benchmarksgame/program/nbody-gcc-6.html
/* The Computer Language Benchmarks Game
 * https://salsa.debian.org/benchmarksgame-team/benchmarksgame/
 *
 * contributed by Christoph Bauer
 * modified by Danny Angelo Carminati Grein
 *  
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define pi 3.141592653589793
#define solar_mass (4 * pi * pi)
#define days_per_year 365.24

struct planet {
  double x, y, z;
  double vx, vy, vz;
  double mass;
};

struct planet* circular_orbits(int n){
    struct planet* particle_buf = malloc((n+1)*sizeof(struct planet));
    particle_buf[0].x = 0.0;
    particle_buf[0].y = 0.0;
    particle_buf[0].z = 0.0;
    particle_buf[0].vx = 0.0;
    particle_buf[0].vy = 0.0;
    particle_buf[0].vz = 0.0;
    particle_buf[0].mass = 1.0;

    for(int i=1;i<=n;i++){
        double d = 0.1 + (i * 5.0 / n);
        double v = sqrt(1.0 / d);
        double theta = rand() * 6.28;
        double x = d * cos(theta);
        double y = d * sin(theta);
        double vx = -v * sin(theta);
        double vy = v * cos(theta);
        particle_buf[i].x = x;
        particle_buf[i].y = y;
        particle_buf[i].z = 0.0;
        particle_buf[i].vx = vx;
        particle_buf[i].vy = vy;
        particle_buf[i].vz = 0.0;
        particle_buf[i].mass = 1.0e-14;
    }
    return particle_buf;
}

void advance(int nbodies, struct planet * bodies, double dt)
{
  int i, j;

  #pragma omp parallel for
  for (i = 0; i < nbodies; i++) {
    struct planet * b = &(bodies[i]);
    for (j = 0; j < nbodies; j++) {
      if (i != j) {
        struct planet * b2 = &(bodies[j]);
        double dx = b->x - b2->x;
        double dy = b->y - b2->y;
        double dz = b->z - b2->z;
        double distanced = dx * dx + dy * dy + dz * dz;
        double distance = sqrt(distanced);
        double mag = dt / (distanced * distance);
        b->vx -= dx * b2->mass * mag;
        b->vy -= dy * b2->mass * mag;
        b->vz -= dz * b2->mass * mag;
      }
    }
    b->x += dt * b->vx;
    b->y += dt * b->vy;
    b->z += dt * b->vz;
  }
}

double energy(int nbodies, struct planet * bodies)
{
  double e;
  int i, j;

  e = 0.0;
  #pragma omp parallel for reduction(+ : e)
  for (i = 0; i < nbodies; i++) {
    struct planet * b = &(bodies[i]);
    e += 0.5 * b->mass * (b->vx * b->vx + b->vy * b->vy + b->vz * b->vz);
    for (j = i + 1; j < nbodies; j++) {
      struct planet * b2 = &(bodies[j]);
      double dx = b->x - b2->x;
      double dy = b->y - b2->y;
      double dz = b->z - b2->z;
      double distance = sqrt(dx * dx + dy * dy + dz * dz);
      e += -(b->mass * b2->mass) / distance;
    }
  }
  return e;
}

int main(int argc, char ** argv)
{
  int steps = atoi(argv[1]);
  int n = atoi(argv[2]);
  int i;

  omp_set_num_threads(atoi(argv[3]));

  struct planet* bodies = circular_orbits(n);

  printf ("%e\n", energy(n+1, bodies));
  for (i = 1; i <= steps; i++) {
    advance(n+1, bodies, 0.001);
  }
  printf ("%e\n", energy(n+1, bodies));
  free(bodies);
  return 0;
}
