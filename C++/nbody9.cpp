// C++ #9 - https://benchmarksgame-team.pages.debian.net/benchmarksgame/program/nbody-gpp-9.html
/* The Computer Language Benchmarks Game
   https://salsa.debian.org/benchmarksgame-team/benchmarksgame/

   contributed by Martin Jambrek
   based off the Java #2 program contributed by Mark C. Lewis and modified slightly by Chad Whipkey
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <omp.h>
#include <array>
#include <vector>
#include <random>


static constexpr double PI = 3.141592653589793;
static constexpr double SOLAR_MASS = 4 * PI * PI;
static constexpr double DAYS_PER_YEAR = 365.24;

struct alignas(32) Particle{
    std::array<double, 3> pos; 
    std::array<double, 3> vel;
    double m;

    constexpr Particle(std::array<double, 3> pos, std::array<double, 3> vel, double m):
        pos(pos),
        vel(vel),
        m(m)
    {
    }
};

std::vector<Particle> circular_orbits(int n){
    std::vector<Particle> particle_buf;
    particle_buf.push_back(Particle{
        {0.0,0.0,0.0},
        {0.0,0.0,0.0},
        1.0
    });

    for(int i=0;i<n;i++){
        double d = 0.1 + (i * 5.0 / n);
        double v = std::sqrt(1.0 / d);
        double theta = rand() * 6.28;
        double x = d * std::cos(theta);
        double y = d * std::sin(theta);
        double vx = -v * std::sin(theta);
        double vy = v * std::cos(theta);
        particle_buf.push_back(Particle{
            {x,y,0.0},
            {vx,vy,0.0},
            1e-14,
        });
    }
    return particle_buf;
}

void advance(std::vector<Particle>& particles, double dt)
{
    #pragma omp parallel for
    for(int i=0;i<particles.size();i++){
        for(int j=0;j<particles.size();j++){
            if(j!=i){
                double dx = particles[i].pos[0] - particles[j].pos[0];
                double dy = particles[i].pos[1] - particles[j].pos[1];
                double dz = particles[i].pos[2] - particles[j].pos[2];

                double dSquared = dx * dx + dy * dy + dz * dz;
                double mag = dt / (dSquared * std::sqrt(dSquared));

                particles[i].vel[0] -= dx * particles[j].m * mag;
                particles[i].vel[1] -= dy * particles[j].m * mag;
                particles[i].vel[2] -= dz * particles[j].m * mag;
            }
        }
        particles[i].pos[0] += dt * particles[i].vel[0];
        particles[i].pos[1] += dt * particles[i].vel[1];
        particles[i].pos[2] += dt * particles[i].vel[2];
    }
}

double energy(const std::vector<Particle>& particles)
{
    double e = 0.0;
    #pragma omp parallel for reduction(+ : e)
    for(int i=0;i<particles.size();i++){
        const Particle& iParticles = particles[i];
        e += 0.5 * iParticles.m * (iParticles.vel[0] * iParticles.vel[0] + iParticles.vel[1] * iParticles.vel[1] + iParticles.vel[2] * iParticles.vel[2]);
        for(int j=0;j<particles.size();j++){
            if(i!=j){
                double dx = iParticles.pos[0] - particles[j].pos[0];
                double dy = iParticles.pos[1] - particles[j].pos[1];
                double dz = iParticles.pos[2] - particles[j].pos[2];

                double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
                e += -(iParticles.m * particles[j].m) / distance;
            }
        }
    }
    return e;
}

int main(int argc, char* argv[])
{
    const auto steps = std::atoi(argv[1]);
    const auto n = std::atoi(argv[2]);

    omp_set_num_threads(std::atoi(argv[3]));

    std::vector<Particle> particles = circular_orbits(n);


    std::printf("%e\n", energy(particles));

    for (size_t i = 0; i < steps; ++i) {
        advance(particles, 0.001);
    }

    std::printf("%e\n", energy(particles));
}

//compile: g++ -std=c++11 -O2 -Wall nbody9.cpp 
//exec: ./nbody9.exe 500
