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

    #pragma omp parallel
    {
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
}

constexpr void offset_momentum(std::vector<Particle>& particles, int n)
{
    #pragma omp parallel
    {
        double px = 0.0;
        double py = 0.0;
        double pz = 0.0;
        for(int i=0;i<n;i++){
            px += particles[i].vel[0] * particles[i].m;
            py += particles[i].vel[1] * particles[i].m;
            pz += particles[i].vel[2] * particles[i].m;
        }
        particles[0].vel[0] = -px / SOLAR_MASS;
        particles[0].vel[1] = -py / SOLAR_MASS;
        particles[0].vel[2] = -pz / SOLAR_MASS;
    }
}

constexpr void advance(std::vector<Particle>& particles, double dt, int n)
{
    #pragma omp parallel
    {
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(j!=i){
                    double dx = particles[i].pos[0] - particles[j].pos[0];
                    double dy = particles[i].pos[1] - particles[j].pos[1];
                    double dz = particles[i].pos[2] - particles[j].pos[2];

                    double dSquared = dx * dx + dy * dy + dz * dz;
                    double mag = dt / (dSquared * std::sqrt(dSquared));

                    particles[i].vel[0] -= dx * particles[j].m * mag;
                    particles[i].vel[1] -= dy * particles[j].m * mag;
                    particles[i].vel[2] -= dz * particles[j].m * mag;

                    particles[j].vel[0] += dx * particles[i].m * mag;
                    particles[j].vel[1] += dy * particles[i].m * mag;
                    particles[j].vel[2] += dz * particles[i].m * mag;
                }
            }
        }
    }

    #pragma omp parallel
    {
        for(int i=0;i<n;i++){
            particles[i].pos[0] += dt * particles[i].vel[0];
            particles[i].pos[1] += dt * particles[i].vel[1];
            particles[i].pos[2] += dt * particles[i].vel[2];
        }
    }
}

double energy(const std::vector<Particle>& particles, int n)
{
    double e = 0.0;
    #pragma omp parallel
    {
        for(int i=0;i<n;i++){
            const Particle& iParticles = particles[i];
            e += 0.5 * iParticles.m * (iParticles.vel[0] * iParticles.vel[0] + iParticles.vel[1] * iParticles.vel[1] + iParticles.vel[2] * iParticles.vel[2]);
            for(int j=0;j<n;j++){
                if(i!=j){
                    double dx = iParticles.pos[0] - particles[j].pos[0];
                    double dy = iParticles.pos[1] - particles[j].pos[1];
                    double dz = iParticles.pos[2] - particles[j].pos[2];

                    double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
                    e -= (iParticles.m * particles[j].m) / distance;
                }
            }
        }
    }
    return e;
}

int main(int argc, char* argv[])
{
    const auto n = std::atoi(argv[1]);

    std::vector<Particle> particles = circular_orbits(n);
    offset_momentum(particles, n);

    std::printf("%e\n", energy(particles, n));

    for (size_t i = 0; i < n; ++i)
        advance(particles, 0.01, n);

    std::printf("%e\n", energy(particles, n));
}

//compile: g++ -std=c++11 -O2 -Wall nbody9.cpp 
//exec: ./nbody9.exe 500