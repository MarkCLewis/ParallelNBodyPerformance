// Rust #3 - https://benchmarksgame-team.pages.debian.net/benchmarksgame/program/nbody-rust-3.html
//
// The Computer Language Benchmarks Game
// https://salsa.debian.org/benchmarksgame-team/benchmarksgame/
//
// Contributed by Alex Drozhak
//

use std::ops::{Add, Sub, Mul, AddAssign, SubAssign};
use std::f64::consts::PI;

use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use rayon::prelude::*;

const N_BODIES: usize = 5;

#[derive(Clone, Copy)]
struct Vec3(pub f64, pub f64, pub f64);

impl Vec3 {
    fn new() -> Self {
        Vec3(0.0, 0.0, 0.0)
    }

    fn norm(&self) -> f64 {
        self.squared_norm().sqrt()
    }

    fn squared_norm(&self) -> f64 {
        self.0 * self.0 + self.1 * self.1 + self.2 * self.2
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Vec3(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2)
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Vec3(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2)
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        *self = Vec3(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2);
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = Vec3(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2);
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Vec3(self.0 * rhs, self.1 * rhs, self.2 * rhs)
    }
}

#[derive(Clone, Copy)]
struct Particle {
    pos: Vec3,
    vel: Vec3,
    mass: f64,
}

fn circular_orbits(n: usize) -> Vec<Particle> {
    let mut particle_buf = vec![];
    particle_buf.push(Particle {
    pos: Vec3(0.0, 0.0, 0.0),
    vel: Vec3(0.0, 0.0, 0.0),
    mass: 1.0,
    });

    for i in 0..n {
        let d = 0.1 + ((i as f64) * 5.0 / (n as f64));
        let v = f64::sqrt(1.0 / d);
        let theta = fastrand::f64() * 2.0 * PI;
        let x = d * f64::cos(theta);
        let y = d * f64::sin(theta);
        let vx = -v * f64::sin(theta);
        let vy = v * f64::cos(theta);
        particle_buf.push(Particle {
            pos: Vec3(x, y, 0.0),
            vel: Vec3(vx, vy, 0.0),
            mass: 1e-14,
        });
    }
    particle_buf
}

struct NBSystem {
    planets: Vec<Particle>,
}

impl NBSystem {
    fn new(n: usize) -> Self {
        NBSystem {
            planets: circular_orbits(n)
        }
    }

    fn offset_momentum(&mut self) {
        let p = self.planets.iter_mut().fold(Vec3::new(), |v, b| v + b.vel * b.mass);
        self.planets[0].vel = p * (-1.0 / self.planets[0].mass);
    }

    fn energy(&self) -> f64 {
        let mut e = 0.0;
        let mut bodies = self.planets.iter();
        while let Some(bi) = bodies.next() {
            e +=
                bi.vel.squared_norm() * bi.mass / 2.0 -
                bi.mass *
                bodies.clone().map(|bj| bj.mass / (bi.pos - bj.pos).norm()).fold(0.0, |a, b| a + b);
        }
        e
    }

    fn advance(&mut self, dt: f64) {
        let planets = &mut self.planets;
       planets.par_iter_mut().enumerate().for_each(|(i, pi)| {
            for j in 0..N_BODIES {
                if j != i {
                    let dp = planets[i].pos - planets[j].pos;

                    let distance = dp.norm();
                    let mag = dt / (distance * distance * distance);
                    let massj = planets[j].mass;

                    planets[i].vel -= dp * massj * mag;
                }
            }
            let dpos = planets[i].vel * dt;
            planets[i].pos += dpos;
        })
    }
}

fn main() {
    let steps = std::env::args_os()
        .nth(1)
        .and_then(|s| s.into_string().ok())
        .and_then(|n| n.parse().ok())
        .unwrap_or(1000);
    let n = std::env::args_os()
        .nth(2)
        .and_then(|s| s.into_string().ok())
        .and_then(|n| n.parse().ok())
        .unwrap_or(1000);
    let mut system = NBSystem::new(n);
    system.offset_momentum();
    println!("{:.9}", system.energy());
    for _ in 0..steps {
        system.advance(0.01);
    }
    println!("{:.9}", system.energy());
}
