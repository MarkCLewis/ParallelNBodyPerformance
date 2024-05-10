module nbodymutparallel
using Printf
# NBody Mutable
# Audrey Tollett

# Constants
const SOLAR_MASS = 4 * pi * pi
const DAYS_PER_YEAR = 365.24
const n = 50000000 #n bodies but hard coded xd

# Body Class
mutable struct Body
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    m::Float64
end

function getDim(bod, dim)
    if dim == 0
        bod.x
    elseif dim == 1
        bod.y 
    else
        bod.x
    end
end

# account for momentum 
function offsetMomentum!(s, bodies)
    px = py = pz = 0.0
    for b in bodies
        px -= b.vx * b.m
        py -= b.vy * b.m
        pz -= b.vz * b.m
    end
    s.vx = px / SOLAR_MASS
    s.vy = py / SOLAR_MASS
    s.vz = pz / SOLAR_MASS
end

# kicks 
function advance!(bodies, dt)
    Threads.@threads for i in 1:length(bodies)
        bodi = bodies[i]
        for j in 1:length(bodies)
            if i != j
                bodj = bodies[j]
                dx = bodi.x - bodj.x
                dy = bodi.y - bodj.y
                dz = bodi.z - bodj.z

                dsq = dx^2 + dy^2 + dz^2
                mag = dt / (dsq*sqrt(dsq))

                bodi.vx -= dx * bodj.m * mag
                bodi.vy -= dy * bodj.m * mag
                bodi.vz -= dz * bodj.m * mag
            end
        end
    end
    for i=1:length(bodies)
        bodi = bodies[i]
        bodi.x += dt * bodi.vx
        bodi.y += dt * bodi.vy
        bodi.z += dt * bodi.vz
    end
end

# energy
function energy(bodies)
    e = 0.0
    for i in 1:length(bodies)
        bodi = bodies[i]
        e += 0.5 * bodi.m * (bodi.vx^2 + bodi.vy^2 + bodi.vz^2)
        for j in (i+1):length(bodies)
            bodj = bodies[j]
            d = sqrt(((bodi.x-bodj.x)^2+(bodi.y-bodj.y)^2+(bodi.z-bodj.z)^2))
            e -= bodi.m * bodies[j].m / d
        end
    end
    e
end

function circular_orbits(n) 
    first = Body(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    bods = [first]
    for i in 1:n
        d = .1 + (i * 5.0 / n)
        v = sqrt(1.0 / d)
        theta = rand(Float64)*2*pi
        x = d * cos(theta)
        y = d * sin(theta)
        vx = -v * sin(theta)
        vy = v * cos(theta)
        temp = Body(x, y, 0.0, vx, vy, 0, 1.0e-7)
        push!(bods, temp)
    end
    bods
    #print(bods)
end

# planets sun - jupiter - saturn - uranus - neptune 
function nbody(steps, numBodies)
    sun = Body(0,0,0,0,0,0, SOLAR_MASS)

    jupiter = Body( 4.84143144246472090e+0,                   # x
                   -1.16032004402742839e+0,                   # y
                   -1.03622044471123109e-1,                   # z
                    1.66007664274403694e-3 * DAYS_PER_YEAR,   # vx
                    7.69901118419740425e-3 * DAYS_PER_YEAR,   # vy
                   -6.90460016972063023e-5 * DAYS_PER_YEAR,   # vz
                    9.54791938424326609e-4 * SOLAR_MASS)      # mass

    saturn = Body( 8.34336671824457987e+0,
                   4.12479856412430479e+0,
                  -4.03523417114321381e-1,
                  -2.76742510726862411e-3 * DAYS_PER_YEAR,
                   4.99852801234917238e-3 * DAYS_PER_YEAR,
                   2.30417297573763929e-5 * DAYS_PER_YEAR,
                   2.85885980666130812e-4 * SOLAR_MASS)

    uranus = Body( 1.28943695621391310e+1,
                  -1.51111514016986312e+1,
                  -2.23307578892655734e-1,
                   2.96460137564761618e-3 * DAYS_PER_YEAR,
                   2.37847173959480950e-3 * DAYS_PER_YEAR,
                  -2.96589568540237556e-5 * DAYS_PER_YEAR,
                   4.36624404335156298e-5 * SOLAR_MASS)

    neptune = Body( 1.53796971148509165e+1,
                   -2.59193146099879641e+1,
                    1.79258772950371181e-1,
                    2.68067772490389322e-3 * DAYS_PER_YEAR,
                    1.62824170038242295e-3 * DAYS_PER_YEAR,
                   -9.51592254519715870e-5 * DAYS_PER_YEAR,
                    5.15138902046611451e-5 * SOLAR_MASS)

    #bods = [jupiter, saturn, uranus, neptune]
    bods = circular_orbits(numBodies)
    #offsetMomentum!(sun, bods)
    # pushfirst!(bods, sun)

    # do advancing stuff
    @printf("%.9f\n", energy(bods))
    for i = 1:steps
        advance!(bods, 0.01)
    end
    @printf("%.9f\n", energy(bods))
end

if !isinteractive()
    nbody(parse(Int, ARGS[1]), parse(Int, ARGS[2]))
end



end
