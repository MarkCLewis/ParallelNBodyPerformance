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

# kicks 
function advance!(bodies::Array{Body}, dt::Float64)
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
    Threads.@threads for i in 1:length(bodies)
        bodi = bodies[i]
        bodi.x += dt * bodi.vx
        bodi.y += dt * bodi.vy
        bodi.z += dt * bodi.vz
    end
end

# energy
function energy(bodies::Array{Body})::Float64
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

function circular_orbits(n::Int64)::Array{Body} 
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
function nbody(steps::Int64, numBodies::Int64)
    bods = circular_orbits(numBodies)

    # do advancing stuff
#    @printf("%e\n", energy(bods))
    for i = 1:steps
        advance!(bods, 0.001)
    end
#    @printf("%e\n", energy(bods))
end

if !isinteractive()
    nbody(parse(Int, ARGS[1]), parse(Int, ARGS[2]))
end



end
