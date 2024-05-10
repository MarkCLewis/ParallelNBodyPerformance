module nbodyimut

using Printf

#Imut should be done and its running about 17s on my mac 
# Constants
const SOLAR_MASS = 4 * pi * pi
const DAYS_PER_YEAR = 365.24
const n = 50000000 # Number of bodies

# Body struct
struct Body
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    m::Float64
end

# Constructor for the Sun
function init_sun()
    Body(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, SOLAR_MASS)
end

# Calculate acceleration due to gravity - these two functions indirectly accounts for energy, mome. 
# calculating the acc due to grav between two bodies  
function calculate_acceleration(body1, body2)
    # calcs diff in pos 
    dx = body1.x - body2.x
    dy = body1.y - body2.y
    dz = body1.z - body2.z
    # calcs dis squared
    dsq = dx^2 + dy^2 + dz^2
    # calc mag of grav 
    mag = (SOLAR_MASS / dsq^1.5) * body2.m
    # calc acc components 
    ax = -dx * mag
    ay = -dy * mag
    az = -dz * mag
    ax, ay, az
end

# Update velocity and position based on acc and time step 
function update(body::Body, ax, ay, az, dt)
    vx_new = body.vx + ax * dt
    vy_new = body.vy + ay * dt
    vz_new = body.vz + az * dt
    x_new = body.x + vx_new * dt
    y_new = body.y + vy_new * dt
    z_new = body.z + vz_new * dt
    Body(x_new, y_new, z_new, vx_new, vy_new, vz_new, body.m)
end

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
        push!(bods, temp) #is this pushing a reference to temp in which case is everything going to explode and how do i fix that
    end
    bods
    #print(bods)
end

# Planets: Sun, Jupiter, Saturn, Uranus, Neptune
function nbody(steps, numBodies)
    sun = init_sun()

    jupiter = Body( 4.84143144246472090e+0,
                   -1.16032004402742839e+0,
                   -1.03622044471123109e-1,
                    1.66007664274403694e-3 * DAYS_PER_YEAR,
                    7.69901118419740425e-3 * DAYS_PER_YEAR,
                   -6.90460016972063023e-5 * DAYS_PER_YEAR,
                    9.54791938424326609e-4 * SOLAR_MASS)

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

    #bodies = [jupiter, saturn, uranus, neptune]
    #pushfirst!(bodies, sun)
    bodies = circular_orbits(numBodies)

    @printf("%.9f\n", energy(bodies))
    
    # main sim loop 
    for i in 1:steps
        # itter over body and calc the acc # this is what you should make parallel 
        #println(Threads.nthreads())
        Threads.@threads for j in 1:length(bodies) 
            ax, ay, az = 0.0, 0.0, 0.0
            # calc acc due to grav from other bodies 
            for k in 1:length(bodies)
                if j != k
                    ax_temp, ay_temp, az_temp = calculate_acceleration(bodies[j], bodies[k])
                    ax += ax_temp
                    ay += ay_temp
                    az += az_temp
                end
            end
            # update veloc and pos of body 
            update(bodies[j], ax, ay, az, 0.01)
        end
    end
    @printf("%.9f\n", energy(bodies))
    # print final pos for all bodies 
    #for body in bodies
       # @printf("x: %.6f, y: %.6f, z: %.6f\n", body.x, body.y, body.z)
    #end
end

if !isinteractive()
    nbody(parse(Int, ARGS[1]), parse(Int, ARGS[2]))
end

end # module
