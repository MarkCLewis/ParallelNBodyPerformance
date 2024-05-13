module nbodyimut

using Printf

#Imut should be done and its running about 17s on my mac 
# Constants

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

# Calculate acceleration due to gravity - these two functions indirectly accounts for energy, mome. 
# calculating the acc due to grav between two bodies  
function calculate_acceleration(body1::Body, body2::Body)
    # calcs diff in pos 
    dx = body1.x - body2.x
    dy = body1.y - body2.y
    dz = body1.z - body2.z
    # calcs dis squared
    dsq = dx^2 + dy^2 + dz^2
    # calc mag of grav 
    mag = body2.m / (dsq * √dsq)
    # calc acc components 
    ax = -dx * mag
    ay = -dy * mag
    az = -dz * mag
    ax, ay, az
end

function update_vel(body::Body, ax::Float64, ay::Float64, az::Float64, dt::Float64)
    vx_new = body.vx + ax * dt
    vy_new = body.vy + ay * dt
    vz_new = body.vz + az * dt
    Body(body.x, body.y, body.z, vx_new, vy_new, vz_new, body.m)
end

function update_pos(body::Body, dt::Float64)
    x_new = body.x + body.vx * dt
    y_new = body.y + body.vy * dt
    z_new = body.z + body.vz * dt
    Body(x_new, y_new, z_new, body.vx, body.vy, body.vz, body.m)
end

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
        push!(bods, temp) #is this pushing a reference to temp in which case is everything going to explode and how do i fix that
    end
    bods
    #print(bods)
end

# Planets: Sun, Jupiter, Saturn, Uranus, Neptune
function nbody(steps::Int64, numBodies::Int64)
    bodies = circular_orbits(numBodies)
    bodies[2] = Body(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1e-7)

#    @printf("%.9f\n", energy(bodies))
    
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
            bodies[j] = update_vel(bodies[j], ax, ay, az, 0.001)
        end
        Threads.@threads for j in 1:length(bodies) 
            bodies[j] = update_pos(bodies[j], 0.001)
        end
    end
#    @printf("%.9f\n", energy(bodies))
end

if !isinteractive()
    nbody(parse(Int, ARGS[1]), parse(Int, ARGS[2]))
end

end # module
