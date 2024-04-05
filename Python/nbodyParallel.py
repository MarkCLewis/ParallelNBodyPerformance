import sys 
import multiprocessing

#=======================================NOTES AND TO DO==================================================
#-Advance needs to be parallel
#-Need to figure out main and execution
#-Need to test code in general, especially concerned about math being accurate
#-Upload to repository once allowed
#========================================================================================================
class Particle:
    def __init__ (self, p, v, r, m):
        self.p = p #position(Vec3): (0.0,0.0,0.0) 
        self.v = v #velocity(Vec3): (0.0,0.0,0.0)
        self.r = r #radius(integer): 0
        self.m = m #mass(integer): 0

def combinations(l):
    result = []
    for x in range(len(l) - 1):
        ls = l[x+1:]
        for y in ls:
            result.append((l[x],y))
    return result

def two_bodies():
    two = []
    twoBodies.append(Particle((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1, 1)) #pos, vel, rad, mass
    twoBodies.append(Particle((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), 1e-4, 1e-20)) #pos, vel, rad, mass
    return twoBodies

#these are pretty much 1-1 with the rust version Dr. Lewis made
def circular_orbits():
    particle_buf = []
    particle_buf.append(Particle((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.00465047, 1.0)) #pos, vel, rad, mass

    for i in range(n):
        d = 0.1 + (i * 5.0 / n)
        v = math.sqrt(1.0 / d)
        theta = random.random() * 6.28
        x = d * math.cos(theta)
        y = d * math.sin(theta)
        vx = -v * math.sin(theta)
        vy = v * math.cos(theta)
        particle_buf.append(Particle([x, y, 0.0], [vx, vy, 0.0], 1e-7, 1e-14))

        return particle_buf #makes a bunch of particles around an origin

bodiesArray = circular_orbits() #this is whats used in all the calculations

def distance_sqr(x1, x2):
    dx = x1[0] - x2[0]
    dy = x1[1] - x2[1]
    dz = x1[2] - x2[2]
    return dx*dx + dy*dy + dz*dz

def distance(x1, x2):
    return math.sqrt(distance_sqr(x1, x2))

#This is where the parallel stuff starts
SYSTEM = list(bodiesArray)
PAIRS = combinations(SYSTEM)


def advance(dt, n, bodies=SYSTEM, pairs=PAIRS):

    for i in range(n):
        for (([x1, y1, z1], v1, m1),
             ([x2, y2, z2], v2, m2)) in pairs: #for each pair of particles...
            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2 #finds distance
            mag = dt * ((dx * dx + dy * dy + dz * dz) ** (-1.5)) #magnitude
            b1m = m1 * mag #body 1 mag
            b2m = m2 * mag #body 2 mag
            v1[0] -= dx * b2m #modify b1's velocity
            v1[1] -= dy * b2m #^
            v1[2] -= dz * b2m #^
            v2[0] += dx * b1m #modify b2's velocity
            v2[1] += dy * b1m #^
            v2[2] += dz * b1m #^
        for (r, [vx, vy, vz], m) in bodies: #for each body...
            r[0] += dt * vx #modify position based on same vel vector
            r[1] += dt * vy #^
            r[2] += dt * vz #^

def modify_i(pairs=PAIRS):
    for (([x1, y1, z1], v1, m1),
             ([x2, y2, z2], v2, m2)) in pairs: #for each pair of particles...
            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2 #finds distance
            mag = dt * ((dx * dx + dy * dy + dz * dz) ** (-1.5)) #magnitude
            b2m = m2 * mag #body 2 mag
            v1[0] -= dx * b2m #modify i's velocity
            v1[1] -= dy * b2m #^
            v1[2] -= dz * b2m #^
            i_vels = v1
    return i_vels

def modify_j(pairs=PAIRS): 
    for (([x1, y1, z1], v1, m1),
             ([x2, y2, z2], v2, m2)) in pairs: #for each pair of particles...
            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2 #finds distance
            mag = dt * ((dx * dx + dy * dy + dz * dz) ** (-1.5)) #magnitude
            b1m = m1 * mag #body 1 mag
            v2[0] -= dx * b1m #modify j's velocity
            v2[1] -= dy * b1m #^
            v2[2] -= dz * b1m #^
            j_vels = v2
    return j_vels #PROBLEM

def advance_parallel(dt, n, bodies=SYSTEM, pairs=PAIRS):
    mp_pool = multiprocessing.Pool()

    i_velocities = mp_pool.map(modify_i, pairs) #PROBLEM
    j_velocities = mp_pool.map(modify_j, pairs) #PROBLEM




#Sequential Energy
def report_energy(bodies=SYSTEM, pairs=PAIRS, e=0.0):

   for (((x1, y1, z1), v1, m1), #for position, velocity(vector?), and mass of 2 particles in pairs
         ((x2, y2, z2), v2, m2)) in pairs:
        dx = x1 - x2 #distance x
        dy = y1 - y2 #distance y
        dz = z1 - z2 #distance z
        e -= (m1 * m2) / ((dx * dx + dy * dy + dz * dz) ** 0.5) #gravitational potential energy, Euclidian distance
        
    for (r, [vx, vy, vz], m) in bodies: #step 2
        e += m * (vx * vx + vy * vy + vz * vz) / 2. #kinetic energy
    print("%.9f" % e)

def potential_energy(pairs=PAIRS):

    for (((x1, y1, z1), v1, m1), #for position, velocity(vector?), and mass of 2 particles in pairs
         ((x2, y2, z2), v2, m2)) in pairs:
        dx = x1 - x2 #distance x
        dy = y1 - y2 #distance y
        dz = z1 - z2 #distance z
        return (m1 * m2) / ((dx * dx + dy * dy + dz * dz) ** 0.5) #gravitational potential energy, Euclidian distance
        
def kinetic_energy(bodies=SYSTEM)
    for (r, [vx, vy, vz], m) in bodies:
        return m * (vx * vx + vy * vy + vz * vz) / 2. #kinetic energy
#Parallel Energy
def report_energy_parallel(bodies=SYSTEM, pairs=PAIRS):
    mp_pool = multiprocessing.Pool()

    potEnergy = mp_pool.map(potential_energy, pairs) #calcualtes all in parallel
    kinEnergy = mp.pool.map(kinetic_energy, bodies) #calculates all in parallel

    sumEnergies = sum(kinEnergy) - sum(potEnergy) #kinetic - potential (sums)

    print("%.9f" % sumEnergies)

def offset_momentum(ref, bodies=SYSTEM, px=0.0, py=0.0, pz=0.0):

    for (r, [vx, vy, vz], m) in bodies:
        px -= vx * m
        py -= vy * m
        pz -= vz * m
    (r, v, m) = ref
    v[0] = px / m
    v[1] = py / m
    v[2] = pz / m

"""def main(n, ref='sun'):
    offset_momentum(bodiesArray[ref])
    report_energy()
    advance(0.01, n)
    report_energy()

if __name__ == '__main__':
    main(int(sys.argv[1]))"""