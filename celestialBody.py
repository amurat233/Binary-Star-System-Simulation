import numpy as np

class celestialBody:
    def __init__(self, mass, radius, initial_velocity, initial_position, G, time_step):
        k = 1 # The constant of luminosity
        self.mass = float(mass)
        # The brightness of the star, proportional to the 3.5th power of the mass.
        self.luminosity = k * self.mass**3.5 
        # Radius of the star, needed when calculating the luminosity of the system.
        self.radius = float(radius)
        # Stars positional properties.
        self.initial_velocity = initial_velocity
        self.velocity = initial_velocity
        self.position = initial_position
        # Constants of the universe.
        self.G = G
        self.time_step = time_step
        # Stores the star's position and velocity at every frame.
        self.history = []
        
    # Creates a unit vector of a vector v.
    def normalize(self, v):
        norm = np.linalg.norm(v)
        if norm == 0: 
            return v
        return v / norm

    # Changes the velocity of the star based on the force acting on it
    # by all the other stars.
    def update_velocity(self, all_bodies):
        # Loops over all the celestial bosies in the simulation
        for body in all_bodies:
            # If the body is not itself
            if body != self:
                # Calculates the acceleration from the gravitational force exerted on it.
                rsquared = (self.position - body.position).dot(self.position - body.position)
                force_dir = self.normalize(self.position - body.position) # The direction of the force
                acceleration = force_dir * self.G * body.mass / rsquared # The acceleration as a vector
                # Updates the velocity using the acceleration.
                self.velocity -= acceleration * self.time_step
    
    # Updates the position by using the velocity.
    def update_position(self):
        # Adds the current position to the history
        self.history.append([self.position[0], self.position[1], np.sqrt(self.velocity.dot(self.velocity)), self.velocity])
        # Updates the position based on the current velocity
        self.position += self.velocity * self.time_step

    # Calculates the area of a traingle with vertices A,B, and C.
    def area_triangle(self, A, B, C):
        return np.absolute((A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]))/2)

    # Calculates the area incribed by the star in a certain time for the entirety
    # of the star's lifetime.
    def area_inscribed(self, time, com):
        areas = []
        curr = 0
        nex = 1
        while True:
            area = 0
            for i in range(time):
                curr_pos = self.history[curr]
                try:
                    nex_pos = self.history[nex]
                except:
                    return areas
                curr += 1
                nex += 1
                area += self.area_triangle(com.position, curr_pos, nex_pos)
            areas.append(area)

    # Gets the distance to apoapsis and periapsis and uses it to calculate eccentricity.
    def get_apep(self, com):
        positions = [np.array(position[:2]) - com.position for position in self.history]
        distances = [np.linalg.norm(position) for position in positions]
        max_d_index = np.argmax(distances)
        min_d_index = np.argmin(distances)
        self.ap_distance = distances[max_d_index]
        self.pe_distance = distances[min_d_index]
        self.e = (self.ap_distance - self.pe_distance)/(self.ap_distance + self.pe_distance)

