import numpy as np

class centreOfMass:
    def __init__(self, bodies):
        self.position = np.array([0,0])
        self.bodies = bodies
        self.mass = sum([body.mass for body in bodies])
        self.update_velocity()
    
    def update_position(self):
        y_moment = sum([body.position[1]*body.mass for body in self.bodies])
        x_moment = sum([body.position[0]*body.mass for body in self.bodies])
        
        self.position[0] = x_moment/self.mass
        self.position[1] = y_moment/self.mass
        
    def update_velocity(self):
        momentum_y = sum([body.velocity[1] * body.mass for body in self.bodies])
        momentum_x = sum([body.velocity[0] * body.mass for body in self.bodies])
        self.velocity = np.array([momentum_x/self.mass, momentum_y/self.mass])
