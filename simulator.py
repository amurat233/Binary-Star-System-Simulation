import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

import functions
from celestialBody import celestialBody
from centreOfMass import centreOfMass

# Gravitational constant of the universe.
G = 500
# Number of time_steps you want to replicate.
simulation_steps = 40000
# The length of each time step, used to speed up/slow down time.
# Lower time_step and increase simulation_steps for more detail.
time_step = 0.01
# The number of times to interpolate the data to make points appear
# more like a line.
times_interp = 0
# The number of time steps over which to calculate the area inscribed
# by the planet to the centre of mass.
time_inscribing = 1
get_stats = True
c = 100
plot_observer = True

theta = np.pi/2
distance_from_stars = 1000
observer_position = np.array([np.cos(theta - np.pi/2)*distance_from_stars,np.sin(theta - np.pi/2)*distance_from_stars])

# Creates star objects
star1 = celestialBody(10.0, 5.0, np.array([10.0,0.0]), np.array([0.0, 200.0]), G, time_step)
star2 = celestialBody(30.0, 10.0, np.array([0.0,0.0]), np.array([0.0, 0.0]), G, time_step)
# An array with all the stars in it. Used to calculate the gracitational forces.
all_stars = [star1, star2]

# Creates a centre of mass object
com =  centreOfMass(all_stars)

# A list that turns the focus function on/off and decides what it should focus on.
focus_on = [True, com]

luminosities = []
#################################################################
# SIMULATION                                                    #
#################################################################
for i in range(simulation_steps):
    # Updates the velocity of all the stars
    for star in all_stars:
        star.update_velocity(all_stars)
    
    #Focuses all the stars if focus_on[0] is true.
    if focus_on[0]:
        functions.focus(all_stars, focus_on[1])
    
    # Updates the position of all the stars
    for star in all_stars:
        star.update_position()

    # Moves COM
    com.update_position()
    com.update_velocity()
    luminosities.append(functions.calculate_luminosity(all_stars, observer_position))
#################################################################

#################################################################
# GETS AREA INSCRIBED, ECCENTRICITY, REDSHIFT, AND BRIGHTNESS   #
#################################################################
if get_stats:
    stars_deviation = []
    stars_redshift = [[]]*len(all_stars)
    for i, star in enumerate(all_stars):
        star_deviation = functions.get_area_deviation(star, com, time_inscribing)
        stars_deviation.append(star_deviation)
        star_redshift = functions.get_redshifts(observer_position, star, c)
        stars_redshift[i] = star_redshift.tolist()
        star.get_apep(com)

    for i in range(len(all_stars)):
        print(f"STAR {i+1}:")
        print(f"Maximum deviation in area inscribed: {stars_deviation[i]}")
        print(f"Eccentricity: {all_stars[i].e} \n")

#################################################################

#################################################################
# PLOTTING STAR PATH                                            #
#################################################################
stars_xyv = []
for star in all_stars:
    star_xyv = []
    for j in range(3):
        star_xyv.append(functions.interpolate([star.history[i][j] for i in range(len(star.history))], times_interp))
    stars_xyv.append(star_xyv)
for j in range(len(all_stars)-1):
    for i in range(3):
        stars_xyv[0][i].extend(stars_xyv[j+1][i])
stars_xyv = stars_xyv[0]


cm = plt.cm.get_cmap('plasma')

plt.scatter(com.position[0], com.position[1], color = "r")
plt.annotate("COM", (com.position[0], com.position[1]))

if plot_observer:
    plt.scatter(observer_position[0], observer_position[1], color = "r")
    plt.annotate("Observer", (observer_position[0], observer_position[1]))

starsx, starsy, starsv = stars_xyv
plt.scatter(starsx, starsy, c=starsv,  edgecolor='none', marker = ".", cmap = cm)

plt.axis('scaled')
cbar = plt.colorbar(label = "Speed")
plt.show()

plt.plot(stars_redshift[0])
plt.plot(stars_redshift[1])
plt.show()

time_steps = list(range(len(luminosities)))

plt.plot(luminosities)
plt.show()
#################################################################
