import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd

import functions
from celestialBody import celestialBody
from centreOfMass import centreOfMass

# Gravitational constant of the universe.
G = 500
# Number of time_steps you want to replicate.
simulation_steps = 20000
# The length of each time step, used to speed up/slow down time.
# Lower time_step and increase simulation_steps for more detail.
time_step = 0.01
# The number of times to interpolate the data to make points appear
# more like a line.
times_interp = 0
# The number of time steps over which to calculate the area inscribed
# by the planet to the centre of mass.
time_inscribing = 2000
get_stats = True
c = 100
plot_observer = True

theta = 0
distance_from_stars = 1000
observer_position = np.array([np.cos(theta - np.pi/2)*distance_from_stars,np.sin(theta - np.pi/2)*distance_from_stars])

# Creates star objects
star1 = celestialBody(10.0, 5.0, np.array([12.0,0.0]), np.array([0.0, 200.0]), G, time_step)
star2 = celestialBody(50.0, 15.0, np.array([0.0,0.0]), np.array([0.0, 0.0]), G, time_step)
# star3 = celestialBody(20.0, 5.0, np.array([-7.0,0.0]), np.array([0.0, -800.0]), G, time_step)
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
    stars_velocities = [[]]*len(all_stars)
    for i, star in enumerate(all_stars):
        star_deviation = functions.get_area_deviation(star, com, time_inscribing)
        stars_deviation.append(star_deviation)
        star_redshift, star_velocities = functions.get_redshifts(observer_position, star, c)
        stars_redshift[i] = star_redshift.tolist()
        stars_velocities[i] = star_velocities
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

if focus_on[1] == com:
    plt.scatter(com.position[0], com.position[1], color = "r")
    plt.annotate("COM", (com.position[0], com.position[1]))

if plot_observer:
    plt.scatter(observer_position[0], observer_position[1], color = "r")
    plt.annotate("Observer", (observer_position[0], observer_position[1]))

starsx, starsy, starsv = stars_xyv
plt.scatter(starsx, starsy, c=starsv,  edgecolor='none', marker = ".", cmap = cm)

plt.axis('scaled')
# cbar = plt.colorbar(label = "Speed")
plt.show()

# cm_redshifts = plt.cm.get_cmap("coolwarm")
# plt.scatter(np.concatenate([np.arange(len(stars_redshift[0])), np.arange(len(stars_redshift[1]))]), stars_redshift[0]+stars_redshift[1], c = stars_velocities[0] + stars_velocities[1], cmap = cm_redshifts, s = 4)
# cbar = plt.colorbar(label = "Recessional Velocity")
# plt.plot(stars_redshift[0])
# plt.plot(stars_redshift[1])
# # print(max(stars_redshift[0]))
# # print(max(stars_redshift[1]))
# plt.ylabel('Redshift')
# plt.xlabel('Time Step')
# plt.show()

# df = pd.read_csv("redshift06.csv")
# df[str(theta*180/np.pi)] = stars_redshift[0]
# df.to_csv("redshift06.csv", index = False)

time_steps = list(range(len(luminosities)))

plt.plot(luminosities)
plt.show()
#################################################################
