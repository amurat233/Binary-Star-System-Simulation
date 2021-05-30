import numpy as np
import copy

# A function that sets the velocities of all the bodies in it relative to the "focus"
# object.
def focus(bodies, focus):
    for body in bodies:
        if body != focus:
            body.velocity = body.velocity - focus.velocity
    focus.velocity = np.array([0.0,0.0])


# Function that prints the maximum deviation in area inscribed by
# a line from the planet to the centre of mass (Kepler's second law)
def get_area_deviation(star, com, time_inscribing):
    areas = star.area_inscribed(time_inscribing, com)
    mean_area = sum(areas)/len(areas)
    max_deviaiton = max(np.absolute(1-(max(areas)/mean_area)), np.absolute(1-(min(areas)/mean_area)))
    return max_deviaiton


# Linear interpolate function that moves adds more points in between
# in order to make the points appear like lines
def interpolate(datapoints, number_times):
    # Doesn't do anything if asked to interpolate 0 times.
    if number_times == 0:
        return datapoints
    # Interpolates number_times times.
    for i in range(number_times):
        # Creates an output list where the interpolated list goes.
        output = []
        # Loops over every point in the list.
        for i in range(len(datapoints)-1):
            # Appends the point itself to the output list.
            output.append(datapoints[i])
            # Appends the average between the point and the next
            # point to the output list.
            output.append((datapoints[i]+datapoints[i+1])/2)
        # Updates the datapoints list so that interpolation can
        # happen again.
        datapoints = output.copy()
    return output

# Calculates the redshift given the spped of light and the recessional
# velocity.
def calculate_redshift(c,v):
    numerator = np.sqrt(c+v)
    denominator = np.sqrt(c-v)
    redshift = numerator/denominator - 1
    return redshift

# Given the position of the observer and information about the star it calculates
# the recessional velocity and uses that to find the redshift.
def get_redshifts(observer_position, star, c = 100):
    star_positions = [[star.history[i][0],star.history[i][1]] for i in range(len(star.history))]
    star_velocities = [star.history[i][3] for i in range(len(star.history))]

    redshifts = np.zeros(len(star_positions))

    for i, position in enumerate(star_positions):
        relative_position = np.array(position) - observer_position
        mag_position = np.linalg.norm(relative_position)
        star_velocity = star_velocities[i]
        recessional_velocity = star_velocity.dot(relative_position)/mag_position
        redshifts[i] = calculate_redshift(c, recessional_velocity)

    return redshifts

# Finds A given a,b,c using cosine rule
def cosine_rule(a,b,c):
    numerator = b**2 + c**2 - a**2
    denominator = 2 * b * c
    A = np.arccos(numerator/denominator)
    return A

# Calculates the brightness of the star system at a given time. BROKEN
def calculate_luminosity(all_stars, observer_position):
    # Finds which star is closest and arranges them in increasing order of distance.
    # Only stores the stars relative position, radius, and luminosity.
    distances = [np.linalg.norm(star.position - observer_position) for star in all_stars]
    sorted_stars = []
    for distance, star in sorted(zip(distances, all_stars)):
        sorted_stars.append([star.position - observer_position, star.radius, star.luminosity])
    distances = sorted(distances)

    # "Resizes" the further star so that it is the same distance as the closer star
    d1_d2 = distances[0]/distances[1] # closer_distance/further_distance
    sorted_stars[1][1] = sorted_stars[1][1]*d1_d2
    sorted_stars[1][0] = sorted_stars[1][0]*d1_d2

    # The distance d between the closer star and the resized star
    d = np.linalg.norm(sorted_stars[0][0]-sorted_stars[1][0])

    # Checks if the closer star is covering the further star,
    # if it isn't then the luminosity doesn't change.
    if d >= sorted_stars[0][1] + sorted_stars[1][1]:
        return sorted_stars[0][2] + sorted_stars[1][2]
    
    # Other values needed to calculate coverage
    r1 = sorted_stars[0][1]
    r2 = sorted_stars[1][1]

    if r1>=r2:
        if d <= r1 - r2:
            return sorted_stars[0][2]
    else:
        if d <= r2-r1:
            proportion_uncovered_2 = 1 - (np.pi * r1**2)/(np.pi * r2**2)
            return sorted_stars[0][2] + sorted_stars[1][2]*proportion_uncovered_2

    # Calculates what proportion of the further star is covered by the
    # closer star.
    theta1 = cosine_rule(r2,r1,d)
    theta2 = cosine_rule(r1,r2,d)

    area_triangle_1 = 0.5*(r1**2)*np.sin(2*theta1)
    area_slice_1 = theta1 * r1**2
    area_section_1 = area_slice_1 - area_triangle_1

    area_triangle_1 = 0.5 * (r2**2) * np.sin(2*theta2)
    area_outer_slice_2 = r2 ** 2 * (np.pi - theta2)
    area_uncovered = area_outer_slice_2 + area_triangle_1 - area_section_1
    proportion_uncovered_2 = area_uncovered/(np.pi * r2**2)

    return sorted_stars[0][2] + sorted_stars[1][2]*proportion_uncovered_2

def calculate_luminosity_simple(all_stars, observer_position):
    # Finds which star is closest and arranges them in increasing order of distance.
    # Only stores the stars relative position, radius, and luminosity.
    distances = [np.linalg.norm(star.position - observer_position) for star in all_stars]
    sorted_stars = []
    for distance, star in sorted(zip(distances, all_stars)):
        sorted_stars.append([star.position - observer_position, star.radius, star.luminosity])
    distances = sorted(distances)

    # "Resizes" the further star so that it is the same distance as the closer star
    d1_d2 = distances[0]/distances[1] # closer_distance/further_distance
    sorted_stars[1][1] = sorted_stars[1][1]*d1_d2
    sorted_stars[1][0] = sorted_stars[1][0]*d1_d2
    # The distance d between the closer star and the resized star
    d = np.linalg.norm(sorted_stars[0][0]-sorted_stars[1][0])

    # Checks if the closer star is covering the further star,
    # if it isn't then the luminosity doesn't change.
    # Otherwise return the luminosity of the closer uncovered star.
    if d >= sorted_stars[0][1] + sorted_stars[1][1]:
        return sorted_stars[0][2] + sorted_stars[1][2]
    
    # Other values needed to calculate coverage
    r1 = sorted_stars[0][1]
    r2 = sorted_stars[1][1]

    if r1>=r2:
        return sorted_stars[0][2]
    else:
        proportion_uncovered_2 = 1 - (np.pi * r1**2)/(np.pi * r2**2)
        return sorted_stars[0][2] + sorted_stars[1][2]*proportion_uncovered_2