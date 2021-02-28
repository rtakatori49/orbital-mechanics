# Planet Specifications
# Ryo Takatori
# 11/18/2020

# Data based from Curtis Textbook
class Planet():
    def __init__(self, list):
        self.radius = list[0] # radius [km]
        self.mass = list[1] # mass [km]
        self.sidereal_rotation_period = list[2] # sidereal rotation period [hour]
        self.inclination_to_equator = list[3] # inclination of equator to orbit plane [deg]
        self.a = list[4] # semi-major axis of orbit around Sun
        self.ecc = list[5] # eccentricity of orbit around Sun
        self.inclination_to_ecliptic = list[6] # inclination of orbit to the ecliptic plane [deg]
        self.orbit_sidereal_period = list[7] # orbit sidereal period [day]
        self.mu = list[8] # gravitational parameter [km^3/s^2]
        self.soi_radius = list[9] # sphere of influence radius

# Solar System
Sun = Planet([696000, 1.989e30, 25.38*24, 7.25, None, None, None, None, 132712000000, None])
Mercury = Planet([2440, 330.2e21, 58.65*24, 0.01, 57.91e6, 0.2056, 7.00, 87.97, 22030, 112000])
Venus = Planet([6052, 4.869e24, 243*24, 177.4, 108.2e6, 0.0067, 3.39, 224.7, 324900, 616000])
Earth = Planet([6378, 5.974e24, 23.9345, 23.45, 149.6e6, 0.0167, 0, 365.256, 398600, 925000])
Moon = Planet([1737, 73.48e21, 27.32*24, 6.68, 384.4e3, 0.0549, 5.145, 27.322, 4903, 66100])
Mars = Planet([3396, 641.9e21, 24.62, 25.19, 227.9e6, 0.0935, 1.850, 1.881*365, 42828, 577000])
Jupiter = Planet([71490, 1.899e27, 9.925, 3.13, 778.6e6, 0.0489, 1.304, 11.86*365, 126686000, 48200000])
Saturn = Planet([60270, 568.5e4, 10.66, 26.73, 1.433e9, 0.0565, 2.485, 29.46*365, 37931000, 54800000])
Uranus = Planet([25560, 86.83e24, 17.24, 97.77, 2.872e9, 0.0457, 0.772, 84.01*365, 5794000, 51800000])
Neptune = Planet([4760, 102.4e24, 16.11, 28.32, 4.495e9, 0.0113, 1.769, 164.8*365, 6835100, 86600000])
Pluto = Planet([1195, 12.5e21, 6.387*24, 122.5, 5.870e9, 0.2444, 17.16, 247.7*365, 830, 3080000])