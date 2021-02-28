# Two-Line Element Reader
# Ryo Takatori
# 11/17/2020

def tle_reader(filename):
    # Open TLE text file
    with open(filename, 'r') as f:
        tle_lines_raw = f.readlines()
    # Preallocate list
    tle_title = []
    tle_1 = []
    tle_2 = []
    tle = []

    # Sort TLE
    tle_title.append(" ".join(tle_lines_raw[0].split())) # title
    tle.append(tle_title)
    # Line 1
    tle_1.append(int(tle_lines_raw[1][0])) # line number
    tle_1.append(int(tle_lines_raw[1][2:7])) # satellite catalog number
    tle_1.append(tle_lines_raw[1][7]) # classification
    # International designator
    tle_1.append(int(tle_lines_raw[1][9:11])) # last two digits of launch year
    tle_1.append(int(tle_lines_raw[1][11:14])) # launch number of the year
    tle_1.append(" ".join(tle_lines_raw[1][14:17].split())) # piece of launch
    #
    tle_1.append(int(tle_lines_raw[1][18:20])) # epoch year (last two digits of launch year)
    tle_1.append(float(tle_lines_raw[1][20:32])) # epoch (day of the year and fractional portion of the day)
    tle_1.append(float(tle_lines_raw[1][33:43])) # first derivative of mean motion (ballistic coefficient)
    tle_1.append(float("0."+" ".join(tle_lines_raw[1][44:50].split()))*10**int(tle_lines_raw[1][50:52])) # second derivative of mean motion
    tle_1.append(float("0."+" ".join(tle_lines_raw[1][53:59].split()))*10**int(tle_lines_raw[1][59:61])) # BSTAR
    tle_1.append(int(tle_lines_raw[1][62])) # ephemeris type
    tle_1.append(int(tle_lines_raw[1][64:68])) # element set number
    tle_1.append(int(tle_lines_raw[1][68])) # check sum
    tle.append(tle_1)
    # Line 2
    tle_2.append(int(tle_lines_raw[2][0])) # line number
    tle_2.append(int(tle_lines_raw[2][2:7])) # satellite catalog number
    tle_2.append(float(tle_lines_raw[2][8:16])) # inclination [deg]
    tle_2.append(float(tle_lines_raw[2][17:25])) # right ascension of ascending node [deg]
    tle_2.append(float("0."+tle_lines_raw[2][26:33])) # eccentricity
    tle_2.append(float(tle_lines_raw[2][34:42])) # argument of perigee [deg]
    tle_2.append(float(tle_lines_raw[2][43:51])) # mean anomaly [deg]
    tle_2.append(float(tle_lines_raw[2][52:63])) # mean motion [rev/day]
    tle_2.append(int(tle_lines_raw[2][63:68])) # revolution number at epoch [rev]
    tle_2.append(int(tle_lines_raw[2][68])) # check sum
    tle.append(tle_2)
    return tle