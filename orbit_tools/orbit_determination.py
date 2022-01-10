# Orbit Determination
# Ryo Takatori
# 2021-06-27

# Module
import copy
import io
import math
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import sys

from . import coordinate_transformation as c_t
from . import equation_of_motion as eom
from . import planet_data as p_d
from . import time_trans as t_t

earth = p_d.Earth()
mue = earth.mu
re = earth.radius

sys.stdout = io.TextIOWrapper(sys.stdout.detach(), encoding = 'utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.detach(), encoding = 'utf-8')

# Initial Orbit Determination

## Gauss' Method
# Initial orbit determination using Gauss' method
# Function takes user input for non-extended and extended
# 0 for non-extend
# 1 for extend

def rltau(lat, long, alt, input_datetime_list, ra_list, dec_list):
    R_temp = []
    L_temp = []
    for input_datetime, ra, dec in zip(input_datetime_list, ra_list, dec_list):
        dec = np.deg2rad(dec)
        ra = np.deg2rad(ra)
        r_temp = c_t.lla_to_eci(lat, long, alt, input_datetime)
        R_temp.append(r_temp) # Site location [km]
        l_temp = [
            np.cos(dec)*np.cos(ra),
            np.cos(dec)*np.sin(ra),
            np.sin(dec)
            ]
        L_temp.append(l_temp) # Observation unit vector
    R = np.transpose(np.array([R_temp[0], R_temp[1], R_temp[2]]))
    L = np.transpose(np.array([L_temp[0], L_temp[1], L_temp[2]]))
    # Time differences [s]
    tau1 = (input_datetime_list[0]-input_datetime_list[1]).total_seconds()
    tau3 = (input_datetime_list[2]-input_datetime_list[1]).total_seconds()
    tau = tau3 - tau1
    return R, L, tau1, tau3, tau

def gauss(lat, long, alt, input_datetime_list, ra_list, dec_list, extend=True):
    def valladosolver(L, R, tau1, tau3, tau):
        # Solve various variables M, a1, a3, a1u, a3u, d1, d2, C
        M = np.dot(np.linalg.inv(L), R)
        a1 = tau3/tau
        a3 = -tau1/tau
        a1u = a1 * (1/6) * (tau**2-tau3**2)
        a3u = a3 * (1/6) * (tau**2-tau1**2)
        d1 = M[1, 0]*a1 - M[1, 1] + M[1, 2]*a3
        d2 = M[1, 0]*a1u + M[1, 2]*a3u
        C = np.dot(L[:, 1], R[:, 1])
        # Solve 8th order polynomial
        coeff = [1, 0, -((d1**2)+(2*C*d1)+(np.linalg.norm(R[:,1])**2)), 0, \
            0, -2*mue*((C*d2)+(d1*d2)), 0, 0, -(mue**2)*(d2**2)]
        solution = np.roots(coeff)

        # Find real positive R2 magnitude and return
        norm_r2 = [sol.real for sol in [complex(sol) for sol in solution] if sol.imag == 0 and sol.real > 0]
        # Calculate u, c
        u = mue/norm_r2[0]**3
        c1 = a1 + a1u*u
        c2 = -1
        c3 = a3 + a3u*u
        c = np.array([c1, c2, c3])
        # Solve for rho magnitudes and r vector
        X = np.dot(M, -c)
        rho = []
        r = []
        for n in range(3):
            rho.append(X[n]/c[n])
            r.append(rho[n]*L[:, n]+R[:, n])
        return r, rho
    # Calculate R, L, and tau
    R, L, tau1, tau3, tau = rltau(
        lat, long, alt, input_datetime_list, ra_list, dec_list)
    # Calculate r
    r, rho = valladosolver(L, R, tau1, tau3, tau)
    # Curtis version
    # [r, rho, D0, D] = curtissolver(L, R, tau1, tau3, tau)
    r2 = r[1]
    norm_r2 = np.linalg.norm(r2)
    # Lagrange coefficient calculation for v2
    f1 = 1 - (0.5*mue*tau1**2)/norm_r2**3
    f3 = 1 - (0.5*mue*tau3**2)/norm_r2**3
    g1 = tau1 - ((1/6)*mue*tau1**3)/norm_r2**3
    g3 = tau3 - ((1/6)*mue*tau3**3)/norm_r2**3
    # v2
    v2 = (1/(f1*g3-f3*g1)) * (-f3*r[0]+f1*r[2])
    # Extended
    # Iterate to improve accuracy
    if extend:
        # Initialize to enter loop
        diff = 1000
        tol = 1e-8 # Tolerance
        rhostore = [rho] # rho values
        j = 0
        while diff > tol:
            norm_r2 = np.linalg.norm(r2) # r2 magntiude [km]
            norm_v2 = np.linalg.norm(v2) # v2 magntiude [km/s]
            alpha = (2/norm_r2) - (norm_v2**2/mue) # Reciprocal of semi-major axis ([1/km])
            v2r2 = np.dot(v2,r2)/norm_r2 # Radial component of v2 [km/s]
            # Store old Lagrange coefficient
            f1_old = f1
            f3_old = f3
            g1_old = g1
            g3_old = g3
            # Universal anomaly [sqrt(km)]
            x1 = eom.kepler_u(tau1, norm_r2, v2r2, alpha, mue)
            x3 = eom.kepler_u(tau3, norm_r2, v2r2, alpha, mue)
            # New Lagrange coefficient
            f1, g1 = eom.f_and_g(x1, tau1, norm_r2, alpha, mue)
            f3, g3 = eom.f_and_g(x3, tau3, norm_r2, alpha, mue)
            # Average Lagrange coefficient
            f1 = (f1+f1_old)/2
            f3 = (f3+f3_old)/2
            g1 = (g1+g1_old)/2
            g3 = (g3+g3_old)/2
            # Calculate u, c
            c1 = g3/(f1*g3-f3*g1)
            c3 = -g1/(f1*g3-f3*g1)
            c2 = -1
            c = np.array([c1, c2, c3])
            # Solve for rho magnitudes and r vector
            M = np.dot(np.linalg.inv(L), R)
            X = np.dot(M, -c)
            rho = []
            r = []
            for n in range(3):
                rho.append(X[n]/c[n])
                r.append(rho[n]*L[:, n]+R[:, n])
            # Update r2 and v2
            r2 = r[1]
            v2 = (1/(f1*g3-f3*g1)) * (-f3*r[0]+f1*r[2])
            # Calculate difference
            diff = abs(rhostore[j][1]-rho[1])
            # Store rho
            rhostore.append(rho)
            j += 1
    return r2, v2


# Double-R
def double_r(lat, long, alt, input_datetime_list, ra_list, dec_list):
    def dr(R, L, tau1, tau3, c1, c2, norm_r1, norm_r2):
        q1 = (-c1+np.sqrt(c1**2-4*(np.dot(R[:, 0], R[:, 0])-norm_r1**2)))/2
        q2 = (-c2+np.sqrt(c2**2-4*(np.dot(R[:, 1], R[:, 1])-norm_r2**2)))/2
        
        r1 = R[:, 0] + q1*L[:, 0]
        norm_r1 = np.linalg.norm(r1)
        r2 = R[:, 1] + q2*L[:, 1]
        norm_r2 = np.linalg.norm(r2)
        
        what = np.cross(r1, r2)/(norm_r1*norm_r2)
        q3 = np.dot(-R[:, 2], what)/np.dot(L[:, 2], what)
        r3 = R[:, 2] + q3*L[:, 2]
        norm_r3 = np.linalg.norm(r3)
        
        # difference angles
        costheta21 = np.dot(r2, r1)/(norm_r2*norm_r1)
        costheta31 = np.dot(r3, r1)/(norm_r3*norm_r1)
        costheta32 = np.dot(r3, r2)/(norm_r3*norm_r2)
        
        sintheta21 = np.sqrt(1-costheta21**2)
        sintheta31 = np.sqrt(1-costheta31**2)
        sintheta32 = np.sqrt(1-costheta32**2)
        
        theta31 = np.arccos(costheta31)
        if theta31 > np.pi:
            c1 = (norm_r2*sintheta32)/(norm_r1*sintheta31)
            c3 = (norm_r2*sintheta21)/(norm_r3*sintheta31)
            p = (c1*norm_r1+c3*norm_r3-norm_r2)/(c1+c3-1)
        elif theta31 <= np.pi:
            c1 = (norm_r1*sintheta31)/(norm_r2*sintheta32)
            c3 = (norm_r1*sintheta21)/(norm_r3*sintheta32)
            p = (c3*norm_r3-c1*norm_r2+norm_r1)/(-c1+c3+1)
        
        ecostheta1 = p/norm_r1 -1
        ecostheta2 = p/norm_r2 -1
        ecostheta3 = p/norm_r3 -1
        
        theta21 = np.arccos(costheta21)
        
        if theta21 != np.pi:
            esintheta2 = (-costheta21*ecostheta2+ecostheta1)/sintheta21
        else:
            esintheta2 = (costheta32*ecostheta2-ecostheta3)/sintheta31
        
        e2 = ecostheta2**2 + esintheta2**2
        a = p/(1-e2)
        if np.sqrt(e2) < 1:
            n = np.sqrt(mue/a**3)
            S = norm_r2/p*np.sqrt(1-e2)*esintheta2
            C = norm_r2/p*(e2+ecostheta2)
            sinE32 = norm_r3/np.sqrt(a*p)*sintheta32 - norm_r3/p*(1-costheta32)*S
            cosE32 = 1 - (norm_r2*norm_r3)/(a*p)*(1-costheta32)
            sinE21 = norm_r1/np.sqrt(a*p)*sintheta21 + norm_r1/p*(1-costheta21)*S
            cosE21 = 1 - (norm_r2*norm_r1)/(a*p)*(1-costheta21)
            M32 = np.arccos(cosE32) + 2*S*np.sin(np.arccos(cosE32)/2)**2 - C*sinE32
            M12 = -np.arccos(cosE21) + 2*S*np.sin(np.arccos(cosE21)/2)**2 + C*sinE21
            F1 = tau1 - M12/n
            F2 = tau3 - M32/n
            f = 1 - a/norm_r2*(1-cosE32)
            g = tau3 - np.sqrt(a**3/mue)*(np.arccos(cosE32)-sinE32)
        else:
            n = np.sqrt(mue/-a**3)
            Sh = norm_r2/p*np.sqrt(e2-1)*esintheta2
            Ch = norm_r2/p*(e2+ecostheta2)
            sinhF32 = norm_r3/np.sqrt(-a*p)*sintheta32 - norm_r3/p*(1-costheta32)*Sh
            F32 = np.log(sinhF32+np.sqrt(sinhF32+1))
            sinhF21 = norm_r1/np.sqrt(-a*p)*sintheta21 + norm_r1/p*(1-costheta32)*Sh
            F21 = np.log(sinhF21+np.sqrt(sinhF21**2+1))
            M32 = -F32 + 2*Sh*np.sinh(F32/2)**2 + Ch*sinhF32
            M12 = F21 + 2*Sh*np.sinh(F21/2)**2 + Ch*sinhF21
            F1 = tau1 - M12/n
            F2 = tau3 - M32/n
            f = 1 - (-a)/norm_r2*(1-np.cosh(F32))
            g = tau3 - np.sqrt((-a)**3/mue)*(F32-sinhF32)
        return F1, F2, f, g, r3 ,r2

    # Calculate R, L, and tau
    R, L, tau1, tau3, _ = rltau(
        lat, long, alt, input_datetime_list, ra_list, dec_list)
    
    # Initial guess for iteration
    norm_r1 = 2*re
    norm_r2 = 2.01*re
    c1 = np.dot(2*L[:, 0], R[:, 0])
    c2 = np.dot(2*L[:, 1], R[:, 1])
    count = 1
    error = 1
    while error > 0.001:
        F1, F2, f, g, r3, r2 = dr(R, L, tau1, tau3, c1, c2, norm_r1, norm_r2)
        dr1 = 0.005*norm_r1
        dr2 = 0.005*norm_r2
        F1r1dr1, F2r1dr1, _, _, _, _ = dr(
            R, L, tau1, tau3, c1, c2, norm_r1+dr1, norm_r2)
        F1r2dr2, F2r2dr2, _, _, _, _ = dr(
            R, L, tau1, tau3, c1, c2, norm_r1, norm_r2+dr2)
        
        dF1dr1 = (F1r1dr1-F1)/dr1
        dF2dr1 = (F2r1dr1-F2)/dr1
        dF1dr2 = (F1r2dr2-F1)/dr2
        dF2dr2 = (F2r2dr2-F2)/dr2
        
        # run the function with the dr conditions
        
        delta = dF1dr1*dF2dr2 - dF2dr1*dF1dr2
        delta1 = dF2dr2*F1 - dF1dr2*F2
        delta2 = dF1dr1*F2 - dF2dr1*F1
        dr1 = -delta1/delta
        dr2 = -delta2/delta
        
        error = (abs(dr1)+abs(dr2))/2
        norm_r1 += dr1
        norm_r2 += dr2
        count += 1
    
    # after convergence
    v2 = (r3-f*r2)/g
    return r2, v2

# Orbit Determination

# Herrick-Gibbs
def herrick_gibbs(r1, r2, r3, mu, jd1, jd2, jd3):
    dt31 = (jd3-jd1)*24*60*60
    dt32 = (jd3-jd2)*24*60*60
    dt21 = (jd2-jd1)*24*60*60
    norm_r1 = np.linalg.norm(r1)
    norm_r2 = np.linalg.norm(r2)
    norm_r3 = np.linalg.norm(r3)
    v2 = -dt32*((1/(dt21*dt31))+(mu/(12*norm_r1**3)))*r1 \
        + (dt32-dt21)*((1/(dt21*dt32))+(mu/(12*norm_r2**3)))*r2 \
        + dt21*((1/(dt32*dt31))+(mu/(12*norm_r3**3)))*r3
    return v2
# Gibbs
def gibbs(r1, r2, r3, mu, jd1=0, jd2=0, jd3=0):
    z12 = np.cross(r1, r2)
    z23 = np.cross(r2, r3)
    z31 = np.cross(r3, r1)
    norm_r1 = np.linalg.norm(r1)
    norm_r2 = np.linalg.norm(r2)
    norm_r3 = np.linalg.norm(r3)
    # Coplanar check
    try:
        alpha_cop = 90 - np.rad2deg(np.arccos(np.dot(z23, r1)
            /(np.linalg.norm(z23)*norm_r1)))
        if alpha_cop > 3:
            raise Exception
    except:
        print("Not coplanar")
        return
    # Check spacing
    try:
        alpha12 = np.rad2deg(np.arccos(np.dot(r1, r2)/(norm_r1*norm_r2)))
        alpha23 = np.rad2deg(np.arccos(np.dot(r2, r3)/(norm_r2*norm_r3)))
        if alpha12 < 1 and alpha23 < 1:
            raise Exception
        n = norm_r1*z23 + norm_r2*z31 + norm_r3*z12
        d = z12 + z23 + z31
        s = (norm_r2-norm_r3)*r1 + (norm_r3-norm_r1)*r2 + (norm_r1-norm_r2)*r3
        b = np.cross(d, r2)
        lg = np.sqrt(mu/(np.linalg.norm(n)*np.linalg.norm(d)))
        v2 = (lg/norm_r2)*b + lg*s
    except:
        print("Spacing too small. Moving into Herrick-Gibbs.")
        v2 = herrick_gibbs(r1, r2, r3, mu, jd1, jd2, jd3)
    
        
    return v2
# Linear Least Squares
def lls(x_o_i, y_o_i, print_result=True):
    n = len(x_o_i) # number of samples
    x_one = np.array(np.ones(n)).T # x full of one size of n
    h = np.array((x_one, np.array(x_o_i))).T # mapping matrix
    p = np.linalg.pinv(np.dot(h.T, h)) # covariance matrix
    x_hat = np.dot(p, np.dot(h.T, np.array(y_o_i))) # solution
    alpha, beta =  x_hat # y-intercept, slope
    ri_list = []
    j_list = []
    for x_o, y_o in zip(x_o_i, y_o_i):
        yc = alpha + beta*x_o # computed value of dependant variable
        ri = y_o - yc # residual
        ri_list.append(ri)
        j = ri**2 # cost function
        j_list.append(j)
    rms = np.sqrt((1/(n-1))*sum(j_list)) # root mean square
    sigma_alpha = np.sqrt(p[0, 0]) # standard deviation of y-intercept
    sigma_beta = np.sqrt(p[1, 1]) # standard deviation of slope
    sigma_x_hat = [sigma_alpha, sigma_beta] # standard deviation solution
    # Display
    if print_result:
        print(f"Solution: y = {alpha:.3f} + {beta:.3f}x")
        print(f"Root mean square: {rms:.3f}")
        print(f"Confidence intervals: \u03B1 = {alpha:.3f} \u00B1 {sigma_alpha:.3f}")
        print(f"Confidence intervals: \u03B2 = {beta:.3f} \u00B1 {sigma_beta:.3f}")
    return x_hat, ri_list, rms, sigma_x_hat

def find_nominal(lat, long, alt, df_raw, fun, mu):
    df_raw["datetime"] = pd.to_datetime(df_raw["datetime"])
    df_obs = pd.DataFrame()
    df_obs["observation"] = df_raw.apply(lambda row:
        c_t.lla_to_eci(lat, long, alt, row["datetime"]), axis=1)
    df_obs["right_ascension"], df_obs["declination"] = df_raw.apply(lambda row:
        c_t.azel_to_radec(row["azimuth"], row["elevation"], lat, long,
        row["datetime"]), axis=1).str
    df_obs["right_ascension"] = np.deg2rad(df_obs["right_ascension"])
    df_obs["declination"] = np.deg2rad(df_obs["declination"])
    df_obs["rho"] = df_raw["range"]*df_obs.apply(lambda row: np.array(
            [np.cos(row["right_ascension"])*np.cos(row["declination"]),
            np.sin(row["right_ascension"])*np.cos(row["declination"]),
            np.sin(row["declination"])]), axis=1)
    df_obs["r"] = df_obs["rho"] + df_obs["observation"]
    df_obs["juliandate"] = df_raw.apply(lambda row: t_t.jd(row["datetime"]), axis=1)
    df_obs["time_diff"] = (df_obs["juliandate"][0]-df_obs["juliandate"])*24*60*60
    def send_raw_observation_to_start(r1, r2, r3, jd1, jd2, jd3, dt, mu):
        v2 = gibbs(r1, r2, r3, mu, jd1=jd1, jd2=jd2, jd3=jd3)
        y0 = np.concatenate((r2, v2))
        t_span = [0, dt]
        sol = solve_ivp(fun, t_span, y0, method='RK45',
            args=(mu,), rtol=1e-8, atol=1e-8)
        r = np.array([val[-1] for val in sol.y[0:3]])
        v = np.array([val[-1] for val in sol.y[3:6]])
        return r, v
    df_grouped = pd.DataFrame()
    obs_to_use = math.floor(len(df_obs["r"])/3)
    temp_r = df_obs["r"][0:obs_to_use*3]
    temp_jd = df_obs["juliandate"][0:obs_to_use*3]
    temp_time_diff = df_obs["time_diff"][0:obs_to_use*3+1]
    for x in range(obs_to_use):
        df_grouped[f"r{x+1}"] = temp_r[x::3].values
        df_grouped[f"jd{x+1}"] = temp_jd[x::3].values
    df_grouped["time_diff"] = temp_time_diff[1::3].values
    if len(df_grouped["time_diff"]) > obs_to_use:
        df_grouped["time_diff"] = df_grouped["time_diff"].drop(-1)
    df_grouped["r"], df_grouped["v"] = df_grouped.apply(lambda row:
        send_raw_observation_to_start(row["r1"], row["r2"], row["r3"],
        row["jd1"], row["jd2"], row["jd3"], row["time_diff"], mu), axis=1).str
    r_nominal = np.array(df_grouped["r"]).mean()
    v_nominal = np.array(df_grouped["v"]).mean()
    return r_nominal, v_nominal, df_obs

# Non-linear Weighted Least Squares
def nlwls(lat, long, alt, df_raw,
    range_err, az_err, el_err, fun, mu):
    r_nominal, v_nominal, df_obs = find_nominal(
        lat, long, alt, df_raw, fun, mu)
    # While loop setup
    diff = 1000
    tol = 1e-3
    k = 1
    x_nominal = np.concatenate((r_nominal, v_nominal))
    df_raw["datetime"] = pd.to_datetime(df_raw["datetime"])
    df_true_obs = pd.concat(
        [df_raw["range"], df_raw["azimuth"], df_raw["elevation"]], axis=1)
    t_span = []
    w = np.array([
        [(1/range_err)**2, 0, 0],
        [0, (1/az_err)**2, 0],
        [0, 0, (1/el_err)**2]
        ])
    prev_rms = 0
    n_obs = len(df_raw["range"])
    n_type_obs = len(df_true_obs.iloc[0,:])
    df_for_loop = pd.concat([df_raw["datetime"], df_obs["time_diff"],
            df_true_obs], axis=1)
    while diff > tol:
        htwh = np.zeros((len(x_nominal), len(x_nominal)))
        htwy = np.zeros(len(x_nominal))
        for _, row in df_for_loop.iterrows():
            if row["time_diff"] == 0:
                r_nom_obs = x_nominal[0:3]
            else:
                t_span = [0, -row["time_diff"]]
                y0 = x_nominal
                sol = solve_ivp(fun, t_span, y0, method='RK45',
                    args=(mu,), rtol=1e-8, atol=1e-8)
                r_nom_obs = np.array([val[-1] for val in sol.y[0:3]])
            rho, az, el = c_t.state_to_razel(
                r_nom_obs, row["datetime"], lat, long, alt)
            x_nom_obs = np.array([rho, az, el])
            y = np.array(list(row[2:])) - x_nom_obs
            h = np.zeros((len(x_nominal), len(x_nom_obs)))
            for x in range(len(x_nominal)):
                x_mod = copy.deepcopy(x_nominal)
                x_mod[x] = x_nominal[x]*(1.001)
                if row["time_diff"] == 0:
                    r_mod_obs = x_mod[0:3]
                else:
                    y0 = x_mod
                    sol = solve_ivp(fun, t_span, y0, method='RK45',
                        args=(mu,), rtol=1e-8, atol=1e-8)
                    r_mod_obs = np.array([val[-1] for val in sol.y[0:3]])
                rho, az, el = c_t.state_to_razel(
                    r_mod_obs, row["datetime"], lat, long, alt)
                h[:][x] = (np.array([rho, az, el])-x_nom_obs)/(x_nominal[x]*0.001)
            htwh += np.dot(h, np.dot(w, h.T))
            htwy +=  np.dot(h, np.dot(w, y))
        p = np.linalg.pinv(htwh)
        dx_hat = np.dot(p, htwy)
        rms = np.sqrt(np.dot(y, np.dot(w, y))/(n_type_obs*(n_obs)))
        if k > 1:
            diff = np.abs((prev_rms-rms)/prev_rms)
        prev_rms = rms
        x_nominal += dx_hat
        k += 1
        if k == 1000:
            raise Exception
    eig_w, eig_v = np.linalg.eig(p)
    confidence = [np.sqrt(p[x, x])*1000 for x in range(3)]
    print(f"Took {k} iterations.")
    print(f"Position vector: {x_nominal[0:3]} [km]")
    print(f"Velocity vector: {x_nominal[3:6]} [km/s]")
    print(f"Confidence in position vector: {confidence} [m]")
    print(f"Error ellpisoid of position vector: {eig_w*1000}")
    print(f"Covariance matrix: {p}")
    print(f"Eigenvalues: {eig_w}")
    print(f"Eigenvectors: {eig_v}")
    