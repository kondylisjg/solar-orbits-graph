from math import *

au = 149600 # in million meters


''' Polar equation of an elliptical orbit with star at the origin
    and the start being the smaller of the two foci
    (i.e. start is the "left" focus):
    
    r = p / (1 - e * cos(theta))
    p is the semi latus-rectum
    e is the eccentricity
'''

class Orbit:
    def __init__(self, period, eccentricity, alpha):
        self.period = period
        self.eccen = eccentricity
        self.alpha = alpha

    # Latus Rectum as function of semi-major axis and eccentricity
    def lr(self):
        return self.alpha * (1.0 - self.eccen**2)

    def beta_from_e(self):
        return self.alpha * sqrt(1.0 - self.eccen**2)

    def focus(self):
        b = self.beta_from_e()
        # Focus if ellipse was centered at origin
        focus = sqrt(self.alpha**2 - b**2)
        return focus

    def eccentricity(self, beta):
        focus = sqrt(self.alpha**2 - beta**2)
        return 1.0 * focus / self.alpha

    '''Angular velocity equation based on Kepler's second Law
       which says that orbit covers equal area in equal time:
       Period * r^2 * omega = 2*pi*alpha*beta
       (Substitute r from polar form of ellipse)
    '''
    def omega(self, theta):
        p = self.lr()
        b = self.beta_from_e()
        v = 2*pi * self.alpha * b * (1 - self.eccen * cos(theta))**2 / (self.period * p**2)
        return v

    ''' Use angular velocity (above) and a standard Runge-Kuta method for solving ODEs
        to find theta as a function of time and thus get polar co-ordinates of orbit 
    '''
    def solve_orbit(self, num_points):
        time = [0]
        theta = [0]
        dt = 1.0*self.period / num_points
        for i in range(1, num_points):
            y = theta[i-1]
            k1 = dt * self.omega(y)
            k2 = dt * self.omega(y + k1/2.0)
            k3 = dt * self.omega(y + k2/2.0)
            k4 = dt * self.omega(y + k3)
            y += 1.0 / 6.0 * (k1 + 2*k2 + 2*k3 + k4)
            print("%f %f" % (i * dt, y))
            theta.append(y)
            time.append(i * dt)

        # Solve for radius
        p = self.lr()
        radius = []
        for i in range(num_points):
            r = p / (1.0 - self.eccen * cos(theta[i]))
            radius.append(r)

        return time, radius, theta

    def rotate_2d(self, phi, x, y):
        phi = phi * pi / 180.0
        x_rot = []; y_rot = []
        for i in range(len(x)):
            x_rot.append(x[i] * cos(phi) - y[i] * sin(phi))
            y_rot.append(x[i] * sin(phi) + y[i] * cos(phi))
        return x_rot, y_rot

    # Incline orbit wrt x-y (earth's) plane
    def tilt_3d(self, inclination, radius, theta):
        x = []; y = []
        for i in range(len(radius)):
            x.append(radius[i] * cos(theta[i]))
            y.append(radius[i] * sin(theta[i]))

        phi = inclination * pi / 180.0
        x_tilt = []; y_tilt = []; z_tilt = []
        for i in range(len(x)):
            x_tilt.append(x[i] * cos(phi))
            y_tilt.append(y[i])
            z_tilt.append(-sin(phi) * x[i])
        return x_tilt, y_tilt, z_tilt



# Abnoba comet orbit parameters
eccentricity = 0.17883
alpha = 2.78907 * au
inclination = 14.43869223703236
rotation = 229.2073001218772
period = 1701.3  # in days

Abnoba = Orbit(period, eccentricity, alpha)
time, radius, theta = Abnoba.solve_orbit(int(period))
Ax, Ay, Az = Abnoba.tilt_3d(inclination, radius, theta)


# Earths' orbit parameters:
eccentricity = 0.0167
alpha = 1.0 * au
period = 365.256

Earth = Orbit(period, eccentricity, alpha)
e_time, e_radius, e_theta = Earth.solve_orbit(int(period))
Ex, Ey, Ez = Earth.tilt_3d(0, e_radius, e_theta)
Ex, Ey = Earth.rotate_2d(-rotation, Ex, Ey)


# Store in Cartesian co-ordinates
with open('earth.txt', 'w') as f:
    for i in range(len(e_time)):
        line = str(e_time[i])+" "+str(Ex[i])+" "+str(Ey[i])+" "+str(Ez[i])+'\n'
        f.write(line)

with open('abnoba2d.txt', 'w') as f:
    for i in range(len(time)):
        line = str(time[i])+" "+str(Ax[i])+" "+str(Ay[i])+" "+str(Az[i])+'\n'
        f.write(line)

