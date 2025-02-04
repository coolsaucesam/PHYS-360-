import numpy as np
import matplotlib.pyplot as plt

#import function for non-linear fitting
from scipy.optimize import curve_fit


X = np.linspace(-100, 100, 100000)

#initial values
wavelength = 594*10**(-9)
a = 80*10**(-6) #slit width
L = 2 #separation between slit and diffraction pattern
k = (2*np.pi)/wavelength #wave number

c = 299.792*10**8 #[m/s] speed of light
epsilon_0 = 8.8542*10**(-12) #[F*m^-1]

I_0 = (c*epsilon_0*a**2)/(2*(k*L)**2)

print(I_0)

def basic_fringe_model(x):
    return x**(-2)*(np.sin(x/2))**2

def fringe_model(x, I, L, k, a):
    return I*((2*L)/(k*a*x))**2*(np.sin((k*a*x)/(2*L)))**2


#plotting data
plt.plot(X, basic_fringe_model(X))
plt.plot(X, fringe_model(X, 1, 1, 1, 1))
plt.xlabel('Displacement (m)')
plt.ylabel('Spring force (N)')
plt.show()