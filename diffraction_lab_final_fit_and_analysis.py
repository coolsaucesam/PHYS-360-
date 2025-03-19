import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from scipy.optimize import curve_fit

#import and initialize data
data = pd.read_csv("single_slit_data.csv", header=0)

# multiplied to get in meters
pos_80 = np.array(data["80mim_pos"])*10**(-3)
volts_80 = np.array(data["80_volts"])

pos_40 = np.array(data["40_pos"])*10**(-3)
volts_40 = np.array(data["40_volts"])

X = np.linspace(0, 0.046, 1000)
Y = np.linspace(0, 0.092, 1000)

# initialize constants
wavelength = 633*10**(-9)
k = (2*np.pi)/wavelength #wave number
a80 = 80*10**(-6) # nominal slit width
a40 = 40*10**(-6) # nominal width
L = 1.04  # measured @ 104cm +- 1mm - separation between slit and diffraction pattern

# lets create a combined constant q_0 that is 2L/ka 
q0_80 = (2*L)/(k * a80) # initial guess value
q0_40 = (2*L)/(k * a40)


I_0 = 3*10**(-3) #maximum intensity 
I0_40 = 3*10**(-2)


# potentially lump L,k, and a into one param
# fringe model for a single slit
# params are x_0: center fringe location
# q: defined above 
# I: maximum intensity. 
def fringe_model(x, x_0, q, I):
    return I * (q/(x-x_0))**2 * (np.sin((x-x_0)/q))**2

# lets calculate the uncertainty in a:
# deltax is the uncertainty in L, and delta y is the uncertainty in q
def calculateUncertainties(deltax, deltay, q):
    return np.sqrt((2/(k*q)*deltax)**2+
                    ((2*L)/(k*q**2)*deltay)**2)

# Now to determine appropriate initial values to create our best fit model
initialGuesses80 = np.array([.023, q0_80, I_0])
initialGuesses40 = np.array([.046, q0_40, I0_40])
# testing out guesses with
# calculatedGuesses80 = fringe_model(pos_80, 23, q_0, I_0)
# calculatedGuess40 = fringe_model(Y, 46, q0_40, I0_40)
# performing the fit
params80Micrometers, cov80Micrometers = curve_fit(fringe_model, pos_80, volts_80, p0=initialGuesses80)
params40Micrometers, cov40Micrometers = curve_fit(fringe_model, pos_40, volts_40, p0=initialGuesses40)

print('---------------------Uncertainty for 80 micrometer slit-------------------\n')
uncertainties80 = np.sqrt(np.diag(cov80Micrometers))
print("x_0:", params80Micrometers[0])
print("q:", params80Micrometers[1])
print()
slit = (2*L)/(abs(params80Micrometers[1])*k)*10**6
print("a:", (2*L)/(abs(params80Micrometers[1])*k)*10**6)
#print("The uncertainties in fit are:", uncertainties80)
print("absolute: +-", calculateUncertainties(.001, uncertainties80[1], params80Micrometers[1])*10**6, "μm -- % uncertainty:", calculateUncertainties(.001,uncertainties80[1], params80Micrometers[1])/((2*L)/(abs(params80Micrometers[1])*k)))
print((80-slit)/80)

print('\n----------------Uncertainty for 40 micrometer slit-----------------\n')
uncertainties40 = np.sqrt(np.diag(cov40Micrometers))
print("x_0:", params40Micrometers[0])
print("q:", params40Micrometers[1])
slit40 = (2*L)/(abs(params40Micrometers[1])*k)*10**6
print("a:", (2*L)/(abs(params40Micrometers[1])*k)*10**6)
#print("The uncertainties in fit are:", uncertainties40)
print("absolute:", calculateUncertainties(.001, uncertainties40[1], params40Micrometers[1])*10**6, "μm -- % uncertainty:", calculateUncertainties(.001, 7.35714755e-05, params40Micrometers[1])/((2*L)/(abs(params40Micrometers[1])*k)))
print('percent difference:', (40-slit40)/40)

#calculating residuals:
residuals80 = []
residuals40 = []
for x in range(len(volts_80)):
    residuals80.append(volts_80[x] - fringe_model(x*10**(-3), params80Micrometers[0], params80Micrometers[1], params80Micrometers[2]))

for y in range(len(volts_40)):
    residuals40.append(volts_40[y] - fringe_model(y*2*10**(-3), params40Micrometers[0], params40Micrometers[1], params40Micrometers[2]))

fig, (ax1, ax2) = plt.subplots(1,2, layout='constrained', figsize=(12,5))


#plotting raw dat9 for 80micrometer slit
ax1.set_title('(a)', loc='left')
ax1.plot(pos_80, volts_80, marker='o', ls='None', color='black', label="data", markersize=5)
ax1.plot(pos_80, fringe_model(pos_80, params80Micrometers[0], params80Micrometers[1], params80Micrometers[2]), color='black', label="best fit curve", linewidth=1.5, linestyle='dashed')
# plt.plot(pos_80, residuals80, linestyle='dashed', label='residuals', color='black')
ax1.set_xticks([ 0.003, .013, .023655, .034, .044], [ '-20', '-10' ,'0', '10','20'])
# plt.axvline(params80Micrometers[0], color='black', linestyle='dotted')
ax1.set_ylabel(" Volts [V]")
ax1.set_xlabel(" Position [mm]")
#plt.legend()
ax1.grid()

#plotting data from 40μm
ax2.set_title('(b)', loc='left')
ax2.plot(pos_40, volts_40, 'o', color='blue', label='data', markersize=5)
ax2.plot(pos_40, fringe_model(pos_40, params40Micrometers[0], params40Micrometers[1], params40Micrometers[2]), color='b', linestyle='dashed', linewidth=1.5,  label='best fit curve')
# plt.plot(pos_40, residuals40, linestyle='dashed', label='residuals', color='blue')
ax2.set_xticks([0.04685-.04, 0.04685-.02, .04685,0.04685+.02, 0.04685+.04], ['-40', '-20' ,'0', '20','40'])
ax2.set_ylabel("Volts [V]")
ax2.set_xlabel("Position [mm]")
#plt.legend()
ax2.grid()
plt.savefig('data_plot.png', dpi=250)
  
# plot of the residuals  
# plt.figure()
# plt.title('residuals')
# plt.plot(pos_80, residuals80, label='80 micrometer', color='black')
# plt.plot(pos_40, residuals40, label='40 micrometer', color='blue')
# plt.savefig("residuals.png", dpi=1000)
# plt.figure()
# plt.plot(Y, fringe_model(Y, 46, q0_40, I0_40))