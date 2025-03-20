import pandas as pd
import numpy as np
from scipy.optimize import curve_fit  
import matplotlib.pyplot as plt 

mk1 = 'o'
ms1 = 50
T_room = 21 + 273 # room temp in C adjusted for Kelvin
T_ice = 4 + 273
T_boil = 70 + 273
T_iso = -77 + 273
T_nitro = -196 + 273

deltaT_room = 1
deltaT_ice = 3 # degrees C
deltaT_boil = 2 # degrees C


# Read in data
data = pd.read_csv("ek_data.csv", header=0)
# room temp data
I_room = data['I']*10**(-6) # current is recorded in microamps
V_room = data['V']*10**(-4) # voltage is recorded in millivolts
# ice bath data
I_ice = data['I_ice']*10**(-6)
V_ice = data['V_ice']*10**(-4)
# boiling water data
I_boil = data['I_boil']*10**(-6)
V_boil = data['V_boil']*10**(-4)
# dry ice isopropanol
I_iso = data['I_iso']*10**(-6) 
V_iso = data['V_iso']*10**(-4)
# liquid nitrogen 
I_nitro = data['I_nitro']*10**(-6)
V_nitro = data['V_nitro']*10**(-4)

# Linear model
def model(x, A, B):
    return A * x + B

def advanced_model(y, A, B, C):
    return A + B*np.log(y) + C*y**(-1)

boltz = 1.3806*10**(-23)
initial_guesses = np.array([0, 0])
advanced_p0 = np.array([0, 100, -11000])

#calculating lines of best fit for each data set

params_room, cov_room = curve_fit(model, V_room, np.log(I_room), p0=initial_guesses)
params_ice, cov_ice = curve_fit(model, V_ice, np.log(I_ice), p0=initial_guesses)
params_boil, cov_boil = curve_fit(model, V_boil, np.log(I_boil), p0=initial_guesses)
params_iso, cov_iso = curve_fit(model, V_iso, np.log(I_iso), p0=initial_guesses)
params_nitro, cov_nitro = curve_fit(model, V_nitro, np.log(I_nitro), p0=initial_guesses)

x = np.linspace(.6, 1.1, 100)
y = np.linspace(.0029, .013, 100)
z = np.linspace(50, 350, 1000)


#calculate e/k:
# slope is e/kT
# so e/k is the slope*T

# Basic average of data
ek_room = params_room[0]*T_room
ek_ice = params_ice[0]*T_ice
ek_boil = params_boil[0]*T_boil
ek_nitro = params_nitro[0]*T_nitro
print('nitro', ek_nitro)

ek_avg = (ek_room + ek_ice + ek_boil)/3
print(round(ek_avg, 1))

#I_0 to get band gap
intercepts = [params_room[1], params_ice[1], params_boil[1], params_iso[1], params_nitro[1]]
temps = [T_room, T_ice, T_boil, T_iso, T_nitro]
inv_temps = [T**(-1) for T in temps]

# params_I0, cov_I0 = curve_fit(model, inv_temps, intercepts, p0=initial_guesses)
params_I0, cov_I0 = curve_fit(advanced_model, temps, intercepts, p0=advanced_p0)

#print(np.sqrt(np.diag(cov_I0)))


#UNCERTAINIES
# grab uncertainties from cov matrix
def CovUnc(cov):
    return np.sqrt(np.diag(cov))

#function to calculate the percent uncertainty in each measurement from slope fit
def SingleUncertaintyEK(fit_slope, delta_slope, temp, delta_temp):
    return np.sqrt((delta_slope/fit_slope)**2+(delta_temp/temp)**2)

# call cov function
unc_room = CovUnc(cov_room)
unc_ice = CovUnc(cov_ice)
unc_boil = CovUnc(cov_boil)

# calculate the uncertainties 
roomEK_unc = SingleUncertaintyEK(params_room[0], unc_room[0], T_room, deltaT_room)
iceEK_unc = SingleUncertaintyEK(params_ice[0], unc_ice[0], T_ice, deltaT_ice)
boilEK_unc = SingleUncertaintyEK(params_boil[0], unc_boil[0], T_boil, deltaT_boil)
print(roomEK_unc)
print(iceEK_unc)
print(boilEK_unc)

# implement average
totalEK_unc = (roomEK_unc + iceEK_unc + boilEK_unc)/3
print(totalEK_unc)




#PLOTTING

fig, ax = plt.subplots()
plt.title(r'LN $I_0$ vs. $T$')
ax.scatter(temps, intercepts, color='black', s=50)
ax.plot(z, advanced_model(z, params_I0[0], params_I0[1], params_I0[2]), color='blue')
plt.ylabel(r'LN $I_0$ [LN A]')
plt.xlabel(r'$T$ [$K$]')
plt.savefig('kircup_model.png', dpi=200)
print('band gap from slope:', (-params_I0[2])*boltz*6.241509074461e18)