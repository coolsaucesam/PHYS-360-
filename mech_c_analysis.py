import pandas as pd
import numpy as np
from scipy.optimize import curve_fit  
import matplotlib.pyplot as plt 

# Read in data
data = pd.read_csv("data_2.csv", header=0)

rev = data["speed"]
deflection = data["deflection"] * 10**(-3) # factor of 10e-3 to convert mm to m

lit_c = 299792458 # [m/s] speed of light in literature.

# Constants and Uncertainties
A = 0.267  # distance between L2 and L1 in mm
deltaA = 0.5 * 10**(-3)  # +- 0.5mm
B = 0.483  # distance between L2 and Mr in mm
deltaB = 0.5 * 10**(-3)  # +- 0.5 mm
D = 9.94  # distance between Mr and Mf
deltaD = 0.1  # uncertainty in D
deltaSPrime = 0.005 * 10**(-3)  # uncertainty in s'
deltaOmega = 1  # uncertainty in omega

# constant combining all constants in equation for c
k = (8 * np.pi * A * D**2) / (B + D)

# Calculate the uncertainty in k
uncertainty_in_k = np.sqrt(
    ( (8*np.pi * D**2) / (B + D) * deltaA)**2 +
    ( (8*np.pi * A * D**2) / (B+D**2) * deltaB)**2 +
    ( (8*np.pi * A * D * (2*B + D))/((B+D)**2) * deltaD)**2)

# Equation from lab
def return_c(rev, deflection):
    return (8 * np.pi * A * D**2 * rev) / ((D + B) * deflection)

# Linear model
def model(x, H, J):
    return H * x + J

# Calculate 'c' values and linear fit data
# Clist = []
# for h in range(len(deflection)):
#     Clist.append(return_c(rev[h], deflection[h]))


# Now, implement curve fitting with uncertainty
initial_guesses = np.array([0, 0])  # Speed of light guess and offset

# Uncertainty vector
# uncertainties = [calculateUncertainty(rev[i], deflection[i]) for i in range(len(rev))]

# Perform curve fitting with uncertainties (use 'sigma' to provide weights for the fit)
params, cov = curve_fit(model, rev, deflection, p0=initial_guesses)

# calculating the uncertainty in c using k/slope
print(np.sqrt(np.diag(cov)))
uncertainty_in_c = np.sqrt(
    (k/((params[0])**2) *(np.sqrt(np.diag(cov)))[0])**2 +
    (1/(params[0]*uncertainty_in_k)**2)
    )

    
print(uncertainty_in_c)
exp_c = k/params[0]
# Print results and analyze
print(f"Fitted Parameters: H = {params[0]:.3e}, J = {params[1]:.3e}")
print(f"Speed of light from fit: {k / params[0]:.3e} plus/minus {uncertainty_in_c:.2e} m/s")
print('percent error in c', round(uncertainty_in_c/exp_c, 5))

# Plot the results with uncertainties
X = np.linspace(min(rev), max(rev), 100)
fig, ax = plt.subplots()
# Plot the linear fit
ax.scatter(rev, deflection, color='blue')
ax.plot(X, model(X, params[0], params[1]), label="Linear Fit", color='orange')
# Labels and title
ax.set_xlabel("Rotation Speed (rev/s)")
ax.set_ylabel("Deflection (m)")
ax.set_title("Linear Fit to Speed of Light Experiment")
ax.grid()
# Show plot
plt.savefig('mech_c_fit.png', dpi=200)
plt.show()
