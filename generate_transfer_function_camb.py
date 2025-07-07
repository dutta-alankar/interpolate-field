import camb
from camb import model, initialpower
import numpy as np

# Set the redshift
z_target = 200.0

# Set up the cosmology parameters (Planck 2018 Lambda-CDM)
h = 0.673
Omega_m			= 0.314
Omega_L			= 0.686
w_0			    = -1.0
w_a			    = 0.0
Omega_b			= 0.049
sigma_8			= 0.812
n_s			    = 0.965

pars = camb.CAMBparams()
pars.set_cosmology(H0=h*100, ombh2=Omega_b*h**2, omch2=(Omega_m-Omega_b)*h**2, mnu=0.06, tau=0.0561)
pars.InitPower.set_params(As=2.105e-9, ns=n_s)
pars.set_matter_power(redshifts=[z_target], kmax=3365.0)  # kmax must be large enough

# Calculate transfer functions
results = camb.get_results(pars)
transfers = results.get_matter_transfer_data()

print(transfers.transfer_data.shape)
# print(results.get_sigma8_0())
# Write the data in MUSIC format
np.savetxt(f"transfer_camb_z{int(z_target)}.txt", transfers.transfer_data[:,:,0].T)
print(f"Transfer functions written for CDM and baryons at z = int(z_target).")
