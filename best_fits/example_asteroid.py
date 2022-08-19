"""
This file documents an example case of using the Asteroid class to perform analysis
on a specific asteroid. Feel free to comment lines out as necessary.
"""
from asteroids import *


"""
To initialize the Asteroid class, the directory containing the asteroid to model will
be used as the constructor argument.
"""
asteroid_ex = Asteroid("asteroid_ex")

"""
The Asteroid class contains a set of call functions which can interface with the direct
output files of the TPM model. These are automatically run in the initializer of the
Asteroid class but can be run individually as shown below.
"""
print("Calling read_csh_inputs function.")
asteroid_ex.read_csh_inputs()
print("Calling retrieve_SED_data function.")
asteroid_ex.retrieve_SED_data()
print("Calling retrieve_MCMC_results function.")
asteroid_ex.retrieve_MCMC_results()

"""
The calling functions make use of the fort.21 and fort.4 files for retrieving albedo
outputs. However, in some rare cases, information can be missed and there may be
a mismatch between the two files. A method is designed to handle this specific
case: `retrieve_unequal_data` however it is currently optimized for being used for
02100 Ra-Shalom and targets specific skip cases when doing a line by line comparison
of the diameter values of both files. This is an experimental method and needs more
work.
"""
#print("Calling retrieve_unequal_data function.")
#asteroid_ex.retrieve_unequal_data()

"""
Next are the plotting capablities of the Asteroid class. A spectral density flux
distribution plot can be generated with the following method. It performs a fit
on each set of epoch flux values in however many wavelength bands were provided,
as well as a set of individual median flux points and uncertainties.
"""
asteroid_ex.generate_SED_plot()

"""
Histograms of each output parameter: (Diameter, Visual Albedo, Thermal Inertia, 
and Period if available) can be generated with corresponding uncertainty lines.
"""
asteroid_ex.generate_histograms()

"""
Fit chi squared labeled scatterplots of each output parameters can be generated
in all possible axis configurations.
"""
asteroid_ex.generate_chi_scatterplots()

"""
Density hexbin plots of each output parameters can be generated
in all possible axis configurations.
"""
asteroid_ex.generate_hexbins()

"""
Finally, chi vs output parameter plots can generated for all available output parameters.
"""
asteroid_ex.generate_chi_plots()

"""
There is also a directory clearing method which will eliminate all current plots from 
a folder.
"""
asteroid_ex.clear_plot_directories()