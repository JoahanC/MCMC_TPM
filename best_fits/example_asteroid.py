from asteroids import *

#Create object
asteroid_ex = Asteroid("asteroid_ex")

#Call functions from asteroids.py
print("Calling read_csh_inputs function. The function has the C shell script serve as an input for the MCMC model, containing information such as H mag, H mag error, and period.")
asteroid_ex.read_csh_inputs()
print("Calling retrieve_SED_data function. The function extracts data from the fort.22 file, to set up the SED plot.")
asteroid_ex.retrieve_SED_data()
print("Calling retrieve_MCMC_results function. The function extracts data from the fort.21 file and obtains diameter, albedo, gamma, chi squared, and period values for the asteroid.")
asteroid_ex.retrieve_MCMC_results()
print("Calling retrieve_unequal_data function. The function attempts to recover data missing from fort.4, fort.3, or fort.21 files if they have missing information.")
asteroid_ex.retrieve_unequal_data()
print("Calling generate_SED_plot function. The function constructs a best fit SED plot, depicting the average flux density during the epochs.")
asteroid_ex.generate_SED_plot()
#print("Calling clear_plot_directories function. The function clears out all plots from the current directory.")
#asteroid_ex.clear_plot_directories()
print("Calling generate_histograms function. This function constructs a histogram for the diameter, albedo, gamma, and period results ")
asteroid_ex.generate_histograms()
print("Calling generate_chi_scatterplots function. This function constructs scatterplots for all characterization solutions.")
asteroid_ex.generate_chi_scatterplots()
print("Calling generate_hexbins function. This function constructs hexbins for all characterization solutions")
asteroid_ex.generate_hexbins()
print("Calling generate_chi_plots function. This function constructs scatterplots demonstrating the relationship between the characterization solutions vs the chi squared values from the MCMC model.")
asteroid_ex.generate_chi_plots()

