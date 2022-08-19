"""
This file declares the Asteroid class which serves as a representation of the
thermophysical modeling outputs and inputs for each object in this repository.
"""
import os
import itertools
import numpy as np
from helpers import *


class Asteroid(object):
    """
    The Asteroid object holds all of the spherical and triaxial
    results following an MCMC run.
    """
    def __init__(self, directory):
        """
        Constructor for the Asteroid class.

        Parameters
        ----------
        directory : str
            The folder name of the asteroid.
        """
        self.directory = directory
        self.is_triaxial = False
        self.packed_name = directory

        if "triaxial_" in directory:
            self.packed_name = directory.replace("triaxial_", '')
            self.is_triaxial = True

        self.csh_file = f"../{directory}/run_{self.packed_name}.csh"
        with open(f"../{directory}/fort.21", 'r') as chain_file:
            chain_lines = chain_file.readlines()
        with open(f"../{directory}/fort.4", 'r') as albedo_file:
            albedo_lines = albedo_file.readlines()
        self.full_lines = False
        if len(chain_lines) == len(albedo_lines) + 1:
            self.full_lines = True

        best_fit_outputs = open(f"{self.directory}.txt", 'r')
        self.b_a_ratio = 0
        self.b_a_pos_sig = 0
        self.b_a_neg_sig = 0
        self.c_b_ratio = 0
        self.c_b_pos_sig = 0
        self.c_b_neg_sig = 0
        self.eigenvalues = {}
        eig_idx = 1
        for line in best_fit_outputs:
            if "dia=" in line:
                diameter_vals = determine_mean_median_vals(line, "multi")
                self.diameter_mean = float(diameter_vals[0])
                self.diameter_mean_16th_percentile = float(diameter_vals[1])
                self.diameter_mean_84th_percentile = float(diameter_vals[1])
                self.diameter_median = float(diameter_vals[2])
                self.diameter_median_16th_percentile = float(diameter_vals[3])
                self.diameter_median_84th_percentile = float(diameter_vals[4])

            if "p_V = " in line:
                p_V_vals = determine_mean_median_vals(line, "single")
                self.albedo_mean = float(p_V_vals[0])
                self.albedo_mean_16th_percentile = float(p_V_vals[1])
                self.albedo_mean_84th_percentile = float(p_V_vals[1])
                self.albedo_median = float(p_V_vals[2])
                self.albedo_median_16th_percentile = float(p_V_vals[3])
                self.albedo_median_84th_percentile = float(p_V_vals[3])

            if "theta1=" in line:
                theta_vals = determine_mean_median_vals(line, "multi")
                self.theta_mean = float(theta_vals[0])
                self.theta_mean_84th_percentile = float(theta_vals[1])
                self.theta_mean_16th_percentile = float(theta_vals[1])
                self.theta_median = float(theta_vals[2])
                self.theta_median_16th_percentile = float(theta_vals[3])
                self.theta_median_84th_percentile = float(theta_vals[4])

            if "Period [h]" in line:
                period_vals = determine_period(line)
                self.period_median = float(period_vals[0])
                self.period_16th_percentile = float(period_vals[1])
                self.period_84th_percentile = float(period_vals[2])

            if "sqrt(kappa*rho*C)" in line:
                gamma_vals = determine_gamma_vals(line)
                self.gamma_median = float(gamma_vals[0])
                self.gamma_16th_percentile = float(gamma_vals[1])
                self.gamma_84th_percentile = float(gamma_vals[2])

            if "crater fraction" in line:
                crater_vals = determine_crater_fraction(line)
                self.crater_fraction_median = float(crater_vals[0])
                self.crater_fraction_16th_percentile = float(crater_vals[1])
                self.crater_fraction_84th_percentile = float(crater_vals[2])

            if "p_IR/p_V=" in line:
                ratio_vals = determine_p_V_ratio(line)
                self.ir_fraction_median = float(ratio_vals[0])
                self.ir_fraction_16th_percentile = float(ratio_vals[1])
                self.ir_fraction_84th_percentile = float(ratio_vals[2])

            if "b/a" in line:
                face_ratios = determine_face_ratios(line)
                self.b_a_ratio = face_ratios[0]
                self.b_a_pos_sig = face_ratios[1]
                self.b_a_neg_sig = face_ratios[2]

            if "c/b" in line:
                face_ratios = determine_face_ratios(line)
                self.c_b_ratio = face_ratios[0]
                self.c_b_pos_sig = face_ratios[1]
                self.c_b_neg_sig = face_ratios[2]

            if "pole peak at = " in line:
                pole_peak_vals = line.split()
                self.pole_peak_ra = float(pole_peak_vals[4])
                self.pole_peak_dec = float(pole_peak_vals[5])

            if "mean pole at = " in line:
                mean_pole_vals = line.split()
                self.pole_mean_ra = float(mean_pole_vals[4])
                self.pole_mean_dec = float(mean_pole_vals[5])
                self.pole_bar = float(mean_pole_vals[7])

            if "Moment eigenvalue=" in line:
                eig_vals = line.split()
                self.eigenvalues[f"eigenvalue_{eig_idx}"] = [float(eig_vals[2]), float(eig_vals[5]), float(eig_vals[6])]
                eig_idx += 1

        best_fit_outputs.close()
        SED_plotters = self.retrieve_SED_data()
        self.esed = SED_plotters[0]
        self.wavelengths = SED_plotters[1]
        self.datelabels = SED_plotters[2]
        self.epoch_count = len(self.datelabels)
        self.data_x = SED_plotters[3]
        self.data_y = SED_plotters[4]
        self.data_y_error = SED_plotters[5]

        MCMC_results = self.retrieve_MCMC_results()
        self.diameters = MCMC_results[0]
        self.albedos = MCMC_results[1]
        self.gammas = MCMC_results[2]
        self.chis = MCMC_results[3]
        self.periods = MCMC_results[4] 
        self.fixed = False
        if self.periods[0] == self.periods[540] and self.periods[9] == self.periods[132]:
            self.fixed = True

        csh_inputs = self.read_csh_inputs()
        self.h_magnitude = csh_inputs["h_mag"]
        self.h_error = csh_inputs["h_error"]
        self.period = csh_inputs["period"]
        self.epochs = {}
        for key in csh_inputs:
            if "epoch" in key:
                self.epochs[key] = csh_inputs[key]

        with open(f"../{directory}/fort.21", 'r') as chain_file:
            chain_lines = chain_file.readlines()
        with open(f"../{directory}/fort.4", 'r') as albedo_file:
            albedo_lines = albedo_file.readlines()

        if "general_plots" not in os.listdir(f"../{self.directory}/"):
            os.mkdir(f"../{self.directory}/general_plots")
        if "diameter_albedo_segments" not in os.listdir(f"../{self.directory}/"):
            os.mkdir(f"../{self.directory}/diameter_albedo_segments")
        if "diameter_period_segments" not in os.listdir(f"../{self.directory}/"):
            os.mkdir(f"../{self.directory}/diameter_period_segments")
        if "diameter_thermal_inertia_segments" not in os.listdir(f"../{self.directory}/"):
            os.mkdir(f"../{self.directory}/diameter_thermal_inertia_segments")


    def read_csh_inputs(self):
        """
        Reads in all the MCMC input information from the object's .csh file
        """
        csh_inputs = {}
        with open(self.csh_file, 'r') as input_file:
            while True:
                header_info = input_file.readline().split(',')
                if len(header_info) == 4:
                    break
            csh_inputs["h_mag"] = float(header_info[0])
            csh_inputs["h_error"] = float(header_info[1])
            csh_inputs["period"] = float(header_info[2])
            
            if self.is_triaxial:
                input_file.readline()
            current_line = input_file.readline()
            index = 1
            while True:
                if "LAST" in current_line:
                    break
                csh_inputs[f"epoch_{index}"] = [index]
                for datum in current_line.split(','):
                    datum.replace('+', '').replace('\n', '')
                    csh_inputs[f"epoch_{index}"].append(float(datum)) 
                current_line = input_file.readline()
                index += 1
            return csh_inputs


    def retrieve_SED_data(self):
        """
        Reads in all the bestfit SED information from the object's fort.22 file
        """
        SED_file = open(f"../{self.directory}/SED_data.out")
        datum_minor = SED_file.readline().rstrip().split()
        if self.is_triaxial:
            SED_file.readline()
        datum_major = SED_file.readline().rstrip().split()
        
        best_fit = {}
        best_fit["pole_ra"] = np.degrees(float(datum_minor[1]))
        best_fit["pole_dec"] = np.degrees(float(datum_minor[2]))
        best_fit["pv_median"] = np.e ** (float(datum_minor[3]))
        best_fit["period"] = np.e ** (float(datum_minor[4]))
        best_fit["gamma"] = np.e ** (float(datum_minor[5]))
        best_fit["pjdfc"] = (float(datum_minor[6])) #period * thermal inertia Gamma * diameter * crater frac * color
        x = float(datum_minor[7])
        best_fit["crater_frac"] = np.e ** (x) / (1 + np.e ** x)
        best_fit["ir_albedo_ratio"] = np.e ** (float(datum_minor[8]))
        best_fit["chisq"] = float(datum_minor[10])
        best_fit["penalties"] = float(datum_minor[9]) - best_fit['chisq']
        best_fit["pv"] = float(datum_major[1])
        best_fit["diameter"] = float(datum_major[2])
        best_fit["theta"] = float(datum_major[3])

        SED_file.readline()
        epoch_condition = {}
        epoch_condition["delta"] = []
        epoch_condition["r_helio"] = []
        epoch_condition["phase"] = []
        epoch_condition["sub_sun_lat"] = []
        epoch_condition["sub_earth_lat"] = []

        epoch = SED_file.readline()
        while epoch[0:4] != "/wvl":
            epoch_data = epoch.split()
            epoch_condition["delta"].append(float(epoch_data[1]))
            epoch_condition["r_helio"].append(float(epoch_data[2]))
            epoch_condition["phase"].append(float(epoch_data[3]))
            epoch_condition["sub_sun_lat"].append(float(epoch_data[4]))
            epoch_condition["sub_earth_lat"].append(float(epoch_data[5]))
            epoch = SED_file.readline()

        # Reads in all wavelength data for corresponding epochs
        wavelengths = []
        wave_line = SED_file.readline()
        while "def" not in wave_line:
            wavelength_data = wave_line.strip().split()
            for wavelength in wavelength_data:
                wavelengths.append(float(wavelength))
            wave_line = SED_file.readline()
        
        color = []
        esed = []
        data_x = []
        data_y = []
        data_y_error = []
        first_line = True

        while True:
            line = SED_file.readline()
            if "SRGB" not in line and first_line:
                color.append("#444444")
                data_x.append([])
                data_y.append([])
                data_y_error.append([])
            # Ends cycle
            if line == "":
                break
            # Gather fluxes for entire epoch
            if "/ft" in line:
                esed.append([])
                line = SED_file.readline()
                while "def doit" not in line:
                    flux = float(line.strip())
                    esed[-1].append(flux)
                    line = SED_file.readline()
            # Acquire x and y data for epoch
            elif "QQ" in line:
                x_datum, error_datum, y_datum = line.rstrip().split()[0:3]
                data_x[-1].append(float(x_datum))
                data_y[-1].append(float(y_datum))
                data_y_error[-1].append(float(error_datum))
            elif "SRGB" in line:
                rr, gg, bb = line.rstrip().split()[0:3]
                red = hex(int(float(rr) * 255))[2:]
                green = hex(int(float(gg) * 255))[2:]
                blue = hex(int(float(bb) * 255))[2:]
                if len(red) == 1:
                    red = '0' + red
                if len(green) == 1:
                    green = '0' + green
                if len(blue) == 1:
                    blue ='0' + blue
                
                color.append('#%2s%2s%2s'%(red, green, blue))
                data_x.append([])
                data_y.append([])
                data_y_error.append([])
            elif "setgray" in line:
                color.append("#444444")
                data_x.append([])
                data_y.append([])
                data_y_error.append([])

            first_line = False

        datelabels = []
        cshfile = os.popen(f"ls ../{self.directory}/run*.csh").readline().rstrip()
        cshlines = open(cshfile)
        split_length = 12
        if self.is_triaxial:
            split_length = 16
        for line in cshlines.readlines():
            epoch_data = line.rstrip().split(',')
            if len(epoch_data) != split_length:
                continue
            mjd = float(epoch_data[0])
            year, month, day = julian_days_utc_converter(2400000.5 + mjd)
            date = "%4i-%02i-%02i"%(year, month, day)
            datelabels.append(date)

        SED_plotters = [esed, wavelengths, datelabels, data_x, data_y, data_y_error]
        return SED_plotters


    def retrieve_MCMC_results(self):
        """
        Generates the diameter and albedo values for an object based on the fort.21 file.
        """   
        output_file = open(f"../{self.directory}/PJDFC.out")
        output_file.readline()
        diameters = []
        chis = []
        gammas = []
        periods = []
        albedos = []

        with open(f"../{self.directory}/fort.4", 'r') as albedo_file:
            for line in albedo_file.readlines():
                albedos.append(float(line.split()[1]))

        if self.is_triaxial:
            for line in output_file.readlines():
                period = np.e ** float(line.strip()[25:34].strip())
                gamma = np.e ** float(line.strip()[34:43].strip())
                diameter = np.e ** float(line.strip()[43:52].strip())
                chi = float(line.strip()[88:105].strip())
                diameters.append(diameter)
                chis.append(chi)
                gammas.append(gamma)
                periods.append(period)
        
        if not self.is_triaxial:
            for line in output_file.readlines():            
                period = np.e ** float(line.strip()[25:34].strip())
                gamma = np.e ** float(line.strip()[34:43].strip())
                diameter = np.e ** float(line.strip()[43:52].strip())
                chi = float(line.strip()[70:87].strip())
                diameters.append(diameter)
                chis.append(chi)
                gammas.append(gamma)
                periods.append(period)

        if not self.full_lines:
            physical_chars = self.retrieve_unequal_data()
            return physical_chars

        physical_chars = [diameters, albedos, gammas, chis, periods]
        return physical_chars


    def retrieve_unequal_data(self):
        """
        In the case that the fort.4, fort.3, or fort.21 files have missing information,
        this code will run to recover as much common data between the files as possible.
        NOTE: this is optimized for 02100 Ra-Shalom currently and needs more work to
        be generalized to work for all MCMC results with this issue.
        """
        output_file = open(f"../{self.directory}/fort.21", 'r')
        albedo_file = open(f"../{self.directory}/fort.4", 'r')
        period_file = open(f"../{self.directory}/fort.3", 'r')
        output_file.readline()
        diameters = []
        diameters_true = []
        chis = []
        gammas = []
        periods = []
        albedos = []
        for line in albedo_file.readlines():
            diameters.append(float(line.strip().split()[0]))
            albedos.append(float(line.strip().split()[1]))
        
        count = 0
        exempt = [6608, 25932, 41328]
        for line in output_file.readlines():
            if round(np.e ** float(line.strip()[44:52].strip()), len(str(diameters[count])) - 2) != diameters[count] and count not in exempt:
                continue
            diameter = np.e ** float(line.strip()[44:52].strip())
            period = np.e ** float(line.strip()[25:34].strip())
            gamma = np.e ** float(line.strip()[34:43].strip())
            chi = float(line.strip()[90:106].strip())
            chis.append(chi)
            gammas.append(gamma)
            diameters_true.append(diameter)
            periods.append(period)
            count += 1
        
        output_file.close()
        period_file.close()
        albedo_file.close()
        physical_chars = [diameters_true, albedos, gammas, chis, periods]
        return physical_chars

    
    def generate_SED_plot(self):
        """
        Generates the bestfit SED plot for this asteroid.
        """
        colors = ["#3498db", "#229954", "#c0392b", "#8e44ad", "#f1c40f", "#ec7063", "#34495e", "#6e2c00"]

        fig, ax = plt.subplots()
        y_low = 1e99
        y_high = 0

        marker = itertools.cycle(('.', 'v', '^', 's', 'D')) 
        lines = itertools.cycle(('-', ':', '--', '-.'))
        for i in range(len(self.esed)):
            # Self corrects x and y limits
            if min(self.esed[i]) < y_low:
                y_low = min(self.esed[i])
            if max(self.esed[i]) > y_high:
                y_high = max(self.esed[i])

            # Plot flux for each wavelength
            ax.plot(self.wavelengths, self.esed[i], label=self.datelabels[i], color=colors[i], lw=0.75, ls=next(lines))    

            # Adjust y values and set up error ranges
            lower_errors = []
            higher_errors = []
            data_y_adjusted=[]
            for j in range(len(self.data_y[i])):
                lower_errors.append(self.data_y[i][j] * (1 - 1 / (10 ** (self.data_y_error[i][j] / 2.5))))
                higher_errors.append(self.data_y[i][j] * (10 ** (self.data_y_error[i][j] / 2.5) - 1))

                if self.data_y_error[i][j] > 0:    
                    data_y_adjusted.append(self.data_y[i][j])
                else:
                    # Measured flux is negative, so plot at an artifically
                    # low point, and make error bars correct
                    data_y_adjusted.append(1e-99)

            # Set scale and include error bars
            ax.set_xscale("log")
            ax.set_yscale("log")
            limtest = ax.get_ylim()
            
            if limtest[0] < 1e-99:
                plt.ylim(y_low / 10., y_high * 10)

            ax.errorbar(self.data_x[i], data_y_adjusted, 
                        yerr=[lower_errors, higher_errors],
                        ecolor=colors[i],
                        fmt=next(marker),
                        elinewidth=0.5,
                        color=colors[i],
                        ms=4,
                        capsize=2)
        ax.legend(loc=2)
        ax.set_xlabel("Wavelength (microns)", fontsize=12)
        ax.set_ylabel(r"$\nu$F$_\nu$", fontsize=12)
        ax.set_title(f"{self.packed_name}", loc="left")
        fig.savefig(f"../{self.directory}/general_plots/bestfit_SED.png", dpi=1000)
        fig.savefig(f"../{self.directory}/general_plots/bestfit_SED.pdf", dpi=1000)
        plt.close(fig)


    def clear_plot_directories(self):
        """
        Clears all plots from this object's current directory.
        """
        print("Clearing all plotting directories.")
        for file in os.listdir(f"../{self.directory}/general_plots/"):
            os.remove(f"../{self.directory}/general_plots/" + file)
        segment_dirs = ["diameter_albedo", "diameter_gamma", "diameter_period"]
        for dir in segment_dirs:
            try:
                for file in os.listdir(f"../{self.directory}/{dir}_segments/"):
                    os.remove(f"../{self.directory}/{dir}_segments/" + file)
            except:
                pass

    def generate_histograms(self):
        """
        Generates histogram plots for diameter, albedo, gamma, and period solutions.
        """
        
        print("Generating histograms.")
        histogram_template(self.directory, self.packed_name, self.diameters, "Diameter", 
                           self.is_triaxial, "km")
        histogram_template(self.directory, self.packed_name, self.albedos, "Albedo", 
                           self.is_triaxial)
        histogram_template(self.directory, self.packed_name, self.gammas,
                           "Thermal Inertia", self.is_triaxial, 
                           r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        if not self.fixed:
            histogram_template(self.directory, self.packed_name, self.periods, "Period", 
                               self.is_triaxial, "hr")

    
    def generate_chi_scatterplots(self):
        """
        Generates fit chi squared labeled scatterplots of all output parameters 
        in all axes configurations.
        """
        print("Generating chi^2 labeled scatterplots.")
        chi_scatterplot_template(self.directory, self.packed_name, self.diameters, 
                                 self.albedos, self.chis, "Diameter", "Albedo", 
                                 self.is_triaxial, unit_x="km")
        chi_scatterplot_template(self.directory, self.packed_name, self.diameters, 
                                 self.gammas, self.chis, "Diameter", "Thermal Inertia", 
                                 self.is_triaxial, unit_x="km", 
                                 unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        chi_scatterplot_template(self.directory, self.packed_name, self.gammas, 
                                 self.albedos, self.chis, "Thermal Inertia", "Albedo", 
                                 self.is_triaxial, unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        if not self.fixed:
            chi_scatterplot_template(self.directory, self.packed_name, self.diameters, 
                                     self.periods, self.chis, "Diameter", "Period", 
                                     self.is_triaxial, unit_x="km", unit_y="hr")
            chi_scatterplot_template(self.directory, self.packed_name, self.gammas, 
                                     self.periods, self.chis, "Thermal Inertia", 
                                     "Period", self.is_triaxial, 
                                     unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$", unit_y="hr")
            chi_scatterplot_template(self.directory, self.packed_name, self.albedos, 
                                     self.periods, self.chis, "Albedo", "Period", 
                                     self.is_triaxial, unit_y="hr")
        
    
    def generate_hexbins(self):
        """
        Generates hexbin density plots for all output parameters in all axes 
        configurations.
        """
        print("Generating hexbin plots.")
        hexbin_template(self.directory, self.packed_name, self.diameters, self.albedos,
                        "Diameter", "Albedo", self.is_triaxial, unit_x="km")
        hexbin_template(self.directory, self.packed_name, self.diameters, self.gammas,
                        "Diameter", "Thermal Inertia", self.is_triaxial, unit_x="km", 
                        unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        hexbin_template(self.directory, self.packed_name, self.gammas, self.albedos,
                            "Thermal Inertia", "Albedo", self.is_triaxial, 
                            unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        if not self.fixed:
            hexbin_template(self.directory, self.packed_name, self.diameters, 
                            self.periods, "Diameter", "Period", self.is_triaxial, 
                            unit_x="km", unit_y="hr")
            hexbin_template(self.directory, self.packed_name, self.gammas, self.periods,
                            "Thermal Inertia", "Period", self.is_triaxial, 
                            unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$", unit_y="hr")
            hexbin_template(self.directory, self.packed_name, self.albedos, self.periods,
                            "Albedo", "Period", self.is_triaxial, unit_y="hr")


    def generate_chi_plots(self):
        """
        Generates scatterplots of all output parameters versus the fit chi squared 
        values for each invididual MCMC solution.
        """
        print("Generating chi^2 plots.")
        chi_plot_template(self.directory, self.packed_name, self.diameters, self.chis, 
                          "Diameter", self.is_triaxial, "km")
        chi_plot_template(self.directory, self.packed_name, self.albedos, self.chis, 
                          "Visual Albedo", self.is_triaxial)
        chi_plot_template(self.directory, self.packed_name, self.gammas, self.chis, 
                          "Thermal Inertia", self.is_triaxial, 
                          r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        if not self.fixed:
            chi_plot_template(self.directory, self.packed_name, self.periods, self.chis, 
                              "Period", self.is_triaxial, "hr")



