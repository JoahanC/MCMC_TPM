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

        packed_name : str
            The packed MPC designation for the object.

        """
        self.directory = directory
        self.is_triaxial = False
        self.packed_name = directory

        if "triaxial_" in directory:
            self.packed_name = directory.replace("triaxial_", '')
            self.is_triaxial = True

        best_fit_outputs = open(f"{self.packed_name}.txt", 'r')
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
                self.albedo_16th_percentile = float(p_V_vals[1])
                self.albedo_84th_percentile = float(p_V_vals[1])
                self.albedo_median = float(p_V_vals[2])
                self.albedo_16th_percentile = float(p_V_vals[3])
                self.albedo_84th_percentile = float(p_V_vals[3])

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

            if "pole peak at = " in line:
                pole_peak_vals = line.split()
                self.pole_peak_ra = float(pole_peak_vals[4])
                self.pole_peak_dec = float(pole_peak_vals[5])

            if "mean pole at = " in line:
                mean_pole_vals = line.split()
                self.pole_mean_ra = float(mean_pole_vals[4])
                self.pole_mean_dec = float(mean_pole_vals[5])
                self.pole_bar = float(mean_pole_vals[7])
        best_fit_outputs.close()

        SED_plotters = self.retrieve_SED_data()

        self.esed = SED_plotters[0]
        self.wavelengths = SED_plotters[1]
        self.datelabels = SED_plotters[2]
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


    def retrieve_SED_data(self):
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

        Parameters
        ----------
        object : str
            The folder name of the object being analyzed.

        Returns
        -------
        A tuple containing the diameters[0], albedos[1], gammas[2], 
        chi^2 values[3], and period[4] values for each solution point. 
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

        physical_chars = [diameters, albedos, gammas, chis, periods]
        return physical_chars

    
    def generate_SED_plot(self):
        colors = ["#3498db", "#229954", "#c0392b", "#8e44ad", "#f1c40f", "#ec7063", "#34495e", "#6e2c00"]

        fig, ax = plt.subplots()
        #plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
        #ax = plt.gca()
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
        histogram_template(self.directory, self.packed_name, self.diameters, "Diameter", self.is_triaxial, "km")
        histogram_template(self.directory, self.packed_name, self.albedos, "Albedo", self.is_triaxial)
        histogram_template(self.directory, self.packed_name, self.gammas, "Thermal Inertia", self.is_triaxial, r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        if not self.fixed:
            histogram_template(self.directory, self.packed_name, self.periods, "Period", self.is_triaxial, "hr")

    
    def generate_chi_scatterplots(self):
        chi_scatterplot_template(self.directory, self.packed_name, self.diameters, self.albedos,
                                 self.chis, "Diameter", "Albedo", self.is_triaxial, unit_x="km")
        chi_scatterplot_template(self.directory, self.packed_name, self.diameters, self.gammas,
                                 self.chis, "Diameter", "Thermal Inertia", self.is_triaxial, unit_x="km", unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        chi_scatterplot_template(self.directory, self.packed_name, self.gammas, self.albedos,
                                    self.chis, "Thermal Inertia", "Albedo", self.is_triaxial, unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        if not self.fixed:
            chi_scatterplot_template(self.directory, self.packed_name, self.diameters, self.periods,
                                    self.chis, "Diameter", "Period", self.is_triaxial, unit_x="km", unit_y="hr")
            chi_scatterplot_template(self.directory, self.packed_name, self.gammas, self.periods,
                                    self.chis, "Thermal Inertia", "Period", self.is_triaxial, unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$", unit_y="hr")
            chi_scatterplot_template(self.directory, self.packed_name, self.albedos, self.periods,
                                    self.chis, "Albedo", "Period", self.is_triaxial, unit_y="hr")
        
    
    def generate_hexbins(self):
        hexbin_template(self.directory, self.packed_name, self.diameters, self.albedos,
                        "Diameter", "Albedo", self.is_triaxial, unit_x="km")
        hexbin_template(self.directory, self.packed_name, self.diameters, self.gammas,
                        "Diameter", "Thermal Inertia", self.is_triaxial, unit_x="km", unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        hexbin_template(self.directory, self.packed_name, self.gammas, self.albedos,
                            "Thermal Inertia", "Albedo", self.is_triaxial, unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        if not self.fixed:
            hexbin_template(self.directory, self.packed_name, self.diameters, self.periods,
                            "Diameter", "Period", self.is_triaxial, unit_x="km", unit_y="hr")
            hexbin_template(self.directory, self.packed_name, self.gammas, self.periods,
                            "Thermal Inertia", "Period", self.is_triaxial, unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$", unit_y="hr")
            hexbin_template(self.directory, self.packed_name, self.albedos, self.periods,
                            "Albedo", "Period", self.is_triaxial, unit_y="hr")


    def generate_chi_plots(self):
        chi_plot_template(self.directory, self.packed_name, self.diameters, self.chis, "Diameter", self.is_triaxial, "km")
        chi_plot_template(self.directory, self.packed_name, self.albedos, self.chis, "Visual Albedo", self.is_triaxial)
        chi_plot_template(self.directory, self.packed_name, self.gammas, self.chis, "Thermal Inertia", self.is_triaxial, r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        if not self.fixed:
            chi_plot_template(self.directory, self.packed_name, self.periods, self.chis, "Period", self.is_triaxial, "hr")



