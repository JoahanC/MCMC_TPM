from asteroids import Asteroid
from helpers import *
import os


class DualPlotter:
    """
    The Dual Plotter class allows for comparison plots to be made from Asteroid
    objects that share the same packed_name but were run with spherical and 
    triaxial models seperately.
    """
    def __init__(self, spherical_asteroid, triaxial_asteroid):
        """
        Constructor for the DualPlotter class.

        Parameters
        ----------
        spherical_asteroid : Asteroid
            An Asteroid object run with the spherical variation of the MCMC code.

        triaxial_asteroid : Asteroid
            An Asteroid object run with the triaxial varation of the MCMC code.
        """

        if spherical_asteroid.packed_name != triaxial_asteroid.packed_name:
            raise ValueError("Both Asteroids must share the same packed_name!\n" +
                             f"{spherical_asteroid.packed_name} and " +
                             f"{triaxial_asteroid.packed_name} do not match!")
        if spherical_asteroid.is_triaxial:
            raise ValueError("Spherical Asteroid is not spherical!")
        if not triaxial_asteroid.is_triaxial:
            raise ValueError("Triaxial Asteroid is not triaxial!")
        
        self.spherical_asteroid = spherical_asteroid
        self.triaxial_asteroid = triaxial_asteroid
        self.packed_name = spherical_asteroid.packed_name
        if self.packed_name not in os.listdir("./comparison_plots/"):
            os.mkdir(f"./comparison_plots/{self.packed_name}")


    def clear_directory(self):
        for file in os.listdir(f"./comparison_plots/{self.packed_name}"):
            os.remove(f"./comparison_plots/{self.packed_name}/{file}")


    def generate_histograms(self):
        """
        Generates comparison histograms for all output parameters.
        """
        print("Generating histograms.")
        if "histograms" not in os.listdir(f"./comparison_plots/{self.packed_name}/"):
            os.mkdir(f"./comparison_plots/{self.packed_name}/histograms/")
        comparison_histogram_template(self.packed_name, self.spherical_asteroid.diameters,
                                      self.triaxial_asteroid.diameters, "Diameter", "km")
        comparison_histogram_template(self.packed_name, self.spherical_asteroid.albedos,
                                      self.triaxial_asteroid.albedos, "Visible Albedo")
        comparison_histogram_template(self.packed_name, self.spherical_asteroid.gammas,
                                      self.triaxial_asteroid.gammas, "Thermal Inertia", r"$J~m^{-2}~s^{-0.5}~K^{-1})$")


    def generate_chi_scatterplots(self):
        """
        Generates fit chi squared labeled comparison scatterplots for all output 
        parameters in all possible configurations.
        """
        print("Generating chi labeled scatterplots.")
        if "scatterplots" not in os.listdir(f"./comparison_plots/{self.packed_name}/"):
            os.mkdir(f"./comparison_plots/{self.packed_name}/scatterplots/")
        comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.diameters, self.spherical_asteroid.albedos, self.spherical_asteroid.chis,
                        self.triaxial_asteroid.diameters, self.triaxial_asteroid.albedos, self.triaxial_asteroid.chis,
                        "Diameter", "Visible Albedo", unit_x="km")
        comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.diameters, self.spherical_asteroid.gammas, self.spherical_asteroid.chis,
                        self.triaxial_asteroid.diameters, self.triaxial_asteroid.gammas, self.triaxial_asteroid.chis,
                        "Diameter", "Thermal Inertia", unit_x="km", unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.gammas, self.spherical_asteroid.albedos, self.spherical_asteroid.chis,
                        self.triaxial_asteroid.gammas, self.triaxial_asteroid.albedos, self.triaxial_asteroid.chis,
                        "Thermal Inertia", "Visible Albedo", unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")     
        comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.albedos, self.spherical_asteroid.diameters, self.spherical_asteroid.chis,
                        self.triaxial_asteroid.albedos, self.triaxial_asteroid.diameters, self.triaxial_asteroid.chis,
                        "Visible Albedo", "Diameter", unit_y="km")
        comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.gammas, self.spherical_asteroid.diameters, self.spherical_asteroid.chis,
                        self.triaxial_asteroid.gammas, self.triaxial_asteroid.diameters, self.triaxial_asteroid.chis,
                        "Thermal Inertia", "Diameter", unit_y="km", unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.albedos, self.spherical_asteroid.gammas, self.spherical_asteroid.chis,
                        self.triaxial_asteroid.albedos, self.triaxial_asteroid.gammas, self.triaxial_asteroid.chis,
                        "Visible Albedo", "Thermal Inertia", unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        
        if not self.spherical_asteroid.fixed and not self.triaxial_asteroid.fixed:
            comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.diameters, self.spherical_asteroid.periods, self.spherical_asteroid.chis,
                            self.triaxial_asteroid.diameters, self.triaxial_asteroid.periods, self.triaxial_asteroid.chis,
                            "Diameter", "Period", unit_x="km", unit_y="hr")
            comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.gammas, self.spherical_asteroid.periods, self.spherical_asteroid.chis,
                            self.triaxial_asteroid.gammas, self.triaxial_asteroid.periods, self.triaxial_asteroid.chis,
                            "Thermal Inertia", "Period", unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$", unit_y="hr")
            comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.albedos, self.spherical_asteroid.periods, self.spherical_asteroid.chis,
                            self.triaxial_asteroid.albedos, self.triaxial_asteroid.periods, self.triaxial_asteroid.chis,
                            "Visible Albedo", "Period", unit_y="hr")
            comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.periods, self.spherical_asteroid.diameters, self.spherical_asteroid.chis,
                            self.triaxial_asteroid.periods, self.triaxial_asteroid.diameters, self.triaxial_asteroid.chis,
                            "Period", "Diameter", unit_y="km", unit_x="hr")
            comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.periods, self.spherical_asteroid.gammas, self.spherical_asteroid.chis,
                            self.triaxial_asteroid.periods, self.triaxial_asteroid.gammas, self.triaxial_asteroid.chis,
                            "Period", "Thermal Inertia", unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$", unit_x="hr")
            comparison_scatterplot_template(self.packed_name, self.spherical_asteroid.periods, self.spherical_asteroid.albedos, self.spherical_asteroid.chis,
                            self.triaxial_asteroid.periods, self.triaxial_asteroid.albedos, self.triaxial_asteroid.chis,
                            "Period", "Visible Albedo", unit_x="hr")


    def generate_hexbins(self):
        """
        Generates comparison hexbin plots for all output parameters in all possible
        configurations.
        """
        print("Generating hexbin plots.")
        if "hexbins" not in os.listdir(f"./comparison_plots/{self.packed_name}/"):
            os.mkdir(f"./comparison_plots/{self.packed_name}/hexbins/")
        comparison_hexbin_template(self.packed_name, self.spherical_asteroid.diameters, self.spherical_asteroid.albedos,
                        self.triaxial_asteroid.diameters, self.triaxial_asteroid.albedos,
                        "Diameter", "Visible Albedo", unit_x="km")
        comparison_hexbin_template(self.packed_name, self.spherical_asteroid.diameters, self.spherical_asteroid.gammas,
                        self.triaxial_asteroid.diameters, self.triaxial_asteroid.gammas,
                        "Diameter", "Thermal Inertia", unit_x="km", unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        comparison_hexbin_template(self.packed_name, self.spherical_asteroid.gammas, self.spherical_asteroid.albedos,
                        self.triaxial_asteroid.gammas, self.triaxial_asteroid.albedos,
                        "Thermal Inertia", "Visible Albedo", unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        comparison_hexbin_template(self.packed_name, self.spherical_asteroid.albedos, self.spherical_asteroid.diameters,
                        self.triaxial_asteroid.albedos, self.triaxial_asteroid.diameters,
                        "Visible Albedo", "Diameter", unit_y="km")
        comparison_hexbin_template(self.packed_name, self.spherical_asteroid.gammas, self.spherical_asteroid.diameters,
                        self.triaxial_asteroid.gammas, self.triaxial_asteroid.diameters,
                        "Thermal Inertia", "Diameter", unit_y="km", unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        comparison_hexbin_template(self.packed_name, self.spherical_asteroid.albedos, self.spherical_asteroid.gammas,
                        self.triaxial_asteroid.albedos, self.triaxial_asteroid.gammas,
                        "Visible Albedo", "Thermal Inertia", unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        
        if not self.spherical_asteroid.fixed and not self.triaxial_asteroid.fixed:
            comparison_hexbin_template(self.packed_name, self.spherical_asteroid.diameters, self.spherical_asteroid.periods,
                            self.triaxial_asteroid.diameters, self.triaxial_asteroid.periods,
                            "Diameter", "Period", unit_x="km", unit_y="hr")
            comparison_hexbin_template(self.packed_name, self.spherical_asteroid.gammas, self.spherical_asteroid.periods,
                            self.triaxial_asteroid.gammas, self.triaxial_asteroid.periods,
                            "Thermal Inertia", "Period", unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$", unit_y="hr")
            comparison_hexbin_template(self.packed_name, self.spherical_asteroid.albedos, self.spherical_asteroid.periods,
                            self.triaxial_asteroid.albedos, self.triaxial_asteroid.periods,
                            "Visible Albedo", "Period", unit_y="hr")
            comparison_hexbin_template(self.packed_name, self.spherical_asteroid.periods, self.spherical_asteroid.diameters,
                            self.triaxial_asteroid.periods, self.triaxial_asteroid.diameters,
                            "Period", "Diameter", unit_y="km", unit_x="hr")
            comparison_hexbin_template(self.packed_name, self.spherical_asteroid.periods, self.spherical_asteroid.gammas,
                            self.triaxial_asteroid.periods, self.triaxial_asteroid.gammas,
                            "Period", "Thermal Inertia", unit_y=r"$J~m^{-2}~s^{-0.5}~K^{-1})$", unit_x="hr")
            comparison_hexbin_template(self.packed_name, self.spherical_asteroid.periods, self.spherical_asteroid.albedos,
                            self.triaxial_asteroid.periods, self.triaxial_asteroid.albedos,
                            "Period", "Visible Albedo", unit_x="hr")


    def generate_chiplots(self):
        """
        Generates comparison fit chi squared scatterplots for all output parameters.
        """
        print("Generating chi^2 fit plots.")
        if "chiplots" not in os.listdir(f"./comparison_plots/{self.packed_name}/"):
            os.mkdir(f"./comparison_plots/{self.packed_name}/chiplots/")
        comparison_chi_template(self.packed_name, self.spherical_asteroid.diameters, self.spherical_asteroid.chis, 
                                self.triaxial_asteroid.diameters, self.triaxial_asteroid.chis,
                                "Diameter", unit_x="km")    
        comparison_chi_template(self.packed_name, self.spherical_asteroid.albedos, self.spherical_asteroid.chis, 
                                self.triaxial_asteroid.albedos, self.triaxial_asteroid.chis,
                                "Visual Albedo")
        comparison_chi_template(self.packed_name, self.spherical_asteroid.gammas, self.spherical_asteroid.chis, 
                                self.triaxial_asteroid.gammas, self.triaxial_asteroid.chis,
                                "Thermal Inertia", unit_x=r"$J~m^{-2}~s^{-0.5}~K^{-1})$")
        if not self.spherical_asteroid.fixed and not self.triaxial_asteroid.fixed:
            comparison_chi_template(self.packed_name, self.spherical_asteroid.periods, self.spherical_asteroid.chis, 
                                    self.triaxial_asteroid.periods, self.triaxial_asteroid.chis,
                                    "Period", unit_x="hr")