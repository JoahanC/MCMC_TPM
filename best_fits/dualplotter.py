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
            #print(spherical_asteroid.is_triaxial)
            raise ValueError("Spherical Asteroid is not spherical!")
        if not triaxial_asteroid.is_triaxial:
            #print(triaxial_asteroid.is_triaxial)
            raise ValueError("Triaxial Asteroid is not triaxial!")
        
        self.spherical_asteroid = spherical_asteroid
        self.triaxial_asteroid = triaxial_asteroid
        self.packed_name = spherical_asteroid.packed_name
        if self.packed_name not in os.listdir("./comparison_plots/"):
            os.mkdir(f"./comparison_plots/{self.packed_name}")


    def generate_histograms(self):
        comparison_histogram_template(self.packed_name, self.spherical_asteroid.diameters,
                                      self.triaxial_asteroid.diameters, "Diameter", "km")

    
        
