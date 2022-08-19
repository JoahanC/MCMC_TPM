"""
This file serves as a tutorial for using the DualPlotter class.
"""
from asteroids import Asteroid
from dualplotter import DualPlotter


"""
The dualplotter class is able to generate comparison plots between objects modeled
by different TPM shape types. To initialize a DualPlotter object, you must first
generate two Asteroid objects of asteroids sharing the same packed name but being
located in separate folders. We'll use Cacus as an example. The spherical model
object must be first in the constructor
"""

spherical_cacus = Asteroid("G1989")
triaxial_cacus = Asteroid("triaxial_G1989")
dual = DualPlotter(spherical_cacus, triaxial_cacus)

"""
In a similar fashion to the Asteroid object, there come plotting capabilities within
the DualPlotter class. We can run a few of these below. The end results will be located
in the comparison plots folder, in the corresponding folder containing Cacus' packed
MPC name.
"""

DualPlotter.generate_histograms()
DualPlotter.generate_chi_scatterplots()
DualPlotter.generate_hexbins()
DualPlotter.generate_chiplots()

"""
The DualPlotter class also comes with its own directory clearing method, which will
clear the appropriate folder in the comparison plots folder.
"""

DualPlotter.clear_directory()