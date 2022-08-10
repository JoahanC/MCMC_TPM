"""
This is the driver files for the mass application of the Asteroids class for
analyis of the MCMC results of correctly stored objects.
"""
import os
from asteroids import Asteroid


files = os.listdir("../")
dirs = []
for file in files:
    if os.path.isdir(f"../{file}"):
        dirs.append(file)
dirs.remove(".git")
dirs.remove("best_fits")
dirs.remove("triaxial_02100")
dirs.remove("02100")


objects = []
print("Generating Asteroid objects.")
for idx, dir in enumerate(dirs):
    objects.append(Asteroid(dir))
print("Finished generating Asteroid objects.\n")

print("Generating histogram plots")
for object in objects:
    object.generate_histograms()

