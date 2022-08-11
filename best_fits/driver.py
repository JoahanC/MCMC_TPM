"""
This is the driver files for the mass application of the Asteroids class for
analyis of the MCMC results of correctly stored objects.
"""
import os
from asteroids import Asteroid
from dualplotter import DualPlotter


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

for idx, object in enumerate(objects):
    if not object.is_triaxial:
        print(f"Processing object {idx + 1}/{len(objects)}: {object.packed_name}")
    elif object.is_triaxial:
        print(f"Processing object {idx + 1}/{len(objects)}: Triaxial {object.packed_name}")
    #object.clear_plot_directories()
    #object.generate_SED_plot()
    object.generate_histograms()
    #object.generate_chi_scatterplots()
    #object.generate_hexbins()

"""pairings = {}
for idx, object in enumerate(objects):
    if object.packed_name not in pairings:
        pairings[object.packed_name] = [object]
    else:
        pairings[object.packed_name].append(object)

dual_objects = []
for code in pairings:
    if pairings[code][0].is_triaxial:
        spherical_ob = pairings[code][1]
        triaxial_ob = pairings[code][0]
        curr_ob = DualPlotter(spherical_ob, triaxial_ob)
        dual_objects.append(curr_ob)
    else:
        spherical_ob = pairings[code][0]
        triaxial_ob = pairings[code][1]
        curr_ob = DualPlotter(spherical_ob, triaxial_ob)
        dual_objects.append(curr_ob)

for object in dual_objects:
    object.generate_histograms()"""
    