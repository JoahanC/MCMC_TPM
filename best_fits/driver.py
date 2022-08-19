"""
This is the driver files for the mass application of the Asteroids class for
analyis of the MCMC results of correctly stored objects.
"""
import os
from asteroids import Asteroid
from .helpers import clear_output_files
from dualplotter import DualPlotter

action = input("This file will run the Asteroid and DualPlotter modeling procedure\n" + 
      "for all objects in this repository. Be aware that this may overwrite any\n" + 
      "material currently present in these folders. Would you like to continue?\n")

if action:
    files = os.listdir("../")
    dirs = []
    for file in files:
        if os.path.isdir(f"../{file}"):
            dirs.append(file)
    dirs.remove(".git")
    dirs.remove("best_fits")

    print("All known asteroids:", dirs)
    clusters = set()
    for dir in dirs:
        name = dir.replace("triaxial_", '')
        clusters.add(name)
    print(clusters)
    objects = []
    print("Generating Asteroid objects.")
    for idx, dir in enumerate(dirs):
        object = Asteroid(dir)
        objects.append(object)
    print("Finished generating Asteroid objects.\n")
    objects = []

    for idx, object in enumerate(objects):
        if not object.is_triaxial:
            print(f"Processing object {idx + 1}/{len(objects)}: {object.packed_name}")
        elif object.is_triaxial:
            print(f"Processing object {idx + 1}/{len(objects)}: Triaxial {object.packed_name}")
        object.generate_SED_plot()
        object.generate_histograms()
        object.generate_chi_scatterplots()
        object.generate_hexbins()
        object.generate_chi_plots()

    pairings = {}
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

    for idx, object in enumerate(dual_objects):
        print(f"Processing object {idx + 1}/{len(dual_objects)}: {object.packed_name}")
        object.clear_directory()
        object.generate_histograms()
        object.generate_chi_scatterplots()
        object.generate_hexbins()
        object.generate_chiplots()

clear_output_files()
    