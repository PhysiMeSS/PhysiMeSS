# PhysiMeSS
PhysiMeSS (PhysiCell Microenvironment Structure Simulation) is a PhysiCell add-on which allows users to simulate ECM components as agents. 

**Version:** 1.0.1
**Base PhysiCell Version:** 1.14.0

**Release date:** 16 September 2024

## Paper
PhysiMeSS paper is available on Gigabyte : [https://gigabytejournal.com/articles/136](https://gigabytejournal.com/articles/136). 

DOI: [10.46471/gigabyte.136](https://doi.org/10.46471/gigabyte.136).

## 3D version
The 3D version of PhysiMeSS (fibres initialised, pushed, rotated and degraded in three dimensions, and the 3D examples) is the work of Ambre Jaeger ([@ambrejaeger](https://github.com/ambrejaeger)).

## Changes of the model
The mechanics of the interaction between cells and fibres has changed since version 1.0.1. Every change, with its justification and its effect on the results, is listed in [addons/PhysiMeSS/MODEL_CHANGES.md](addons/PhysiMeSS/MODEL_CHANGES.md).

## Dedicated sample project
PhysiMeSS comes with a dedicated sample project, called **physimess-sample**. To build it, go to the root directory and use : 

```
    make physimess-sample
    make
```
This will load all the files and directories contained in ```./sample_projects/physimess/config``` in the ```./config``` directory, replace the makefile in your PhysiCell directory by the one in ```./sample_projects/physimess``` and replace the ```./custom_modules``` directory by the ```./sample_projects/physimess/custom_modules```
## Pre-loaded examples 
The following example directories are populated in config directory once the **physimess-sample** project is loaded as above. The previous commands compiles and creates the .exe file by default named project in your working directory. To run the following example, you need to specify the location of .xml setting file populated in the config folder. You can as well overwrite the output folder of the simulation set in the .xml setting file directly from the command line. For instance in Linux, for the Fibre_Degradation_3D:

```
./project -s "./config/Fibre_Degradation_3D/mymodel_fibredegradation3D.xml" -o "./output/Fibre_Degradation_3D"
```
After ```-o``` or ```--output```: specify the path to the desired output folder. If not specified the output folder is the one specified in your setting file
After ```-s``` or ```--settings```: specify the path to your setting .xml file. If not specified the default is ```./config/PhysiCell_settings.xml```
### Examples in 2D
#### Fibre_Initialisation
The directory Fibre_Initialisation contains simple examples in which you can initialise ECM fibres in the domain. Fibres are cylindrical agents described by their centre, radius, length and orientation. The centre of each fibre is prescribed either from a csv file or at random (as per cells in PhysiCell). The other attributes can be altered via user parameters in the xml or GUI. The following default parameters are found in ```mymodel_initialisation.xml```

```
number_of_fibres = 2000
anisotropic_fibres = false
fibre_length = 75.0 (microns)
length_normdist_sd  = 0.0 (microns)
fibre_radius = 2.0 (microns)
fibre_angle = 0.0 (radians)
angle_normdist_sd = 0.0 (radians)
```

Using these parameters you can set up a domain with 2000 fibres randomly positioned and randomly aligned with a length of 75 microns and radius of 2 microns. Note, since we remove any fibres which overlap the boundaries of the domain, 1931 fibres remain after initialisation. Although fibres are cylinders, they are visualised in the domain and plot legend as lines. Note, agent names with the following strings are presumed to be fibres: **ecm**; **fibre**; **fiber**; **rod**; **matrix**.

A second xml file ```mymodel_initialisation_maze.xml``` along with the csv ```initialfibres.csv``` allows you to create a maze with horizontal and vertical fibres. Agent names **fibre_horizontal** and **fibre_vertical** are reserved for creating horizontal and vertical fibres, respectively. In this way fibre agents can be used to create walls in your domain.

#### Fibre_Degradation 
The directory Fibre_Degradation contains examples which show how fibre degradation affects simulations of cell migration and cell proliferation.

* The xml file ```mymodel_fibre_degradation.xml``` with csv file ```cells_and_fibres_attractant.csv``` simulates the migration of a single cell towards an attractant through a mesh of fibres. Fibres are assigned a position, radius and length but their orientation is random. By turning fibre degradation on (or off) and adjusting parameters you can control whether the cell navigates to the attractant.

Default parameters are:

```
fibre_length = 40.0 (microns)
fibre_degradation = true 
fibre_deg_rate = 0.01 (1/min)
fibre_stuck = 10.0 (mechanics timesteps)
```

Note the parameter ```fibre_stuck``` determines how many mechanics time steps a cell needs to have been stuck before it can possibly degrade a fibre at the rate ```fibre_deg_rate```. When a fibre is degraded it is immediately removed from the domain.

* The xml file ```mymodel_matrix_degradation.xml``` with csv file ```cells_and_fibres.csv``` simulates cell proliferation and the forming of a cell mass within a mesh of fibres. The fibrous mesh is as above. By turning fibre degradation on (or off) and adjusting parameters you can control the development of the growing mass of cells.


#### Cell_Fibre_Mechanics
The directory Cell_Fibre_Mechanics contains examples which demonstrate cell-fibre mechanics. 

* The xml file ```mymodel_fibremaze.xml``` with csv file  ```fibre_maze.csv``` simulates a single cell moving within a maze made of fibres towards an attractant secreting a nutrient. By adjusting parameters you can control whether the cell navigates to the attractant.

Default parameters are:

```
fibre_length = 60.0 (microns)
vel_adhesion = 0.6 
vel_contact = 0.1
cell_velocity_max = 1.0
```

A fibre within reach of a cell accelerates it along the fibre, in the direction in which the cell moves, by up to ```vel_adhesion``` (less as the speed of the cell approaches ```cell_velocity_max```, and not at all above it), and reduces by up to ```vel_contact``` (between 0 and 1) the part of its motility which is across the fibre.

* The xml file ```mymodel_pushing.xml``` with csv file  ```snowplough.csv``` simulates a single cell pushing a single free/non-crosslinked fibre out of the way to access an attractant. By turning fibre pushing on (or off) you can control whether the cell can push fibres. Note, since cells can only push non-crosslinked fibres by increasing the ```fibre_length``` in this example and thereby creating crosslinks you can prevent cells from pushing fibres.

Default parameters are:

```
fibre_length = 40.0 (microns)
fibre_pushing = true
```

* The xml file ```mymodel_rotating.xml``` with csv file ```snowplough.csv``` simulates a single cell rotating free/non-crosslinked fibres to access an attractant. By turning fibre rotating on (or off) you can control whether the cell can push fibres. Note the parameter ```fibre_sticky``` modulates how much a cell can rotate a fibre. A cell rotates a fibre about its centre, or about its crosslink point if it has one, in two ways: by contact (the repulsion between the cell and the fibre, as for pushing) and by traction (a migrating cell drags the fibre with its motility velocity multiplied by ```fibre_sticky```). ```fibre_sticky``` is a property of the cell: it represents its adhesion to fibres, and can differ between cell types. With 1 the cell drags the fibre with the force which would move it at the speed of the cell; with 0 the cell does not adhere and only the contact rotates the fibre. A fibre is harder to rotate about a crosslink far from its centre.

Default parameters are:

```
fibre_length = 40.0 (microns)
fibre_sticky = 1.0
fibre_rotation = true
```

* The xml file ```mymodel_hinge.xml``` with csv file ```hinge.csv``` simulates two crosslinked fibres being rotated at their crosslink point (hinge) by a single cell in order for the cell to navigated towards an attractant. By turning fibre rotating on (or off)  and adjusting parameters you can control whether the cell successfully navigates towards the attractant.

### Examples in 3D
In the following examples ```use_2D``` is set to false in the xml file: fibres are oriented, pushed and rotated in three dimensions. With anisotropic fibres, the orientation is set by the two angles ```fibre_angle``` (in the xy plane) and ```fibre_angle_phi``` (from the z axis), in radians. When ```fibre_angle_phi``` is not defined, pi/2 is used: the fibres are in the xy plane. As in 2D, fibres of cell definitions named ```fibre_vertical``` and ```fibre_horizontal``` are always along y and x, whatever their angles.

#### Cell_Fibre_Mechanics_3D
The directory Cell_Fibre_Mechanics_3D contains the 3D counterparts of the cell-fibre mechanics examples.

* The xml file ```mymodel_pushing.xml``` with csv file ```push.csv``` simulates a single cell pushing a free fibre on its way to an attractant.

* The xml file ```mymodel_rotating.xml``` with csv file ```rotate.csv``` simulates a single cell rotating a free fibre on its way to an attractant.

* The xml files ```mymodel_squish.xml``` and ```mymodel_squish_2.xml``` with csv files ```squish.csv``` and ```squish_2.csv``` simulate two cells moving in opposite directions, each towards its own attractant, with a single fibre between them which they both push and rotate.

* The xml file ```mymodel_hinge.xml``` with csv file ```hinge.csv``` simulates a single cell moving towards an attractant placed out of the plane of two crosslinked fibres. The fibre in its way is pushed and rotated out of the plane, about its crosslink with the other fibre.

* The xml file ```barrier_fibre.xml``` with csv file ```barrier.csv``` simulates a single cell moving towards an attractant through a barrier of fibres placed in a 3D domain. Fibre pushing and rotation are off by default.

#### Fibre_Degradation_3D
The directory Fibre_Degradation_3D contains examples of a growing mass of cells within a 3D mesh of fibres. The cell is placed from the csv file ```cells.csv``` while the fibres are created at random: in the ```initial_conditions``` of the xml file,

```
<fibres_from_file>
    <enable>false</enable>
    <ID>1</ID>
    <relative_fibre_volume>0.001</relative_fibre_volume>
</fibres_from_file>
```

creates fibres of the cell definition ```ID``` until they occupy the fraction ```relative_fibre_volume``` of the domain. When ```enable``` is true (the default when the element is absent), fibres are instead read from the csv file as in the 2D examples. Fibres which overlap a cell at initialisation are removed.

* The xml file ```mymodel_fibredegradation3D.xml``` is the base example: a dense mesh (```relative_fibre_volume``` of 0.0094) of loosely aligned fibres of length 50 microns, degraded by non motile cells.

* The xml file ```deg3D_dlow_ani.xml``` has longer (100 microns), more aligned fibres at low density (```relative_fibre_volume``` of 0.001).

* The xml file ```deg3Dmove_dlow_ani.xml``` is the same with motile cells which push and rotate the fibres.

Other conditions can be explored from these files by changing ```relative_fibre_volume``` (fibre density), ```anisotropic_fibres``` and ```angle_normdist_sd``` (alignment of the fibres), ```fibre_degradation```, and ```fibre_custom_degradation``` (degradation depending on pressure).

## Visualization in ParaView
Cells and fibres can be displayed in 3D with [ParaView](https://www.paraview.org/download/): cells as spheres, fibres as cylinders. This requires Python 3.6+ and ParaView (tested with ParaView 6.2).

**1. Save the full data.** Check that ```<full_data>``` ```<enable>``` is true in your .xml settings file, so that the simulation writes .mat and .xml files in the output folder.

**2. Convert the output.** From the PhysiCell root directory,
```
    make vtu output="./output"
```
converts the output folder given in ```output``` (```./output``` by default). The first time, it installs the required python packages (numpy, scipy and vtk) in a virtual environment (physimessenv). The script ```paraview_viewer.py``` (in your config folder) then writes, in the output folder, one .vtu file per saved time step and a file ```simulation.pvd``` which lists them. Each agent is a point, located at the centre of the cell or fibre, with all the data saved by PhysiCell (cell type, radius, orientation, custom data, ...) and a flag ```is_fibre```.

**3. Display cells and fibres.** In ParaView:
1. File > Open, select ```simulation.pvd``` in your output folder, and click on Apply in the properties panel. You should see points.
2. View > Python Shell, click on Run Script and select ```cells_and_fibres.py``` in your config folder. You should see spheres and cylinders, coloured by cell type.
3. Use the animation controls of the toolbar to go through the time steps.

Cells are drawn with their radius, and fibres with the length and radius of their cell definition (```fibre_length```, ```fibre_radius```) along their orientation. Fibres are the agents whose cell definition name contains ecm, matrix, fiber, fibre or rod, as in the simulation. The script creates the filters ```cells```, ```spheres```, ```fibres```, ```fibre_axes``` and ```cylinders``` in the pipeline browser; select one of them to change its colour or its properties.

**Images without the graphical interface.** With the ```pvpython``` program installed with ParaView:
```
    pvpython ./config/cells_and_fibres.py ./output/simulation.pvd --png image.png
    pvpython ./config/cells_and_fibres.py ./output/simulation.pvd --frames ./output/frames
```
The first command saves an image of the last time step (or of the time step given with ```--step```), the second one image per time step.

## Future work
The next development of PhysiMeSS concerns the mechanics of the fibre network. In this version a crosslink between two fibres is a contact: a fibre with one crosslink rotates about it, a fibre with two or more crosslinks does not move, and no force is transmitted from one fibre to another. The plan is to make crosslinks persistent bonds with a force between the two fibres, so that cells can deform and remodel a network of fibres instead of only degrading it. The current rules and their limits are described in [addons/PhysiMeSS/MODEL_CHANGES.md](addons/PhysiMeSS/MODEL_CHANGES.md).
