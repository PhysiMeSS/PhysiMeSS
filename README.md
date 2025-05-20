# PhysiMeSS
PhysiMeSS (PhysiCell Microenvironment Structure Simulation) is a PhysiCell add-on which allows users to simulate ECM components as agents. 

**Version:** 1.0.1
**Base PhysiCell Version:** 1.14.0

**Release date:** 16 September 2024

## Paper
PhysiMeSS paper is available on Gigabyte : [https://gigabytejournal.com/articles/136](https://gigabytejournal.com/articles/136). 

DOI: [10.46471/gigabyte.136](https://doi.org/10.46471/gigabyte.136).

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
./project -s "./config/Fibre_Degradation_3D/PhysiCell_settings.xml" -o "./output/Fibre_Degradation"
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

* The xml file ```mymodel_pushing.xml``` with csv file  ```snowplough.csv``` simulates a single cell pushing a single free/non-crosslinked fibre out of the way to access an attractant. By turning fibre pushing on (or off) you can control whether the cell can push fibres. Note, since cells can only push non-crosslinked fibres by increasing the ```fibre_length``` in this example and thereby creating crosslinks you can prevent cells from pushing fibres.

Default parameters are:

```
fibre_length = 40.0 (microns)
fibre_pushing = true
```

* The xml file ```mymodel_rotating.xml``` with csv file ```snowplough.csv``` simulates a single cell rotating free/non-crosslinked fibres to access an attractant. By turning fibre rotating on (or off) you can control whether the cell can push fibres. Note the parameter ```fibre_sticky``` modulates how much a cell can rotate a fibre. Values less than 1.0 decrease how much the cell rotates the fibre per timestep and values greater than 1.0 increase how much the cell rotates the fibre per timestep.

Default parameters are:

```
fibre_length = 40.0 (microns)
fibre_sticky = 1.0
fibre_rotation = true
```

* The xml file ```mymodel_hinge.xml``` with csv file ```hinge.csv``` simulates two crosslinked fibres being rotated at their crosslink point (hinge) by a single cell in order for the cell to navigated towards an attractant. By turning fibre rotating on (or off)  and adjusting parameters you can control whether the cell successfully navigates towards the attractant.

### Examples in 3D

## Requirements for 3D vizualisation 
If you wish to generate file using VTK readable directly in parawiew, you will need: 
- Python 3.6+
- Paraview : https://www.paraview.org/download/ 

Tests were only performed using Paraview 5.11.0.
Paraview documentation and tutorial are available [here](https://docs.paraview.org/en/v5.11.0/UsersGuide/index.html).

## Visualization in Paraview
When running your simulation check that save full data is enabled in your .xml setting file. This will generate .mat and .xml files during your simulation that will be used by the ```paraview_viewer.py``` (located ```./sample_projects/physimess/config```) script to create .vtu files and a .pvd file that can be read in paraview. To run the script from the command line do:
```
    make vtu output="./output"
```
Replace ``` ./output``` by your output folder where the .mat and .xml files are located. If no argument is provided, ```"./output"``` will be set as default. 

This command will install the required pyhton packages (scipy, numpy and vtk) in a virtual environment (physimessenv), and activate the environment before running the python script.

Then open paraview. 

**Brief overview of paraview functioning** 

The python script paraview_viewer.py generates .vtu files containing points (coordinate in space) associated to a number of properties that are saved during PhysiCell simulation when in the simulation settings file, \<full_data>\<enable> is set to true. This points correspond to the center of cells and fibres. In paraview, we can apply filters on those points to visualize them as different geometrical shape of various sizes and so on. Paraview works with layer. On the input file we are going to apply multiple filters on points based on the properties associated to them. These filters to apply on input data can be saved in a state file (.pvsm) and thus reuse on other input data.

**Load data and a state file in paraview**

In paraview: 
1. In the top menu bar, go to File > Load State ... select the appropriate .pvsm file (it should be loade in your config folder). 
2. In the next window select "Choose file name in the scroll menu. Then select the .vtu files you wish to visualize in your output folder.
3. Click on Apply at the top of the properties menu. You should see spheres and cylinders.

**Discriminate appropriately between spheres and cylinders**

By default, agents of type 0 in your simulation are displayed as sphere and agent of type 1 are displayed as cylinders.The type is the number and name associated with the cell definition in your .xml setting file, if you are unsure of how your types are set, you can check the value in the setting files that is saved in your output folder when you run a simulation. If in your simulation the type differ from this default:
1. In the pipeline browser, select type_0 (this is a Threshold filter that select all agents with cell type equal to 0).
2. If your cells are of cell type 2 for instance, in the properties menu, modify Lower Threshold and Upper Threshold to 2.
3. Click on Apply. You should see spheres appear.

Proceed the same way if your fibres are not of type 1.

**Modify the default size parameters of spheres and cylinders:**

By default, the sphere radius is set as 8.413 according to PhysiCell parameters, the radius of fibres is set at 2 and the length is set at 40. If the parameters differ in your simulation:
1. In the pipeline browser, select spheres or cylinders
2. In the properties menu, in the section Glyph Source, you can modify radius or height (which corresponds to length) as you wish. 

