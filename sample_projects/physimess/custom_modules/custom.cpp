/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"



void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 

	cell_defaults.functions.instantiate_cell = instantiate_physimess_cell;	
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = physimess_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
    
    cell_defaults.functions.plot_agent_SVG = fibre_agent_SVG;
	cell_defaults.functions.plot_agent_legend = fibre_agent_legend;
	
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
    
	for (auto* pCD: *getFibreCellDefinitions()){
		Variable crosslink_count = Variable();
		crosslink_count.name = "crosslink_count";
		crosslink_count.value = 0.0;
		pCD->custom_data.add_variable(crosslink_count);
		pCD->functions.instantiate_cell = instantiate_physimess_fibre;
		pCD->functions.plot_agent_SVG = fibre_agent_SVG;
		pCD->functions.plot_agent_legend = fibre_agent_legend;
		
		}

	for (auto* pCD: cell_definitions_by_index){
		if (!isFibre(pCD) && pCD->custom_data.find_variable_index("fibre_custom_degradation") > 0){	
			if (pCD->custom_data["fibre_custom_degradation"] > 0.5)
				pCD->functions.instantiate_cell = instantiate_physimess_cell_custom_degrade;	
		}
	}
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	return; 
}
void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	std::cout << " x min :" << microenvironment.mesh.bounding_box[0] << " y min :" << microenvironment.mesh.bounding_box[1] << " z min :" << microenvironment.mesh.bounding_box[2] << std::endl;
	std::cout << " x max :" << microenvironment.mesh.bounding_box[3] << " y max :" << microenvironment.mesh.bounding_box[4] << " z max :" << microenvironment.mesh.bounding_box[5] << std::endl;
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin;
	/*Agents are set from a file or at random in number specified in the setting file*/
    if (load_cells_from_pugixml())
	{
		bool isFibreFromFile = true;
		// new fibre related parameters and bools
		isFibreFromFile = read_isFibreFromFile_status();
		
		std::cout << "Fibre from file: " << isFibreFromFile << std::endl;
		for( int i=0; i < (*all_cells).size(); i++ ) {
			
			if (isFibre((*all_cells)[i])) {
				
				if (!isFibreFromFile) {
					std::cerr << "Error: You have specified fibres position in your .csv to initialize agents position.\n"
					<< "If you wish to initialize fibre position from your .csv file modify IsFibreFromFile\n"
					<< "and set enable to true in your initial conditions in your xml setting file.\n"
					<< "Otherwise please choose a .csv position file without fibre data." << std::endl;
					exit(EXIT_FAILURE);
				}
				// Original fiber orientation assignment code
				static_cast<PhysiMeSS_Fibre*>((*all_cells)[i])->assign_fibre_orientation();
			}
		}

		/* Fibres agents are not added from the file but are created at random 
		given a relative volume index between 0 and 1 */

		if(!isFibreFromFile) {
			int my_type = read_FibreID(); 
			Cell_Definition* pCD = find_cell_definition( my_type );
			if( pCD != NULL )
			{
				
				double relative_fibres_volume = read_RelativeFibreVolume();
				double fibre_volume = pCD->custom_data["fibre_length"]*pCD->custom_data["fibre_radius"]*pCD->custom_data["fibre_radius"]*3.14;
				double microenvironment_volume = Xrange*Yrange*Zrange;
				int fibres_count = 0;
				while (fibres_count*fibre_volume/microenvironment_volume < relative_fibres_volume ) {
					Cell* pC;
					std::vector<double> position = {0, 0, 0};

					position[0] = Xmin + UniformRandom() * Xrange;
					position[1] = Ymin + UniformRandom() * Yrange;
					position[2] = Zmin + UniformRandom() * Zrange;

					pC = create_cell(*pCD);

					static_cast<PhysiMeSS_Fibre*>(pC)->assign_fibre_orientation();
					static_cast<PhysiMeSS_Fibre*>(pC)->check_out_of_bounds(position);

					pC->assign_position(position);
					fibres_count++;
				}
				
			}
	
		}
	}
	else {
		Cell* pC;
        std::vector<double> position = {0, 0, 0};

        for( int k=0; k < cell_definitions_by_index.size() ; k++ ) {

            Cell_Definition *pCD = cell_definitions_by_index[k];
            // std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
            
            if (!isFibre(pCD))
            {
                for (int n = 0; n < parameters.ints("number_of_cells"); n++) {

                    position[0] = Xmin + UniformRandom() * Xrange;
                    position[1] = Ymin + UniformRandom() * Yrange;
                    position[2] = Zmin + UniformRandom() * Zrange;

                    pC = create_cell(*pCD);
                                        
                    pC->assign_position(position);
                }
            } 
            
            else 
            {
                for ( int nf = 0 ; nf < parameters.ints("number_of_fibres") ; nf++ ) {

                    position[0] = Xmin + UniformRandom() * Xrange;
                    position[1] = Ymin + UniformRandom() * Yrange;
                    position[2] = Zmin + UniformRandom() * Zrange;

                    pC = create_cell(*pCD);

                    static_cast<PhysiMeSS_Fibre*>(pC)->assign_fibre_orientation();
                    static_cast<PhysiMeSS_Fibre*>(pC)->check_out_of_bounds(position);

                    pC->assign_position(position);
                }
            }
		}
	}
	remove_physimess_out_of_bounds_fibres();
	
	//adding a part to remove fibres that overlap with cells initialize at time zero
	for (auto cell: *all_cells)
	{
		if (!isFibre(cell))
		{
			for (auto* neighbor: cell->get_container()->agent_grid[cell->get_current_mechanics_voxel_index()])
        	{
				if(isFibre(neighbor))
				{
					std::vector<double> disp(3,0.0);
					if (norm(static_cast<PhysiMeSS_Fibre*>(neighbor)->nearest_point_on_fibre(cell->position, disp)) < (*cell).phenotype.geometry.radius + static_cast<PhysiMeSS_Fibre*>(neighbor)->mRadius)
					{
						delete_cell(neighbor);
					}
				}
			}
			for (auto neighbor_voxel_index: cell->get_container()->underlying_mesh.moore_connected_voxel_indices[cell->get_current_mechanics_voxel_index()])
        	{
				if(!is_neighbor_voxel(cell, cell->get_container()->underlying_mesh.voxels[cell->get_current_mechanics_voxel_index()].center, cell->get_container()->underlying_mesh.voxels[neighbor_voxel_index].center, neighbor_voxel_index))
					continue;
			
				for (auto* neighbor: cell->get_container()->agent_grid[neighbor_voxel_index])
				{
					if(isFibre(neighbor))
					{
						
						std::vector<double> disp(3,0.0);
						if (norm(static_cast<PhysiMeSS_Fibre*>(neighbor)->nearest_point_on_fibre(cell->position, disp)) < (*cell).phenotype.geometry.radius + static_cast<PhysiMeSS_Fibre*>(neighbor)->mRadius)
						{
							delete_cell(neighbor);
						}
						
					}
				}
			}
		}
	}
}

std::vector<std::string> paint_by_cell_pressure( Cell* pCell ){

	std::vector< std::string > output( 0);
	int color = (int) round( ((pCell->state.simple_pressure) / 10) * 255 );
	if(color > 255){
		color = 255;
	}
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,0,%u)", color, 255 - color);
	output.push_back( std::string("black") );
	output.push_back( szTempString );
	output.push_back( szTempString );
	output.push_back( szTempString );
	return output;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ 
	if (parameters.bools("color_cells_by_pressure")){
		return paint_by_cell_pressure(pCell); 
	} else {
		return paint_by_number_cell_coloring(pCell);
	}
}
std::string my_coloring_function_for_substrate( double concentration, double max_conc, double min_conc )
{ return paint_by_density_percentage( concentration,  max_conc,  min_conc); }

void my_cellcount_function(char* string){
	int nb_fibres = 0;
	for (Cell* cell : *all_cells) {
		if (isFibre(cell))
			nb_fibres++;
	}

	sprintf( string , "%lu cells, %u fibres" , all_cells->size()-nb_fibres, nb_fibres ); 
}

void custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

Cell* instantiate_physimess_cell() { return new PhysiMeSS_Cell; }
Cell* instantiate_physimess_fibre() { return new PhysiMeSS_Fibre; }
Cell* instantiate_physimess_cell_custom_degrade() { return new PhysiMeSS_Cell_Custom_Degrade; }


void PhysiMeSS_Cell_Custom_Degrade::degrade_fibre(PhysiMeSS_Fibre* pFibre)
{
	// Here this version of the degrade function takes cell pressure into account in the degradation rate
    double distance = 0.0;
    pFibre->nearest_point_on_fibre(position, displacement);
    for (int index = 0; index < 3; index++) {
        distance += displacement[index] * displacement[index];
    }
    distance = std::max(sqrt(distance), 0.00001);
    
    
        // Fibre degradation by cell - switched on by flag fibre_degradation
        double stuck_threshold = this->custom_data["fibre_stuck_time"];
        double pressure_threshold = this->custom_data["fibre_pressure_threshold"];
        if (this->custom_data["fibre_degradation"] > 0.5 && (stuck_counter >= stuck_threshold
                                                        || state.simple_pressure > pressure_threshold)) {
            // if (stuck_counter >= stuck_threshold){
            //     std::cout << "Cell " << ID << " is stuck at time " << PhysiCell::PhysiCell_globals.current_time
            //                 << " near fibre " << pFibre->ID  << std::endl;;
            // }
            // if (state.simple_pressure > pressure_threshold){
            //     std::cout << "Cell " << ID << " is under pressure of " << state.simple_pressure << " at "
            //                 << PhysiCell::PhysiCell_globals.current_time << " near fibre " << pFibre->ID  << std::endl;;
            // }
            displacement *= -1.0/distance;
            double dotproduct = dot_product(displacement, phenotype.motility.motility_vector);
            if (dotproduct >= 0) {
                double rand_degradation = PhysiCell::UniformRandom();
                double prob_degradation = this->custom_data["fibre_degradation_rate"];
                if (state.simple_pressure > pressure_threshold){
                    prob_degradation *= state.simple_pressure;
                }
                if (rand_degradation <= prob_degradation) {
                    //std::cout << " --------> fibre " << (*other_agent).ID << " is flagged for degradation " << std::endl;
                    // (*other_agent).parameters.degradation_flag = true;
                    pFibre->flag_for_removal();
                    // std::cout << "Degrading fibre agent " << pFibre->ID << " using flag for removal !!" << std::endl;
                    stuck_counter = 0;
                }
            }
        }
    // }
}


bool read_isFibreFromFile_status(pugi::xml_node config_root) {
	pugi::xml_node node;
	// Assign values to node and enable_vtk_saves inside a function
	node = xml_find_node(config_root, "initial_conditions");
	node = xml_find_node(node, "fibres_from_file");
	if (node) {
		return xml_get_bool_value(node, "enable");
	}
	else{
		return true; // by default fibres are thought to be defined in the file
	}
 }
 
 bool read_isFibreFromFile_status(void) {
	return read_isFibreFromFile_status(physicell_config_root);
 }


int read_FibreID(pugi::xml_node config_root) {
	pugi::xml_node node;
	// Assign values to node and enable_vtk_saves inside a function
	node = xml_find_node(config_root, "initial_conditions");
	node = xml_find_node(node, "fibres_from_file");
	if (!node.child("ID")) {
        return 1;
    }
	return xml_get_int_value(node, "ID");
 }
 
 int read_FibreID(void) {
	return read_FibreID(physicell_config_root);
 }

 double read_RelativeFibreVolume(pugi::xml_node config_root) {
	pugi::xml_node node;

	// Assign values to node and enable_vtk_saves inside a function
	node = xml_find_node(config_root, "initial_conditions");
	node = xml_find_node(node, "fibres_from_file");
	if (!node.child("relative_fibre_volume")) {
        return 0.0;
    }
	return xml_get_double_value(node, "relative_fibre_volume");;
 }
 
 double read_RelativeFibreVolume(void) {
	return read_RelativeFibreVolume(physicell_config_root);
 }