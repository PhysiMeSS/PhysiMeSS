#include "PhysiMeSS_fibre.h"
#include "PhysiMeSS_cell.h"
#include <algorithm>

bool isFibre(PhysiCell::Cell* pCell) 
{
    const auto agentname = std::string(pCell->type_name);
    const auto ecm = std::string("ecm");
    const auto matrix = std::string("matrix");
    const auto fiber = std::string("fiber");
    const auto fibre = std::string("fibre");
    const auto rod = std::string("rod");

    return (agentname.find(ecm) != std::string::npos ||
        agentname.find(matrix) != std::string::npos ||
        agentname.find(fiber) != std::string::npos ||
        agentname.find(fibre) != std::string::npos ||
        agentname.find(rod) != std::string::npos
    );
}

bool isFibre(PhysiCell::Cell_Definition * cellDef)
{
    const auto agentname = std::string(cellDef->name);
    const auto ecm = std::string("ecm");
    const auto matrix = std::string("matrix");
    const auto fiber = std::string("fiber");
    const auto fibre = std::string("fibre");
    const auto rod = std::string("rod");

    return (agentname.find(ecm) != std::string::npos ||
        agentname.find(matrix) != std::string::npos ||
        agentname.find(fiber) != std::string::npos ||
        agentname.find(fibre) != std::string::npos ||
        agentname.find(rod) != std::string::npos
    );
}

std::vector<PhysiCell::Cell_Definition*>* getFibreCellDefinitions() {
    std::vector<PhysiCell::Cell_Definition*>* result = new std::vector<PhysiCell::Cell_Definition*>();
    PhysiCell::Cell_Definition* pCD;
    
    
    for (auto& cd_name: PhysiCell::cell_definitions_by_name) {
        if (isFibre(cd_name.second)) {
            result->push_back(cd_name.second);
        }
    }
    
    return result;
}

PhysiMeSS_Fibre::PhysiMeSS_Fibre() 
{
    fibres_crosslinkers.clear();
    fibres_crosslink_point.clear();
    X_crosslink_count = 0;
    fail_count = 0;
}

void PhysiMeSS_Fibre::assign_fibre_orientation() 
{ 
    mLength = PhysiCell::NormalRandom(this->custom_data["fibre_length"], this->custom_data["length_normdist_sd"]) / 2.0;//mLength initialization
    mRadius = this->custom_data["fibre_radius"];
    this->assign_orientation();
    if (default_microenvironment_options.simulate_2D) {
        if (this->custom_data["anisotropic_fibres"] > 0.5){
            double theta = PhysiCell::NormalRandom(this->custom_data["fibre_angle"], this->custom_data["angle_normdist_sd"]);
            this->state.orientation[0] = std::cos(theta);
            this->state.orientation[1] = std::sin(theta);
        }
        else{
            this->state.orientation = PhysiCell::UniformOnUnitCircle();
        }
        this->state.orientation[2] = 0.0;
    }
    else {
        if (this->custom_data["anisotropic_fibres"] > 0.5) {
            double theta = PhysiCell::NormalRandom(this->custom_data["fibre_angle"], this->custom_data["angle_normdist_sd"]);
            double phi = PhysiCell::NormalRandom(this->custom_data["fibre_angle_phi"], this->custom_data["angle_normdist_sd"]);

            this->state.orientation[0] = std::cos(theta)*std::sin(phi);
            this->state.orientation[1] = std::sin(theta)*std::sin(phi);
            this->state.orientation[2] = std::cos(phi);
            normalize(this->state.orientation);
        }
        else {
            this->state.orientation = PhysiCell::UniformOnUnitSphere();
        }
        
    }
    //###########################################//
    //   this bit a hack for PacMan and maze	 //
    //###########################################//
    if (this->type_name == "fibre_vertical") {
        this->state.orientation[0] = 0.0;
        this->state.orientation[1] = 1.0;
        this->state.orientation[2] = 0.0;
    }
    if (this->type_name == "fibre_horizontal") {
        this->state.orientation[0] = 1.0;
        this->state.orientation[1] = 0.0;
        this->state.orientation[2] = 0.0;
    }
    //###########################################// 
}

void PhysiMeSS_Fibre::check_out_of_bounds(std::vector<double>& position)
{
    double Xmin = BioFVM::get_default_microenvironment()->mesh.bounding_box[0]; 
	double Ymin = BioFVM::get_default_microenvironment()->mesh.bounding_box[1]; 
	double Zmin = BioFVM::get_default_microenvironment()->mesh.bounding_box[2]; 

	double Xmax = BioFVM::get_default_microenvironment()->mesh.bounding_box[3]; 
	double Ymax = BioFVM::get_default_microenvironment()->mesh.bounding_box[4]; 
	double Zmax = BioFVM::get_default_microenvironment()->mesh.bounding_box[5]; 
    
    // start and end points of a fibre are calculated from fibre center
    double xs = position[0] - this->mLength * this->state.orientation[0];
    double xe = position[0] + this->mLength * this->state.orientation[0];
    double ys = position[1] - this->mLength * this->state.orientation[1];
    double ye = position[1] + this->mLength * this->state.orientation[1];
    double zs = position[2] - this->mLength * this->state.orientation[2];
    double ze = position[2] + this->mLength * this->state.orientation[2];
   
    /* check whether a fibre end point leaves the domain and if so initialise fibre again
                assume user placed the centre of fibre within the domain so reinitialise orientation,
                break after 10 failures
     */
    
    if (default_microenvironment_options.simulate_2D) {
        zs = 0;
        ze = 0;
    }
    // Anisotropic fibres have similar orientation, therefore if one of the endpoint is out of bound, no fibres will be set at this position
    if (this->custom_data["anisotropic_fibres"]) {
        if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
            ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax ||
            zs < Zmin || ze > Zmax || ze < Zmin || zs > Zmax ) {
            fail_count = 10;
        }
        
    }
    //Otherwise, if fibre orientation is not constrained an other random orientation is set if out of bound
    else{
        while (fail_count < 10) {
            if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax ||
                zs < Zmin || ze > Zmax || ze < Zmin || zs > Zmax ) {
                    fail_count++;

                    this->state.orientation = PhysiCell::UniformOnUnitSphere();
                    
                    xs = position[0] - mLength * this->state.orientation[0];
                    xe = position[0] + mLength * this->state.orientation[0];
                    ys = position[1] - mLength * this->state.orientation[1];
                    ye = position[1] + mLength * this->state.orientation[1];
                    zs = position[2] - mLength * this->state.orientation[2];
                    ze = position[2] + mLength * this->state.orientation[2];

                    if (default_microenvironment_options.simulate_2D) {
                        this->state.orientation[2] = 0;
                        zs = 0;
                        ze = 0;
                        
                    }
            }
            else {
                break;
            }
        }
    }
}


void PhysiMeSS_Fibre::add_potentials_from_cell(PhysiMeSS_Cell* cell) {
    if ( !(this->phenotype.motility.is_motile) || this->X_crosslink_count >= 2) {
        return;
    }

    nearest_point_on_fibre(cell->position, this->displacement);
    if (default_microenvironment_options.simulate_2D) {
        add_potentials_from_cell_2D(cell);
    }

    else if (!default_microenvironment_options.simulate_2D) {
        add_potentials_from_cell_3D(cell);
    }
}



void PhysiMeSS_Fibre::add_potentials_from_cell_2D(PhysiMeSS_Cell* cell) {
    double distance = 0.0;
    this->displacement[2] = 0.0;
    distance = std::max(norm(this->displacement), 0.00001);
    double R = cell->phenotype.geometry.radius + this->mRadius;
    
    if (distance > R){
        return; // no interaction
    }
    
    else if (distance <= R) { 

        std::vector<double> point_of_impact(3, 0.0); // displacement vector is the vector from cell center to the nearest position on the fibre,
        // point of impact is the point on the fibre vector nearest to the cell
        std::transform((*cell).position.begin(), (*cell).position.begin() + 2,
                        this->displacement.begin(),
                        point_of_impact.begin(),
                        [](double pos, double disp) { return pos - disp; });
        
         //center of rotation
         std::vector<double> center(3,0.0); 
         if ( X_crosslink_count == 0 ) {
             std::copy(this->position.begin(),
                     this->position.begin() + 2,
                     center.begin()); 
         }
         else if ( X_crosslink_count == 1) {
             std::copy(this->fibres_crosslink_point.begin(),
                 this->fibres_crosslink_point.begin() + 2,
                 center.begin()); 
         }
         std::vector<double> center_to_impact_vector(3,0.0);
         std::transform(center.begin(), center.begin() + 2,
                         point_of_impact.begin(),
                         center_to_impact_vector.begin(),
                         [](double pos, double imp) { return imp - pos; });
         double norm_center_to_impact = sqrt(dot_product(center_to_impact_vector, center_to_impact_vector));
         
         std::vector<double> old_orientation(3,0.0);
         std::transform(center_to_impact_vector.begin(), center_to_impact_vector.begin() + 2,
                         old_orientation.begin(),
                         [norm_center_to_impact](double pos) { return pos/norm_center_to_impact; });
        
        if (cell->custom_data["fibre_rotation"] > 0.5) {
            
            if (norm_center_to_impact < 1e-10) { //1e-10
                return; 
            }
            
            std::vector<double> torque(3,0.0);
            torque = cross_product(center_to_impact_vector,(*cell).phenotype.motility.motility_vector);

            // moment_arm_magnitude = length of pivot if orientation is a unit vector
            double torque_magnitude = norm(torque);
            double fibre_length = 2 * mLength;

            double moment_of_inertia = fibre_length * fibre_length / 12; //Formula for a rod
            double angular_acceleration = custom_data["fibre_sticky"] * torque_magnitude / moment_of_inertia;
            double angle = angular_acceleration;   
            if( angular_acceleration < 1e-16 ) {
                // if angular acceleration is small no tilt
                return;
            }

            //torque normalization
            std::transform(torque.begin(), torque.end(),
                        torque.begin(),
                        [torque_magnitude](double tor) { return tor/torque_magnitude; });
            

            double c = std::cos(angle);
            double s = std::sin(angle);
            
            // (k x v)
            std::vector<double> v_cross_k(3,0.0);
            v_cross_k = cross_product(torque,old_orientation);

            // Because pivot dot v=0 if pivot is perpendicular
            // the full Rodrigues formula simplifies (https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula):
            // v_new = vc + (k x v)s

            std::vector<double> new_orientation(3,0.0);
            std::transform(old_orientation.begin(), old_orientation.end(),
                            v_cross_k.begin(),
                            new_orientation.begin(),
                            [c, s](double or_val, double v_val) { return or_val * c + v_val * s; });

            state.orientation = new_orientation;
            normalize(&state.orientation);

            //CASE WHERE THE CROSSLINK POINT DIFFERENT FROM CENTER IS THE CENTER OF ROTATION
            double distance_center_position = PhysiCell::dist(this->position, center);
            
            std::vector<double> center_to_position_vector(3,0.0);
            std::transform(this->position.begin(), this->position.begin() + 2,
                            center.begin(),
                            center_to_position_vector.begin(),
                            [](double p, double c) { return p - c; });
            
            if ( distance_center_position > 1e-10) {
                for (int i=0; i<2; i++ ) {
                    this-> position[i] += (-old_orientation[i] + state.orientation[i])*distance_center_position;
                }
            }
             
        }

        if ( X_crosslink_count == 0 || 
            ((this->fibres_crosslinkers).size() == 1 && static_cast<PhysiMeSS_Fibre*>((this->fibres_crosslinkers[0]))->fibres_crosslinkers.size() < 2)) {       
            if (cell->custom_data["fibre_pushing"] > 0.5) {
                // as per PhysiCell
                static double simple_pressure_scale = 0.027288820670331;
                double temp_r = 1 - distance/R;
               
                temp_r *= temp_r;
                // add the relative pressure contribution 
                state.simple_pressure += (temp_r / simple_pressure_scale);

                double effective_repulsion = sqrt(this->phenotype.mechanics.cell_cell_repulsion_strength *
                (*cell).phenotype.mechanics.cell_cell_repulsion_strength);
                temp_r *= effective_repulsion;
                if (fabs(temp_r) > 1e-16) {// if temp_r very small <=> almost no overlap between fibre and cell, the overlap must be "significant" for rotation to occur
                    temp_r /= distance; //temp_r = (1.0 - distance/R)²/distance 
                    naxpy(&velocity, temp_r, displacement); // Performs the operation: this->velocity = this->velocity - temp_r * displacement
                }
            }
        }
    }
}

/*
In 3D simulation, this functions modifies the velocity and orientation of a fibre in function of the cells in contact with it.
*/
void PhysiMeSS_Fibre::add_potentials_from_cell_3D(PhysiMeSS_Cell* cell) {
    double distance = 0.0;
    distance = std::max(norm(this->displacement), 0.00001);
    double R = cell->phenotype.geometry.radius + this->mRadius;
    if (distance > R){
        return; // no interaction
    }
    else if (distance <= R) { 
        std::vector<double> point_of_impact(3, 0.0); // displacement vector is the vector from cell center to the nearest position on the fibre,
        // point of impact is the point on the fibre vector nearest to the cell
        std::transform((*cell).position.begin(), (*cell).position.end(),
                        this->displacement.begin(),
                        point_of_impact.begin(),
                        [](double pos, double disp) { return pos - disp; });
        
         //center of rotation
         std::vector<double> center(3,0.0); 
         if ( X_crosslink_count == 0 ) {
             std::copy(this->position.begin(),
                     this->position.end(),
                     center.begin()); 
         }
         else if ( X_crosslink_count == 1) {
             std::copy(this->fibres_crosslink_point.begin(),
                 this->fibres_crosslink_point.end(),
                 center.begin()); 
         }
         
         std::vector<double> center_to_impact_vector(3,0.0);
         std::transform(center.begin(), center.end(),
                         point_of_impact.begin(),
                         center_to_impact_vector.begin(),
                         [](double pos, double imp) { return imp - pos; });
         double norm_center_to_impact = sqrt(dot_product(center_to_impact_vector, center_to_impact_vector));
         
         std::vector<double> old_orientation(3,0.0);
         std::transform(center_to_impact_vector.begin(), center_to_impact_vector.end(),
                         old_orientation.begin(),
                         [norm_center_to_impact](double pos) { return pos/norm_center_to_impact; });
        
        if (cell->custom_data["fibre_rotation"] > 0.5) {
            if (norm_center_to_impact < 1e-10) { //1e-10
                return; 
            }
            
            //torque = orientation x motility
            std::vector<double> torque(3,0.0);
            torque = cross_product(center_to_impact_vector,(*cell).phenotype.motility.motility_vector);

            // moment_arm_magnitude = length of pivot if orientation is a unit vector
            double torque_magnitude = sqrt(dot_product(torque,torque));
            double fibre_length = 2 * mLength;

            double moment_of_inertia = fibre_length * fibre_length / 12; //Formula for a rod
            double angular_acceleration = custom_data["fibre_sticky"] * torque_magnitude / moment_of_inertia;
            double angle = angular_acceleration;
                
            if( angular_acceleration < 1e-16 ) {
                // if angular acceleration is small no tilt
                return;
            }

            //torque normalization
            std::transform(torque.begin(), torque.end(),
                        torque.begin(),
                        [torque_magnitude](double tor) { return tor/torque_magnitude; });
            

            double c = std::cos(angle);
            double s = std::sin(angle);
            
            // (k x v)
            std::vector<double> v_cross_k(3,0.0);
            v_cross_k = cross_product(torque,old_orientation);

            // Because pivot dot v=0 if pivot is perpendicular
            // the full Rodrigues formula simplifies (https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula):
            // v_new = vc + (k x v)s

            std::vector<double> new_orientation(3,0.0);
            std::transform(old_orientation.begin(), old_orientation.end(),
                            v_cross_k.begin(),
                            new_orientation.begin(),
                            [c, s](double or_val, double v_val) { return or_val * c + v_val * s; });

            state.orientation = new_orientation;
            normalize(&state.orientation);

            //CASE WHERE THE CROSSLINK POINT DIFFERENT FROM CENTER IS THE CENTER OF ROTATION
            double distance_center_position = PhysiCell::dist(this->position, center);
            
            std::vector<double> center_to_position_vector(3,0.0);
            std::transform(this->position.begin(), this->position.end(),
                            center.begin(),
                            center_to_position_vector.begin(),
                            [](double p, double c) { return p - c; });
            
            if ( distance_center_position > 1e-10) {
                for (int i=0; i<3; i++ ) {
                    this-> position[i] += (-old_orientation[i] + state.orientation[i])*distance_center_position;
                }
            } 
        }

        if ( X_crosslink_count == 0 || 
            ((this->fibres_crosslinkers).size() == 1 && static_cast<PhysiMeSS_Fibre*>((this->fibres_crosslinkers[0]))->fibres_crosslinkers.size() < 2)) {       
            if (cell->custom_data["fibre_pushing"] > 0.5) {
                // as per PhysiCell
                static double simple_pressure_scale = 0.027288820670331;//what does this do precisely?
                double temp_r = 1 - distance/R;
               
                temp_r *= temp_r;
                // add the relative pressure contribution NOT SURE IF NEEDED
                state.simple_pressure += (temp_r / simple_pressure_scale);

                double effective_repulsion = sqrt(this->phenotype.mechanics.cell_cell_repulsion_strength *
                (*cell).phenotype.mechanics.cell_cell_repulsion_strength);
                temp_r *= effective_repulsion;
                if (fabs(temp_r) > 1e-16) {// if temp_r very small <=> almost no overlap between fibre and cell, the overlap must be "significant" for rotation to occur
                    temp_r /= distance; //temp_r = (1.0 - distance/R)²/distance   
                    naxpy(&velocity, temp_r, displacement); // Performs the operation: this->velocity = this->velocity - temp_r * displacement
                }
            }
        }
    }
}


void PhysiMeSS_Fibre::add_potentials_from_fibre(PhysiMeSS_Fibre* other_fibre) 
{
    /* probably want something here to model tension along fibres
        * this will be strong tension along the fibre for those fibres with a crosslink
        * and weak tension with background ECM */
    return;
}

void PhysiMeSS_Fibre::register_fibre_voxels() {

    int voxel;
    int voxel_size = this->get_container()->underlying_mesh.dx; // note this must be the same as the mechanics_voxel_size
    int test = 2.0 * this->mLength / voxel_size; //allows us to sample along the fibre

    std::vector<double> fibre_start(3, 0.0);
    std::vector<double> fibre_end(3, 0.0);
    for (unsigned int i = 0; i < 3; i++) {
        fibre_start[i] = this->position[i] - this->mLength * this->state.orientation[i];
        fibre_end[i] = this->position[i] + this->mLength * this->state.orientation[i];
    }
    // first add the voxel of the fibre end point
    voxel = this->get_container()->underlying_mesh.nearest_voxel_index(fibre_end);
    physimess_voxels.push_back(voxel);
    if (std::find(this->get_container()->agent_grid[voxel].begin(),
                    this->get_container()->agent_grid[voxel].end(),
                    this) == this->get_container()->agent_grid[voxel].end()) {
        this->get_container()->agent_grid[voxel].push_back(this);
    }
    // then walk along the fibre from fibre start point sampling and adding voxels as we go
    std::vector<double> point_on_fibre(3, 0.0);
    for (unsigned int j = 0; j < test + 1; j++) {
        for (unsigned int i = 0; i < 3; i++) {
            point_on_fibre[i] = fibre_start[i] + j * voxel_size * this->state.orientation[i];
        }
        voxel = this->get_container()->underlying_mesh.nearest_voxel_index(point_on_fibre);
        physimess_voxels.push_back(voxel);
        if (std::find(this->get_container()->agent_grid[voxel].begin(),
                        this->get_container()->agent_grid[voxel].end(),
                        this) == this->get_container()->agent_grid[voxel].end()) {
            this->get_container()->agent_grid[voxel].push_back(this);
        }
    }

    physimess_voxels.sort();
    physimess_voxels.unique();
}

void PhysiMeSS_Fibre::deregister_fibre_voxels() 
{
    int centre_voxel = this->get_container()->underlying_mesh.nearest_voxel_index(this->position);
    for (int voxel: physimess_voxels) {
        if (voxel != centre_voxel) {
            this->get_container()->remove_agent_from_voxel(this, voxel);
        }
    }
}


/*
This function modifies computes the point on the central axis of the fibre closest to the argument point
*/
std::vector<double> PhysiMeSS_Fibre::nearest_point_on_fibre(std::vector<double> point, std::vector<double> &displacement) 
{
    // don't bother if the "fibre_agent" is not a fibre
    if (!isFibre(this)) { return displacement; }

    double fibre_length = 2 * this->mLength;
    // vector pointing from one endpoint of "fibre_agent" to "point"
    std::vector<double> fibre_to_agent(3, 0.0);
    // scalar product fibre_to_agent * fibre_vector
    double fibre_to_agent_dot_endpoint_to_center = 0;
    std::vector<double> endpoint(3,0.0);
    double distance = 0;

    for (unsigned int i = 0; i < 3; i++) {
        endpoint[i] = this->position[i] - this->mLength * this->state.orientation[i];
        fibre_to_agent[i] = point[i] - endpoint[i];
        fibre_to_agent_dot_endpoint_to_center += fibre_to_agent[i] * 2 * (position[i] - endpoint[i]);//Magnitude of the vector fiber agent projected on the vector orientation of the fiber
    }
    // First 2 cases: one of the endpoint is the closest point on the fibre from the point
    // "point" is closest to the selected endpoint of this fibre
    if (fibre_to_agent_dot_endpoint_to_center < 0.) {
        for (int i = 0; i < 3; i++) {
            displacement[i] = fibre_to_agent[i];
        }
    }
        // “point” is closest to the other endpoint of this fibre
    else if (fibre_to_agent_dot_endpoint_to_center > fibre_length * fibre_length) {
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = point[i] - (this->position[i]
                                            + this->mLength * this->state.orientation[i]);
        }
    }
       
       // Final case: “point” is closest to a point along this fibre
    else {
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = fibre_to_agent[i] - (position[i] - endpoint[i])/mLength * fibre_to_agent_dot_endpoint_to_center/fibre_length;
        }
    }
    // the function returns the displacement vector
    return displacement;
}

void PhysiMeSS_Fibre::check_fibre_crosslinks(PhysiMeSS_Fibre *fibre_neighbor) {
    
    if (this == fibre_neighbor) { return; }
    else if (isFibre(this) && isFibre(fibre_neighbor)) {

            //computing distance between fibre and its neighbor
            double co_length = this->mLength + fibre_neighbor->mLength;
            std::vector<double> centre_to_centre(3, 0.0);
            
            for (int i = 0; i < 3; i++) {
                // vector between fibre centres
                centre_to_centre[i] = fibre_neighbor->position[i] - this->position[i];
            }

            if (norm(centre_to_centre) <= co_length) { //if it is possible for the fibre to touch, check if there is a crosslink
            std::vector<double> p1(3, 0.0);
            std::vector<double> p2(3, 0.0);
            std::vector<double> p3(3, 0.0);
            std::vector<double> p4(3, 0.0);
            for (int i = 0; i < 3; i++) {
                // endpoints of "this" fibre
                p1[i] = this->position[i] - mLength * this->state.orientation[i];
                p2[i] = this->position[i] + mLength * this->state.orientation[i];
                
                // endpoints of "neighbor" fibre
                p3[i] = fibre_neighbor->position[i] - fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
                p4[i] = fibre_neighbor->position[i] + fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
            }

            //vectors between fibre endpoints
            std::vector<double> p1_to_p2(3, 0.0);
            std::vector<double> p3_to_p4(3, 0.0);
            std::vector<double> p1_to_p3(3, 0.0);
            std::vector<double> p3_to_p1(3, 0.0);
            
            for (int i = 0; i < 3; i++) {
                // "fibre" fibre vector
                p1_to_p2[i] = p2[i] - p1[i];
                // "neighbor" fibre vector
                p3_to_p4[i] = p4[i] - p3[i];
                // vector from "fibre" to "neighbor"
                p1_to_p3[i] = p3[i] - p1[i];
                p3_to_p1[i] = p1[i] - p3[i];
                // vector between fibre centres
                //centre_to_centre[i] = fibre_neighbor->position[i] - this->position[i];
            }

            std::vector<double> disp1(3, 0.0);
            std::vector<double> disp2(3, 0.0);
            std::vector<double> disp3(3, 0.0);
            std::vector<double> disp4(3, 0.0);
            //std::cout << "This is running" << std::endl;
            //QUESTIONS: If fibres are perfectly parallel, 2 tests are enough, if of course we don't allow for tolerance i the colinearity
            fibre_neighbor->nearest_point_on_fibre(p1, disp1);
            double test_point1 = norm(disp1);
            fibre_neighbor->nearest_point_on_fibre(p2, disp2);
            double test_point2 = norm(disp2);
            this->nearest_point_on_fibre(p3, disp3);
            double test_point3 = norm(disp3);
            this->nearest_point_on_fibre(p4, disp4);
            double test_point4 = norm(disp4);

            double co_radius = this->mRadius + fibre_neighbor->mRadius;
            double co_length = this->mLength + fibre_neighbor->mLength;

            double distance = sqrt(dot_product(centre_to_centre, centre_to_centre)); //this computes the norm of centre_to_centre
              
            normalize(&centre_to_centre);

            std::vector<double> n = cross_product(p1_to_p2, p3_to_p4);
            double n_norm = sqrt(dot_product(n, n));
            
            if (n_norm < co_radius) {
                // CASE OF COLINEAR  FIBRE VECTORS 
                if (centre_to_centre == this->state.orientation || centre_to_centre == -1.0 * this->state.orientation) {
                    if (distance <= co_length) {
                        this->fibres_crosslinkers.push_back(fibre_neighbor);
                        std::vector<double> fibres_crosslink_point(3, 0.0);
                        for (int i=0; i<3;i++) {
                            this->fibres_crosslink_point[i] = this->position[i] + this->mLength * centre_to_centre[i]; //the chosen crosslink point is the endpoint of this fibre
                            //fibre_neighbor->fibres_crosslink_point[i] = fibre_neighbor->position[i] - fibre_neighbor->mLength * centre_to_centre[i];
                        }
                        return;
                    }
                    else { 
                        return; 
                    }
                }
                    

               else if ((test_point1 <= co_radius ||
                    test_point2 <= co_radius || test_point3 <= co_radius || test_point4 <= co_radius)) {
                    
                    this->fibres_crosslinkers.push_back(fibre_neighbor);
                    //fibre_neighbor->fibres_crosslinkers.push_back(this);
        
                    if (test_point1 < test_point2) {
                        this->fibres_crosslink_point = p1;
                        //fibre_neighbor->fibres_crosslink_point = p1; //QUESTIONS:
                        }
                    else {
                        this->fibres_crosslink_point = p2;
                        //fibre_neighbor->fibres_crosslink_point = p2;
                    }
                }

                else {return;}
            }
                    
            else {
                //CASE WHERE VECTORS ARE NOT COLINEAR
                // Check if lines are coplanar
                double triple_product = dot_product(p1_to_p3, n);
                // The distance between lines is |triple_product|/|cross|
                double line_distance = abs(triple_product) / n_norm;
               
                if (line_distance > co_radius) {
                    return; // Lines are further apart than co_radius
                }
                std::vector<double> d1 = {p1_to_p2[0]/norm(p1_to_p2), p1_to_p2[1]/norm(p1_to_p2), p1_to_p2[2]/norm(p1_to_p2)};
                std::vector<double> d2 = {p3_to_p4[0]/norm(p3_to_p4), p3_to_p4[1]/norm(p3_to_p4), p3_to_p4[2]/norm(p3_to_p4)};

                std::vector<double> n2 = cross_product(d2, n);
                std::vector<double> n1 = cross_product(d1, n);
                std::vector<double> c1(3, 0.0);
                std::vector<double> c2(3, 0.0);
                for (int i=0; i<3; i++) {
                    c1[i] = p1[i] + dot_product(p1_to_p3, n2)/dot_product(d1, n2) * d1[i];
                    c2[i] = p3[i] + dot_product(p3_to_p1, n1)/dot_product(d2, n1) * d2[i];

                }
                  
                // Verify the distance is within tolerance
                std::vector<double> diff = {c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2]};
                
                if (sqrt(dot_product(diff,diff)) <= co_radius) {
                    if (dot_product(p1_to_p3, n2)/dot_product(d1, n2) <= (2*this->mLength) && 0 <= dot_product(p1_to_p3, n2)/dot_product(d1, n2) && dot_product(p3_to_p1, n1)/dot_product(d2, n1) <= (2*fibre_neighbor->mLength) && 0 <= dot_product(p3_to_p1, n1)/dot_product(d2, n1)) {
                        this->fibres_crosslinkers.push_back(fibre_neighbor);
                        //fibre_neighbor->fibres_crosslinkers.push_back(this);
    
                        this->fibres_crosslink_point = c1;
                        //fibre_neighbor->fibres_crosslink_point = c2;
                        return;
                   }
                   
                    else if (test_point1 < co_radius || test_point2 < co_radius ) {
                        this->fibres_crosslinkers.push_back(fibre_neighbor);
                        //fibre_neighbor->fibres_crosslinkers.push_back(this);
                        if (test_point1 < test_point2) {
                            this->fibres_crosslink_point = p1;
                        }
                        else {
                            this->fibres_crosslink_point = p2;
                        }
                        
                        //fibre_neighbor->fibres_crosslink_point = c2;
                        return;

                    }
                    else if ( test_point3 < co_radius || test_point4 < co_radius ) {
                        this->fibres_crosslinkers.push_back(fibre_neighbor);
                        //fibre_neighbor->fibres_crosslinkers.push_back(this);
                        if (test_point3 < test_point4) {
                            this->fibres_crosslink_point = p3;
                        }
                        else {
                            this->fibres_crosslink_point = p4;
                        }
                        
                        //fibre_neighbor->fibres_crosslink_point = c2;
                        return;

                    }
                    else {return;}
                }
            }
        }
    }
            
     
    else { return; }
}

/*Establishes one crosslink point and the fibre crosslinkers as well as for it's neighbor
before check_fibre_crosslinks is run, crosslikers and crosslink point are cleared for ALL fibre 
therefore if there is a crosslink, we add crosslinker and crosslink point to it's neighbor as well
we don't run the code if the crosslink has already been established ! */
/*
/*void PhysiMeSS_Fibre::check_fibre_crosslinks(PhysiMeSS_Fibre *fibre_neighbor) {

    //std::cout << "check_fibre_crosslinks" << std::endl;
    if (this == fibre_neighbor) { return; }

    else if (isFibre(this) && isFibre(fibre_neighbor)) {
        //if the crosslink has aleady been established, we don't consider it 
        std::cout << "The timestep is : " << PhysiCell::PhysiCell_globals.current_time << " we consider the fibre " << this->ID <<  " and its neighbor " << fibre_neighbor->ID << std::endl;
        if (std::find(this->fibres_crosslinkers.begin(), 
        this->fibres_crosslinkers.end(),fibre_neighbor) == this->fibres_crosslinkers.end()) {
            for (int i = 0; i < this->fibres_crosslinkers.size(); i++)
            {std::cout << "crosslinkers : " << (this->fibres_crosslinkers[i])->ID << std::endl;}
            // fibre endpoints
            std::vector<double> p1(3, 0.0);
            std::vector<double> p2(3, 0.0);
            std::vector<double> p3(3, 0.0);
            std::vector<double> p4(3, 0.0);
            for (int i = 0; i < 3; i++) {
                // endpoints of "this" fibre
                p1[i] = this->position[i] - mLength * this->state.orientation[i];
                p2[i] = this->position[i] + mLength * this->state.orientation[i];
                
                // endpoints of "neighbor" fibre
                p3[i] = fibre_neighbor->position[i] - fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
                p4[i] = fibre_neighbor->position[i] + fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
            }

            //vectors between fibre endpoints
            std::vector<double> p1_to_p2(3, 0.0);
            std::vector<double> p3_to_p4(3, 0.0);
            std::vector<double> p1_to_p3(3, 0.0);
            std::vector<double> p3_to_p1(3, 0.0);
            std::vector<double> centre_to_centre(3, 0.0);
            for (int i = 0; i < 3; i++) {
                // "fibre" fibre vector
                p1_to_p2[i] = p2[i] - p1[i];
                // "neighbor" fibre vector
                p3_to_p4[i] = p4[i] - p3[i];
                // vector from "fibre" to "neighbor"
                p1_to_p3[i] = p3[i] - p1[i];
                p3_to_p1[i] = p1[i] - p3[i];
                // vector between fibre centres
                centre_to_centre[i] = fibre_neighbor->position[i] - this->position[i];
            }

            std::vector<double> disp1(3, 0.0);
            std::vector<double> disp2(3, 0.0);
            std::vector<double> disp3(3, 0.0);
            std::vector<double> disp4(3, 0.0);
            //QUESTIONS: If fibres are perfectly parallel, 2 tests are enough, if of course we don't allow for tolerance i the colinearity
            fibre_neighbor->nearest_point_on_fibre(p1, disp1);
            double test_point1 = norm(disp1);
            fibre_neighbor->nearest_point_on_fibre(p2, disp2);
            double test_point2 = norm(disp2);
            this->nearest_point_on_fibre(p3, disp3);
            double test_point3 = norm(disp3);
            this->nearest_point_on_fibre(p4, disp4);
            double test_point4 = norm(disp4);

            double co_radius = this->mRadius + fibre_neighbor->mRadius;
            double co_length = this->mLength + fibre_neighbor->mLength;

            double distance = sqrt(dot_product(centre_to_centre, centre_to_centre)); //this computes the norm of centre_to_centre
            normalize(this->state.orientation);
            normalize(fibre_neighbor->state.orientation);
            normalize(&centre_to_centre);

            /* test if fibres intersect
                (1) if fibres are coplanar and parallel:
                the cross product of the two fibre vectors is zero
                [(P2 - P1) x (P4 - P3)].[(P2 - P1) x (P4 - P3)] = 0 
            std::vector<double> n = cross_product(p1_to_p2, p3_to_p4);
            double n_norm = sqrt(dot_product(n, n));
            /* coplanar parallel fibres could intersect if colinear
                i.e. the orientation of the fibres are parallel or
                antiparallel to the centre_to_centre vector and
                distance between fibre centres is less than their co_length 
            if (n_norm < co_radius) {
                // CASE OF COLINEAR  FIBRE VECTORS 
                if (centre_to_centre == this->state.orientation || centre_to_centre == -1.0 * this->state.orientation) {
                    if (distance <= co_length) {
                        //std::cout << "This is happening " << std::endl;
                        this->fibres_crosslinkers.push_back(fibre_neighbor);
                        fibre_neighbor->fibres_crosslinkers.push_back(this);
                        for (int i=0; i<3;i++) {
                            this->fibres_crosslink_point[i] = this->position[i] + this->mLength * centre_to_centre[i]; //the chosen crosslink point is the endpoint of this fibre
                            fibre_neighbor->fibres_crosslink_point[i] = fibre_neighbor->position[i] - fibre_neighbor->mLength * centre_to_centre[i];
                        }
                        return;
                    }
                    else { 
                        return; 
                    }
                }
                    
                /* (2) parallel fibres may sit on top of one another
                    we check the distance between fibre end points and
                    the nearest point on neighbor fibre to see if they do 

               else if ((std::abs(test_point1) <= co_radius ||
                    std::abs(test_point2) <= co_radius) &&
                    centre_to_centre != this->state.orientation &&
                    centre_to_centre != -1.0 * this->state.orientation) {
                        
                    this->fibres_crosslinkers.push_back(fibre_neighbor);
                    fibre_neighbor->fibres_crosslinkers.push_back(this);
        
                    if (test_point1 < test_point2) {
                        this->fibres_crosslink_point = p1;
                        fibre_neighbor->fibres_crosslink_point = p1; //QUESTIONS:
                        }
                    else {
                        this->fibres_crosslink_point = p2;
                        fibre_neighbor->fibres_crosslink_point = p2;
                    }
                    /*
                    if (test_point3 < test_point4) {
                        fibre_neighbor->fibres_crosslink_point = p3; //QUESTIONS:
                    }
                    else {
                        fibre_neighbor->fibres_crosslink_point = p4;
                    }
                        
                    return;
                }

                else {return;}
            }
                    
            /*  (3) if fibres are skew (in parallel planes):
                the scalar triple product (P3 - P1) . [(P2 - P1) x (P4 - P3)] != 0
                so intersecting fibres require (P3 - P1) . [(P2 - P1) x (P4 - P3)] = 0
                we include a tolerance on this to allow for fibre radius 
            else {
                //CASE WHERE VECTORS ARE NOT COLINEAR
                // Check if lines are coplanar
                double triple_product = dot_product(p1_to_p3, n);
                // The distance between lines is |triple_product|/|cross|
                double line_distance = abs(triple_product) / n_norm;
    
                if (line_distance > co_radius) {
                    return; // Lines are further apart than co_radius
                }
                std::vector<double> d1 = {p1_to_p2[0]/norm(p1_to_p2), p1_to_p2[1]/norm(p1_to_p2), p1_to_p2[2]/norm(p1_to_p2)};
                std::vector<double> d2 = {p3_to_p4[0]/norm(p3_to_p4), p3_to_p4[1]/norm(p3_to_p4), p3_to_p4[2]/norm(p3_to_p4)};

                std::vector<double> n2 = cross_product(d2, n);
                std::vector<double> n1 = cross_product(d1, n);
                std::vector<double> c1(3, 0.0);
                std::vector<double> c2(3, 0.0);
                for (int i=0; i<3; i++) {
                    c1[i] = p1[i] + dot_product(p1_to_p3, n2)/dot_product(d1, n2) * d1[i];
                    c2[i] = p3[i] + dot_product(p3_to_p1, n1)/dot_product(d2, n1) * d2[i];

                }
                  
                // Verify the distance is within tolerance
                std::vector<double> diff = {c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2]};
                
                if (sqrt(dot_product(diff,diff)) <= co_radius) {
                    if (dot_product(p1_to_p3, n2)/dot_product(d1, n2) <= (2*this->mLength) && 0 <= dot_product(p1_to_p3, n2)/dot_product(d1, n2) && dot_product(p3_to_p1, n1)/dot_product(d2, n1) <= (2*fibre_neighbor->mLength) && 0 <= dot_product(p3_to_p1, n1)/dot_product(d2, n1)) {
                        this->fibres_crosslinkers.push_back(fibre_neighbor);
                        fibre_neighbor->fibres_crosslinkers.push_back(this);
    
                        this->fibres_crosslink_point = c1;
                        fibre_neighbor->fibres_crosslink_point = c2;
                        std::cout << "This is going on" << std::endl;
                        return;
                   }
                   
                    else if (test_point1 < co_radius || test_point2 < co_radius || test_point3 < co_radius || test_point3 < co_radius ) {
                        //std::cout << co_radius <<", " << test_point1 << ", " << test_point2 <<  ", " << test_point3 << ", " << test_point4 << std::endl;
                        this->fibres_crosslinkers.push_back(fibre_neighbor);
                        fibre_neighbor->fibres_crosslinkers.push_back(this);
                        //std::cout << "This fibre: " << this->ID << " has a crosslink with its neighbor: " << fibre_neighbor->ID << std::endl;
    
                        this->fibres_crosslink_point = c1;
                        fibre_neighbor->fibres_crosslink_point = c2;
                        //std::cout << "This is going on 2" << std::endl;
                        return;

                    }
                }
            }
        }
    } 
    else { return; }
}
*/


void PhysiMeSS_Fibre::add_crosslinks() 
{
    for (auto* neighbor : physimess_neighbors)
    {
           
            this->check_fibre_crosslinks(static_cast<PhysiMeSS_Fibre*>(neighbor));
    }
}