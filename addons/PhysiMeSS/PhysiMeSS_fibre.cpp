#include "PhysiMeSS_fibre.h"
#include "PhysiMeSS_cell.h"
#include "PhysiMeSS_geometry.h"
#include <algorithm>
#include <cmath>

// agents whose name contains one of these strings are fibres
static bool isFibreName(const std::string& agentname)
{
    static const std::vector<std::string> fibre_keywords = {"ecm", "matrix", "fiber", "fibre", "rod"};

    for (const auto& keyword : fibre_keywords) {
        if (agentname.find(keyword) != std::string::npos) {
            return true;
        }
    }
    return false;
}

bool isFibre(PhysiCell::Cell* pCell)
{
    return isFibreName(pCell->type_name);
}

bool isFibre(PhysiCell::Cell_Definition * cellDef)
{
    return isFibreName(cellDef->name);
}

std::vector<PhysiCell::Cell_Definition*>* getFibreCellDefinitions() {
    std::vector<PhysiCell::Cell_Definition*>* result = new std::vector<PhysiCell::Cell_Definition*>();

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
    fibres_crosslink_point.assign(3, 0.0);
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

/*
    Effect of a cell in contact with this fibre: the cell can rotate the fibre
    (about its centre, or about its crosslink point if it has exactly one) and
    push it. Fibres with two or more crosslinks do not move.

    In 2D simulations everything is computed in the z=0 plane.
*/
void PhysiMeSS_Fibre::add_potentials_from_cell(PhysiMeSS_Cell* cell)
{
    if ( !(this->phenotype.motility.is_motile) || this->X_crosslink_count >= 2) {
        return;
    }

    const bool planar = default_microenvironment_options.simulate_2D;
    const int dimensions = planar ? 2 : 3;

    // displacement: vector from the nearest point on the fibre axis to the cell centre
    nearest_point_on_fibre(cell->position, this->displacement);
    if (planar) {
        this->displacement[2] = 0.0;
    }

    double distance = std::max(norm(this->displacement), 0.00001);
    double R = cell->phenotype.geometry.radius + this->mRadius;
    if (distance > R) {
        return; // no interaction
    }

    if (cell->custom_data["fibre_rotation"] > 0.5) {
        // the cell acts on the fibre at the point of the fibre axis nearest to it
        std::vector<double> point_of_impact(3, 0.0);
        // a free fibre rotates about its centre, a fibre with one crosslink about the crosslink
        std::vector<double> pivot(3, 0.0);
        const std::vector<double>& pivot_position = (X_crosslink_count == 1) ? this->fibres_crosslink_point : this->position;
        for (int i = 0; i < dimensions; i++) {
            point_of_impact[i] = cell->position[i] - this->displacement[i];
            pivot[i] = pivot_position[i];
        }

        PhysiMeSS_Geometry::FibreRotation rotation = PhysiMeSS_Geometry::rotate_fibre_about_pivot(
            pivot, this->position, point_of_impact,
            cell->phenotype.motility.motility_vector,
            this->mLength, this->custom_data["fibre_sticky"]);

        if (!rotation.rotated) {
            return;
        }

        this->state.orientation = rotation.orientation;
        for (int i = 0; i < dimensions; i++) {
            this->position[i] = rotation.centre[i];
        }
    }

    // only free fibres, or pairs of fibres crosslinked to each other alone, can be pushed
    if ( X_crosslink_count == 0 ||
        ((this->fibres_crosslinkers).size() == 1 && static_cast<PhysiMeSS_Fibre*>((this->fibres_crosslinkers[0]))->fibres_crosslinkers.size() < 2)) {
        if (cell->custom_data["fibre_pushing"] > 0.5) {
            // repulsion as per PhysiCell
            static double simple_pressure_scale = 0.027288820670331;
            double temp_r = 1 - distance/R;
            temp_r *= temp_r; // (1 - d/R)^2

            // add the relative pressure contribution
            state.simple_pressure += (temp_r / simple_pressure_scale);

            double effective_repulsion = sqrt(this->phenotype.mechanics.cell_cell_repulsion_strength *
                                              cell->phenotype.mechanics.cell_cell_repulsion_strength);
            temp_r *= effective_repulsion;
            if (fabs(temp_r) > 1e-16) { // no push when there is almost no overlap between fibre and cell
                temp_r /= distance;
                naxpy(&velocity, temp_r, displacement); // velocity -= temp_r * displacement
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
    Computes the vector from the point of the fibre axis nearest to "point", to "point".
    The result is stored in "displacement" and returned.
*/
std::vector<double> PhysiMeSS_Fibre::nearest_point_on_fibre(const std::vector<double>& point, std::vector<double> &displacement)
{
    // don't bother if the "fibre_agent" is not a fibre
    if (!isFibre(this)) { return displacement; }

    displacement = PhysiMeSS_Geometry::displacement_from_fibre(point, this->position, this->state.orientation, this->mLength);
    return displacement;
}

/*
    Tests whether this fibre is crosslinked with a neighbouring fibre, i.e. whether
    the two fibres touch. If so, the neighbour is added to the crosslinkers of this
    fibre and the crosslink point (on this fibre, or the touching endpoint) is recorded.
*/
void PhysiMeSS_Fibre::check_fibre_crosslinks(PhysiMeSS_Fibre *fibre_neighbor)
{
    if (this == fibre_neighbor) { return; }
    if (!isFibre(this) || !isFibre(fibre_neighbor)) { return; }

    double co_length = this->mLength + fibre_neighbor->mLength;
    double co_radius = this->mRadius + fibre_neighbor->mRadius;

    // vector between fibre centres
    std::vector<double> centre_to_centre(3, 0.0);
    for (int i = 0; i < 3; i++) {
        centre_to_centre[i] = fibre_neighbor->position[i] - this->position[i];
    }
    double distance = norm(centre_to_centre);

    // the fibres are too far apart to touch
    if (distance > co_length) { return; }

    std::vector<double> p1(3, 0.0); // endpoints of "this" fibre
    std::vector<double> p2(3, 0.0);
    std::vector<double> p3(3, 0.0); // endpoints of "neighbor" fibre
    std::vector<double> p4(3, 0.0);
    for (int i = 0; i < 3; i++) {
        p1[i] = this->position[i] - mLength * this->state.orientation[i];
        p2[i] = this->position[i] + mLength * this->state.orientation[i];
        p3[i] = fibre_neighbor->position[i] - fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
        p4[i] = fibre_neighbor->position[i] + fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
    }

    std::vector<double> p1_to_p2(3, 0.0); // "this" fibre vector
    std::vector<double> p3_to_p4(3, 0.0); // "neighbor" fibre vector
    std::vector<double> p1_to_p3(3, 0.0); // from "this" fibre to "neighbor"
    std::vector<double> p3_to_p1(3, 0.0);
    for (int i = 0; i < 3; i++) {
        p1_to_p2[i] = p2[i] - p1[i];
        p3_to_p4[i] = p4[i] - p3[i];
        p1_to_p3[i] = p3[i] - p1[i];
        p3_to_p1[i] = p1[i] - p3[i];
    }

    // distance from each endpoint to the other fibre
    std::vector<double> disp(3, 0.0);
    double test_point1 = norm(fibre_neighbor->nearest_point_on_fibre(p1, disp));
    double test_point2 = norm(fibre_neighbor->nearest_point_on_fibre(p2, disp));
    double test_point3 = norm(this->nearest_point_on_fibre(p3, disp));
    double test_point4 = norm(this->nearest_point_on_fibre(p4, disp));

    normalize(&centre_to_centre);

    std::vector<double> n = cross_product(p1_to_p2, p3_to_p4);
    double n_norm = norm(n);

    if (n_norm < co_radius) {
        // (1) the fibres are (nearly) parallel

        if (centre_to_centre == this->state.orientation || centre_to_centre == -1.0 * this->state.orientation) {
            // colinear fibres touching end to end (distance <= co_length was checked above):
            // the crosslink point is the endpoint of this fibre
            this->fibres_crosslinkers.push_back(fibre_neighbor);
            for (int i = 0; i < 3; i++) {
                this->fibres_crosslink_point[i] = this->position[i] + this->mLength * centre_to_centre[i];
            }
        }
        else if (test_point1 <= co_radius || test_point2 <= co_radius ||
                 test_point3 <= co_radius || test_point4 <= co_radius) {
            // parallel fibres lying on top of one another
            this->fibres_crosslinkers.push_back(fibre_neighbor);
            this->fibres_crosslink_point = (test_point1 < test_point2) ? p1 : p2;
        }
        return;
    }

    // (2) the fibres are not parallel: the distance between their supporting lines
    //     is |(P3 - P1) . n| / |n|  with n = (P2 - P1) x (P4 - P3)
    double line_distance = std::fabs(dot_product(p1_to_p3, n)) / n_norm;
    if (line_distance > co_radius) {
        return;
    }

    std::vector<double> d1(3, 0.0); // unit vectors along each fibre
    std::vector<double> d2(3, 0.0);
    double p1_to_p2_norm = norm(p1_to_p2);
    double p3_to_p4_norm = norm(p3_to_p4);
    for (int i = 0; i < 3; i++) {
        d1[i] = p1_to_p2[i] / p1_to_p2_norm;
        d2[i] = p3_to_p4[i] / p3_to_p4_norm;
    }

    // closest points of the two lines: c1 = p1 + t1 d1 on "this", c2 = p3 + t2 d2 on "neighbor"
    std::vector<double> n2 = cross_product(d2, n);
    std::vector<double> n1 = cross_product(d1, n);
    double t1 = dot_product(p1_to_p3, n2) / dot_product(d1, n2);
    double t2 = dot_product(p3_to_p1, n1) / dot_product(d2, n1);
    std::vector<double> c1(3, 0.0);
    std::vector<double> c2(3, 0.0);
    std::vector<double> diff(3, 0.0);
    for (int i = 0; i < 3; i++) {
        c1[i] = p1[i] + t1 * d1[i];
        c2[i] = p3[i] + t2 * d2[i];
        diff[i] = c2[i] - c1[i];
    }

    if (norm(diff) > co_radius) {
        return;
    }

    if (0 <= t1 && t1 <= 2 * this->mLength && 0 <= t2 && t2 <= 2 * fibre_neighbor->mLength) {
        // the closest points of the lines are on both fibres
        this->fibres_crosslinkers.push_back(fibre_neighbor);
        this->fibres_crosslink_point = c1;
    }
    else if (test_point1 < co_radius || test_point2 < co_radius) {
        // an endpoint of this fibre touches the neighbor
        this->fibres_crosslinkers.push_back(fibre_neighbor);
        this->fibres_crosslink_point = (test_point1 < test_point2) ? p1 : p2;
    }
    else if (test_point3 < co_radius || test_point4 < co_radius) {
        // an endpoint of the neighbor touches this fibre
        this->fibres_crosslinkers.push_back(fibre_neighbor);
        this->fibres_crosslink_point = (test_point3 < test_point4) ? p3 : p4;
    }
}

void PhysiMeSS_Fibre::add_crosslinks()
{
    for (auto* neighbor : physimess_neighbors) {
        this->check_fibre_crosslinks(static_cast<PhysiMeSS_Fibre*>(neighbor));
    }
}
