#ifndef __PhysiMeSS_cell_h__
#define __PhysiMeSS_cell_h__

#include "PhysiMeSS_agent.h"
#include "PhysiMeSS_fibre.h"

class PhysiMeSS_Fibre;

class PhysiMeSS_Cell : public PhysiMeSS_Agent 
{
    private:
    public:
    
    int stuck_counter;
    int unstuck_counter;

    // strongest damping by the contact with a fibre during the current mechanics step, and
    // orientation of that fibre: applied once to the motility vector, in physimess_update_cell_velocity
    double fibre_contact_damping;
    std::vector<double> fibre_contact_orientation;
    
    PhysiMeSS_Cell();
    ~PhysiMeSS_Cell() {};
    
    virtual void add_potentials_from_fibre(PhysiMeSS_Fibre* fibre);
    
    void register_fibre_voxels();
    void deregister_fibre_voxels();
    
    void force_update_motility_vector(double dt_);
    
    virtual void degrade_fibre(PhysiMeSS_Fibre* pFibre); 

};

#endif