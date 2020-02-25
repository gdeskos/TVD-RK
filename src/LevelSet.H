#ifndef LevelSet_H_
#define LevelSet_H_

#include <string>
#include <limits>
#include <memory>

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_BCRec.H>

using namespace amrex;

class LevelSet
    : public amrex::AmrCore
{
public:
    ////////////////////////////////////////////////////////////////
    // public member function
    //
    // constructor - reads in parameters from inputs fil
    //             - sizes multilevel arrays and data structures
    ////////////////////////////////////////////////////////////////
    LevelSet ();

private:

    ////////////////////////////////////////////////////////////////
    // private member functions
    ////////////////////////////////////////////////////////////////
    
    // read in some parameters from inputs file
    void ReadParameters();

    ////////////////////////////////////////////////////////////////
    // private data members
    ////////////////////////////////////////////////////////////////
    amrex::Vector<int> istep; //which step ?
    amrex::Vector<int> nsubsteps; //how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    ////////////////////////////////////////////////////////////////
    // runtime parameters
    ////////////////////////////////////////////////////////////////
    
    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // if >=0 we start from a checkpoint
    std::string restart_chkfile = "";

    // advective cfl number - dt = cfl*dx/umax
    amrex::Real cfl = 0.7;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // hyperbolic refluxing as part of multilevel synchronization
    int do_reflux = 1;

    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    int plot_int = -1;

    // checkpoint prefix and frequency
    std::string chk_file {"chk"};
    int chk_int = -1;

};

#endif
