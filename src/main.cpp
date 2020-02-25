#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include "LevelSet.H"

using namespace amrex;

void initialize_fields();

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    initialize_fields();

    amrex::Finalize();
    return 0;
}

void initialize_fields()
{
    Real strt_time = amrex::second();
    amrex::Print() << "Initialising fields for level-set"<< std::endl;
    Real stop_time = amrex::second() - strt_time;
    
    amrex::Print() << "total time of initialization "<<stop_time<< std::endl;
}
