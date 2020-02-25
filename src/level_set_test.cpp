#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include<AMReX_Geometry.H>
#include<AMReX_MultiFab.H>
#include<AMReX_Array.H>

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
