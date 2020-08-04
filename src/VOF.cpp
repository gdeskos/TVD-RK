#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include "VOF.H"

using namespace amrex;

LevelSet::LevelSet ()
{
    ReadParameters();

}

// read in some parameters from inputs fil
void
LevelSet::ReadParameters()
{
    {
        ParmParse pp; //Traditionally, max_strep and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time",stop_time);
    }

    {
        ParmParse pp("amr"); // Added the amr prefix

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        pp.query("restart",restart_chkfile);
    }

    {
        ParmParse pp("levelset"); // Added the levelset prefix

        pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
    }
}

