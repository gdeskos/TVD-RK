#include "AmrGVOF.H"
#include "Kernels.H"

#include <AMReX_MultiFabUtil.H>

using namespace amrex;

// advance all levels for a single time step
void
AmrGVOF::LagrangianSplitAdvectionAllLevels (Real time, Real dt_lev, int iteration)
{
    constexpr int num_grow = 3;

    // Define the volume fraction fluxes (this is NOT actual fluxes a.k.a u_mac*phi*A)
    Vector< Array<MultiFab,3> > fluxes(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(idim);
            fluxes[lev][idim] = MultiFab(ba, dmap[lev], 1, 0);
        }
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
        std::swap(phi_old[lev], phi_new[lev]);
        t_old[lev] = t_new[lev];
        t_new[lev] += dt_lev;

        const Real old_time = t_old[lev];
        const Real new_time = t_new[lev];
        const Real ctr_time = 0.5*(old_time+new_time);

        const auto dx = geom[lev].CellSizeArray();
        GpuArray<Real, AMREX_SPACEDIM> dtdx;
        for (int i=0; i<AMREX_SPACEDIM; ++i)
            dtdx[i] = dt_lev/(dx[i]);
    
        const Real* prob_lo = geom[lev].ProbLo();

        // State with ghost cells
        MultiFab Sborder(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
        FillPatch(lev, time, Sborder, 0, Sborder.nComp());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {


            for (MFIter mfi(phi_new[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                GpuArray<Array4<Real>, AMREX_SPACEDIM> vel{ AMREX_D_DECL( facevel[lev][0].array(mfi),
                                                                          facevel[lev][1].array(mfi),
                                                                          facevel[lev][2].array(mfi)) };
                const Box& bx = mfi.tilebox();
                const Box& gbx = amrex::grow(bx, 1);

                Array4<Real> statein  = Sborder.array(mfi);
                Array4<Real> stateout = phi_new[lev].array(mfi);

                GpuArray<Array4<Real>, 3> flux{ AMREX_D_DECL(fluxes[lev][0].array(mfi),
                                                             fluxes[lev][1].array(mfi),
                                                             fluxes[lev][2].array(mfi)) };

                // Define Volume Fraction and Volume Fraction BCs
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                // Sweeping over the x-direction
                lagrangian_advection(i,j,k,0,dtdx[0],flux[0],vel[0],stateout);
                });
            }
        } //end omp
    } // end lev
}