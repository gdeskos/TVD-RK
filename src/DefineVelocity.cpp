#include "AmrGVOF.H"
#include "Kernels.H"

#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void
AmrGVOF::DefineVelocityAllLevels (Real time)
{
    for (int lev = 0; lev <= finest_level; ++lev)
        DefineVelocityAtLevel(lev,time);

    // =======================================================
    // Average down face velocities before using them
    // =======================================================
    for (int lev = finest_level; lev > 0; lev--)
    {
        average_down_faces(amrex::GetArrOfConstPtrs(facevel[lev  ]),
                           amrex::GetArrOfPtrs     (facevel[lev-1]),
                           MaxRefRatio(lev-1), 0);
    }
}

void
AmrGVOF::DefineVelocityAtLevel (int lev, Real time)
{
    const auto dx = geom[lev].CellSizeArray();
    const Real* prob_lo = geom[lev].ProbLo();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(phi_new[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {

        // ======== GET FACE VELOCITY =========
            GpuArray<Box, AMREX_SPACEDIM> nbx;
            AMREX_D_TERM(nbx[0] = mfi.nodaltilebox(0);,
                         nbx[1] = mfi.nodaltilebox(1);,
                         nbx[2] = mfi.nodaltilebox(2););

            AMREX_D_TERM(const Box& ngbxx = amrex::grow(mfi.nodaltilebox(0),1);,
                         const Box& ngbxy = amrex::grow(mfi.nodaltilebox(1),1);,
                         const Box& ngbxz = amrex::grow(mfi.nodaltilebox(2),1););

            GpuArray<Array4<Real>, AMREX_SPACEDIM> vel{ AMREX_D_DECL( facevel[lev][0].array(mfi),
                                                                      facevel[lev][1].array(mfi),
                                                                      facevel[lev][2].array(mfi)) };

            auto prob_lo = geom[lev].ProbLoArray();
            auto dx = geom[lev].CellSizeArray();

            AMREX_D_TERM(
                         amrex::ParallelFor(ngbxx,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_x(i, j, k, vel[0], prob_lo, dx, time, problem_type); 
                         });,

                         amrex::ParallelFor(ngbxy,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_y(i, j, k, vel[1], prob_lo, dx, time, problem_type);
                         });,

                         amrex::ParallelFor(ngbxz,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_z(i, j, k, vel[2], prob_lo, dx, time, problem_type);
                         });
                        );
        }
    }
}
