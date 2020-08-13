#ifndef LAGRANGIAN_ADVECTION_K_H_
#define LAGRANGIAN_ADVECTION_K_H_


#include "VolumeFractions_K.H"
#include <AMReX_FArrayBox.H>
#include <cmath>

/** Split advection of the interfae along the x,y,z directions
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void lagrangian_advection(
    int i,
    int j,
    int k,
    int direction,
    amrex::Real dtdx,
    amrex::Array4<amrex::Real> const& flux,
    amrex::Array4<amrex::Real> const& vel,
    amrex::Array4<amrex::Real> const& cc) noexcept
{
    amrex::Real qL, qR;
    amrex::Real const_tiny = std::numeric_limits<amrex::Real>::epsilon();

    qL = vel(i, j, k) * dtdx;
    qR = vel(i + 1, j, k) * dtdx;

    // 3 cases: 1: DEFAULT (c=0.0 and fluxes=0.0); 2: c=1.0; 3: c>0.0
    flux(i, j, k,0) = 0.0;
    flux(i, j, k,1) = 0.0;
    flux(i, j, k,2) = 0.0;

    if (cc(i, j, k) > 1.0 - const_tiny) {
        flux(i, j, k,0) = std::max(-qL, 0.0);
        flux(i, j, k,1) = 1.0 - std::max(qL, 0.0) + std::min(qR, 0.0);
        flux(i, j, k,2) = std::max(qR, 0.0);
    } else if (cc(i, j, k) > 0.0) {
        // (1) normal vector: mx,my,mz;
        //mixed_Youngs_centered<amr_wind::StencilInterior>(i, j, k, dx,dy,dz,cc, mxyz);
        amrex::Real mx, my, mz, invx, invy, invz, mmag;
        compute_normals(i,j,k,cc,mx,my,mz);

        // (2) mx,my,mz>0. and mx+my+mz = 1.0;
        invx = 1;
        invy = 1;
        invz = 1;
        if (mx < 0.0) {
            mx = -mx;
            invx = -1;
        }
        if (my < 0.0) {
            my = -my;
            invy = -1;
        }
        if (mz < 0.0) {
            mz = -mz;
            invz = -1;
        }
        amrex::Real mm2 = mx + my + mz;
        mx = mx / mm2;
        my = my / mm2;
        mz = mz / mm2;
        
        // (3) Get alpha;
        amrex::Real alpha = AL3D(i, j, k, mx, my, mz, cc);
        // (4) Back to original plane
        mx = invx * mx;
        my = invy * my;
        mz = invz * mz;
        alpha = alpha + std::min(0.0, mx) 
                      + std::min(0.0, my) 
                      + std::min(0.0, mz);

        // (5) Lagrangian Advection
        amrex::Real mm1 = std::max(qL, 0.0);
        mm2 = 1.0 - mm1 + std::min(0.0, qR);
        
        if (direction==0){
            mx = mx / (1.0 - qL + qR);
            alpha = alpha + mx * qL;    
            if (qL < 0.0) {
                flux(i, j, k,0) = FL3D(i, j, k, mx, my, mz, alpha, qL, -qL);
            }
            if (qR > 0.0) {
                flux(i, j, k, 2) = FL3D(i, j, k, mx, my, mz, alpha, 1.0, qR);
                flux(i, j, k, 1) = FL3D(i, j, k, mx, my, mz, alpha, mm1, mm2);
            }
        } else if (direction == 1) {
            my = my / (1.0 - qL + qR);
            alpha = alpha + my * qL;
            if (qL < 0.0) {
                flux(i, j, k, 0) = FL3D(i, j, k, my, mz, mx, alpha, qL, -qL);
            }
            if (qR > 0.0) {
                flux(i, j, k, 2) = FL3D(i, j, k, my, mz, mx, alpha, 1.0, qR);
                flux(i, j, k, 1) = FL3D(i, j, k, my, mz, mx, alpha, mm1, mm2);
            }
        } else if (direction == 2) {
            mz = mz / (1.0 - qL + qR);
            alpha = alpha + mz * qL;
            if (qL < 0.0) {
                flux(i, j, k, 0) = FL3D(i, j, k, mz, mx, my, alpha, qL, -qL);
            }
            if (qR > 0.0) {
                flux(i, j, k, 2) = FL3D(i, j, k, mz, mx, my, alpha, 1.0, qR);
                flux(i, j, k, 1) = FL3D(i, j, k, mz, mx, my, alpha, mm1, mm2);
            }
        }  
    }
}
// This applies boundary condition to the volume fluxes
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void bc_fluxes(
    int i,
    int j,
    int k,
    int d,
    amrex::Array4<amrex::Real> const& flux) noexcept
{
    if(d==0){
        flux(0,j,k,2)=0.0;
    }else if (d==1){
        flux(i,0,k,2)=0.0; 
    }else if (d==2){
        flux(i,j,0,2)=0.0;
    }
}

// This applies boundary condition to the volume fluxes
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void apply_fluxes(
    int i,
    int j,
    int k,
    int d,
    amrex::Array4<amrex::Real> const& flux,
    amrex::Array4<amrex::Real> const& cc) noexcept
{

    // Compute new values of cc and then do clipping
    if (d == 0) {
        cc(i, j, k) =
            flux(i - 1, j, k, 2) + flux(i, j, k, 1) + flux(i + 1, j, k, 0);
    } else if (d == 1) {
        cc(i, j, k) =
            flux(i, j - 1, k, 2) + flux(i, j, k, 1) + flux(i, j + 1, k, 0);
    } else if (d == 2) {
        cc(i, j, k) =
            flux(i, j, k - 1, 2) + flux(i, j, k, 1) + flux(i, j, k + 1, 0);
    }
    cc(i, j, k) = std::max(0.0, std::min(1.0, cc(i, j, k)));

    // Apply BCs to the volume fractions cc
}

#endif // LAGRANGIAN_ADVECTION_K.H
