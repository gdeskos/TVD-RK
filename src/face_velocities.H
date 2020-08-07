#ifndef FACE_VELOCITY_H_
#define FACE_VELOCITY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_x(int i, int j, int k,
                         Array4<Real> const& vx,
                         GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                         GpuArray<Real, AMREX_SPACEDIM> dx,
                         const Real time,
                         int prob_type)
{   
    if(prob_type==1){
    Real x = dx[0]*i + prob_lo[0]; // At the face
    Real y = dx[1]*(0.5+j) + prob_lo[1];
    Real z = dx[2]*(0.5+k) + prob_lo[2];
    vx(i,j,k) = 2.0*pow(sin(M_PI*x),2)*sin(2.0*M_PI*y)*sin(2.0*M_PI*z)*cos(M_PI*time/3.0);
    }else if(prob_type==0){
        amrex::Abort("Problem type not available");
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_y(int i, int j, int k,
                         Array4<Real> const& vy,
                         GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                         GpuArray<Real, AMREX_SPACEDIM> dx,
                         const Real time,
                         int prob_type)
{
    if(prob_type==1){
    Real x = dx[0]*(0.5+i) + prob_lo[0]; // At the face
    Real y = dx[1]*j + prob_lo[1];
    Real z = dx[2]*(0.5+k) + prob_lo[2];
    vy(i,j,k) =  -pow(sin(M_PI*y),2)*sin(2.0*M_PI*z)*sin(2.0*M_PI*x)*cos(M_PI*time/3.0); 
        }else if(prob_type==0){
        amrex::Abort("Problem type not available");
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_z(int i, int j, int k,
                         Array4<Real> const& vz,
                         GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                         GpuArray<Real, AMREX_SPACEDIM> dx,
                         const Real time,
                         int prob_type)
{
    if(prob_type==1){
    Real x = dx[0]*(0.5+i) + prob_lo[0]; // At the face
    Real y = dx[1]*(0.5+j) + prob_lo[1];
    Real z = dx[2]*0.5 + prob_lo[2];
    vz(i,j,k) =  -pow(sin(M_PI*z),2)*sin(2.0*M_PI*x)*sin(2.0*M_PI*y)*cos(M_PI*time/3.0); 
    }else if(prob_type==0){
        amrex::Abort("Problem type not available");
    }
}

#endif
