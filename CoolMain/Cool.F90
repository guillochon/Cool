!!****f* source/physics/sourceTerms/Cool/Cool
!!
!! NAME
!!
!!  Cool
!!
!! SYNOPSIS
!!
!!  Cool(integer(IN) :: blockCount
!!       integer(IN) :: blockList(blockCount),
!!          real(IN) :: dt,
!!          real(IN) :: time)
!!
!!
!!
!! DESCRIPTION
!!  Apply a cooling operator on the list of blocks provided as input
!!
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
!!  time : the current time
!!
!!***



subroutine Cool(blockCount,blockList,dt, time)

  use Grid_interface, ONLY : Grid_fillGuardCells, &
       Grid_getBlkIndexLimits, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped
  use Grid_data, ONLY : gr_smalle
  use Multispecies_interface, ONLY : Multispecies_getSumInv
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Simulation_data, ONLY : sim_smallT, sim_xCenter, sim_yCenter, sim_zCenter, &
                              sim_ptMass, sim_tAmbient, obj_mu, sim_condCoeff
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Cool_data, ONLY : cool_useCool

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"

  integer, intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real,intent(IN) :: dt, time

  integer :: thisBlock, blockID, i, j, k, sizeX, sizeY, sizeZ, istat
  logical :: cooledZone
  real :: sdot, sgamma, ek, ei, rho, mp, me, abar, xx, yy, zz, dist
  real :: rsc, T0, Tback, temp, kb, newton, softening_radius
  real, pointer, dimension(:,:,:,:)            :: solnData
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  double precision, dimension(SPECIES_BEGIN:SPECIES_END) :: xn
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  logical :: gcell = .true.

  if (.not. cool_useCool) return

  call PhysicalConstants_get("proton mass", mp)
  call PhysicalConstants_get("electron mass", me)
  call PhysicalConstants_get("Boltzmann", kb)
  call PhysicalConstants_get("Newton", newton)

  call RuntimeParameters_get("sink_softening_radius", softening_radius)

  sgamma = 3.e-22 !From Burkert 2012

  rsc = 1.2d16
  T0 = 0.4d0*obj_mu*newton*sim_ptMass*mp/kb/rsc ! Anninos 2012

  ! make sure that guardcells are up to date
  call Grid_fillGuardCells(CENTER, ALLDIR)

  ! loop over list of blocks passed in
  do thisBlock = 1, blockCount

     blockID = blockList(thisBlock)
     cooledZone = .FALSE.

     ! get dimensions/limits and coordinates
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Get a pointer to solution data 
     call Grid_getBlkPtr(blockID,solnData)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     allocate(xCoord(sizeX),stat=istat)
     sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     allocate(yCoord(sizeY),stat=istat)
     sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
     allocate(zCoord(sizeZ),stat=istat)

     if (NDIM == 3) call Grid_getCellCoords&
                         (KAXIS, blockID, CENTER, gcell, zCoord, sizeZ)
     if (NDIM >= 2) call Grid_getCellCoords&
                         (JAXIS, blockID, CENTER,gcell, yCoord, sizeY)
     call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCoord, sizeX)

     ! now guaranteed that tmp, rho, etc. exist
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        zz = zCoord(k) - sim_zCenter
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           yy = yCoord(j) - sim_yCenter
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              xx = xCoord(i) - sim_xCenter
              dist = dsqrt( xx**2 + yy**2 + zz**2 )
              cooledZone = .true.

              rho  = solnData(DENS_VAR,i,j,k)
              temp  = solnData(TEMP_VAR,i,j,k)

              xn = solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)

              call Multispecies_getSumInv(A, abar, xn)
              abar = 1.d0 / abar

              sdot = -sgamma*rho/(abar*mp)**2

              ! Simple radial conduction, saturated value from Cowie & McKee assuming an effective area ~ r.
              Tback = max (T0*(dist/rsc)**(-1.d0), sim_tAmbient)
              if (temp .lt. 0.5d0*Tback .and. dist .gt. softening_radius) then
                  sdot = sdot + sim_condCoeff*0.4d0*sqrt(2.d0*kb*Tback/PI/me)*kb*Tback/mp/dist
              endif

              ! kinetic energy
              ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  & 
                   solnData(VELY_VAR,i,j,k)**2 +  & 
                   solnData(VELZ_VAR,i,j,k)**2)

              ! internal energy, add on nuclear rate*timestep
              !ei = solnData(ENER_VAR,i,j,k) - ek
              ei = max(solnData(EINT_VAR,i,j,k) + dt*sdot, gr_smalle)
                
#ifdef EINT_VAR
              solnData(EINT_VAR,i,j,k) = ei
#endif
              solnData(ENER_VAR,i,j,k) = ei + ek
              solnData(ECOO_VAR,i,j,k) = sdot
           enddo
        enddo
     enddo

     ! we've altered the EI, let's equilabrate
     if (cooledZone) then
        call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 solnData(TEMP_VAR,i,j,k) = max(solnData(TEMP_VAR,i,j,k), sim_smallT)
              enddo
           enddo
        enddo
        call Eos_wrapped(MODE_DENS_TEMP,blkLimits,blockID)
     end if

     call Grid_releaseBlkPtr(blockID,solnData)
     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)
  end do

  return
end subroutine Cool
