!!****f* source/physics/sourceTerms/Cool/Cool_computeDt
!!
!!
!! NAME
!!  
!!  Cool_computeDt
!!
!!
!! SYNOPSIS
!! 
!!  Cool_computeDt ( integer(IN) : blockID, 
!!                   
!!                  real,pointer :  solnData(:,:,:,:),   
!!                  real,(INOUT):   dt_cool, 
!!                  integer(INOUT): dt_minloc(:) )
!!  
!! DESCRIPTION
!!
!!  Computes the timestep limiter for heating source term solver.
!! 
!!
!!
!! ARGUMENTS
!!
!!  blockID        local block ID
!!  
!!  solnData        the physical, solution data from grid
!!  dt_check        variable to hold timestep constraint
!!  dt_minloc(5)    array to hold limiting zone info:  zone indices
!!
!!***

subroutine Cool_computeDt (blockID, &
                              blkLimits,blkLimitsGC,        &
                              solnData,   &
                              dt_cool, dt_minloc )
  use cool_data, ONLY: cool_coolDtFactor, cool_usecool, cool_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash.h"

  !! arguments
  integer, intent(IN)   :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer           :: solnData(:,:,:,:) 
  real, intent(INOUT)     :: dt_cool
  integer, intent(INOUT)  :: dt_minloc(5)

  !! local variables
  real              :: dt_temp, dt_tempInv
  integer           :: temploc(5)
  integer           :: i, j, k

  real, PARAMETER :: SMALL = TINY(1.0)
  real :: eint_zone, energyRatioInv

!!===================================================================

  ! initialize the timestep from this block to some obscenely high number

  if (.not. cool_usecool)  return

  dt_temp = HUGE(0.0)
  dt_tempInv = SMALL

  ! loop over all of the zones and compute the minimum eint/enuc
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

#ifdef EINT_VAR
           ! compute the internal energy in the zone
           eint_zone = solnData(EINT_VAR,i,j,k) 
#else
           eint_zone = solnData(ENER_VAR,i,j,k) - &
                0.5*(solnData(VELX_VAR,i,j,k)**2 + &
                solnData(VELY_VAR,i,j,k)**2 + &
                solnData(VELZ_VAR,i,j,k)**2)
#endif

           ! compute the ratio.  Note, it is the absolute value that matters.
           ! Also prevent a divide by zero by first computing and comparing
           ! the inverse of what we want, and then only (un)invert that inverse
           ! if it is a reasonable number.
           energyRatioInv = abs(solnData(ECOO_VAR,i,j,k)) / eint_zone

           if (energyRatioInv > dt_tempInv) then
              dt_tempInv = energyRatioInv
              dt_temp = 1.0 / energyRatioInv
              temploc(1) = i
              temploc(2) = j
              temploc(3) = k
              temploc(4) = blockID
              temploc(5) = cool_meshMe
           endif

        enddo
     enddo
  enddo


  ! Set the timestep from this block.
  ! A little bit of trickery to avoid multiplying HUGE by something that is > 1. - KW
  dt_temp = min( dt_temp, HUGE(0.0)/max(1.0,cool_coolDtFactor) )
  dt_temp = cool_coolDtFactor*dt_temp

  if (dt_temp < dt_cool) then
     dt_cool = dt_temp
     dt_minloc = temploc
  endif

  if(dt_cool <= 0.0) call Driver_abortFlash("[cool]: computed dt is not positive! Aborting!")

  return
end subroutine Cool_computeDt


