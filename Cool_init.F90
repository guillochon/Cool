!!****f* source/physics/sourceTerms/Cool/Cool_init
!!
!! NAME
!!
!!  Cool_init
!!
!!
!! SYNOPSIS
!!
!!  Cool_init()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***



subroutine Cool_init()
  use Cool_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype

  implicit none

#include "constants.h"
#include "Flash.h"

  call RuntimeParameters_get("useCool", cool_useCool)
  call RuntimeParameters_get("coolDtFactor", cool_coolDtFactor)

  call Driver_getMype(MESH_COMM, cool_meshMe)

  if (.not. cool_useCool) then
     write(6,*)'WARNING:  You have included the Cool unit but have set '
     write(6,*)'   the runtime parameter useCool to FALSE'
     write(6,*)'   No cooling will occur but Cool_init will continue.'
  end if
   
  return
end subroutine Cool_init
