!!****ih* source/physics/sourceTerms/Cool/CoolMain/Cool_data
!!
!! NAME
!!  Cool_data
!!
!! SYNOPSIS
!!
!!  use Cool_data
!!
!! DESCRIPTION 
!!  Cool_data is a fortran module that holds variables with
!!  the Cool Unit scope.  All variables located in Cool_data are
!!  accessible to subroutines in the Cool unit only.  (Note this 
!!  is a convention, there is nothing in the fortran language that
!!  prevents another routine from using Cool_data.  It is the FLASH3
!!  architecture that tries to enforce these rules.)
!!  
!!  All variables located in the Cool_data fortran module start with 
!! "cool_".  This is to indicate where they come from, and to make it easier
!! on the developer and user to see which variables are local variables and 
!! which belong to the Cool unit.
!!
!!***

Module Cool_data
  implicit none

  logical, save :: cool_useCool
  integer, save :: cool_meshMe
  real, save    :: cool_coolDtFactor

end Module Cool_data
