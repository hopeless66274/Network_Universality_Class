
  !! THIS IS PROJECT-DEPENDENT
  SUBROUTINE markov
    implicit none
    integer :: i
!   call Swendsen_WangV2
		call NewmanZiff
    return
  END SUBROUTINE markov
  !===================================================================
! include "Simu/NewmanZiff.f90" 
! include "Simu/NewmanZiffV3.f90" 
 include "Simu/NewmanZiffV4.f90" 

