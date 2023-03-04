!---------------------------------------------------------!
! Modified according to NGAsoot code for fdfsolver        !
! Wang Han Jul. 30 2016, advisored by Venkat Raman        !
!---------------------------------------------------------!
module soot ! harding code for Adelaide sooting flame
  use soot_mod
  implicit none

end module soot

! Part I:   Initialize the soot module - soot_init
! Part II:  Soot source Routine - soot_source
! Part III: Soot rhodot - soot_rhodot
! Part IV:  Soot poststep - soot_poststep

! ======================================================================== !
! Part I: Initialize the soot module                                       !
! ======================================================================== !
subroutine soot_init(nscalar_,Y_,nmom_,rho_)
  use soot
!  use parser
  implicit none
  integer, intent(in) :: nscalar_,nmom_,rho_
  double precision, dimension(nscalar_) :: Y_


!  call parser_read('Use mom',use_mom,2)  !0-SEMI 1-MOMIC,2-HMOM,3-DQMOM,4-CQMOM,5-ECQMOM
  use_mom = 2

!  call parser_read('Number of composition',nComposition,46) !ns+1
!  call parser_read('Number of moments',nMoments,4)
  nComposition = nscalar_
  nMoments = nmom_ !2-SEMI 6-MOMIC 7-HMOM
  nEquations = nMoments

  ! Solve lumped PAH equation
!  call parser_read('Solve PAH equation',use_pah,.false.)
  use_pah = .false.
  if (use_pah) nEquations = nEquations + 1

  ! With or without nucleation
!  call parser_read('Soot nucleation',use_nucleation,.true.)
  use_nucleation = .true.

  ! With or without coagulation
!  call parser_read('Soot coagulation',use_coagulation,.true.)
  use_coagulation = .true.

  ! With or without condensation
!  call parser_read('Soot condensation',use_condensation,.true.)
  use_condensation = .true.

  ! With or without surfacereaction
!  call parser_read('Soot surfacereaction',use_surfacereaction,.true.)
  use_surfacereaction = .true. 

  ! With or without oxidation
!  call parser_read('Soot oxidation',use_surfaceoxidation,.true.)
  use_surfaceoxidation    = .true.
  use_surfaceoxidation_OH = .true.
  use_surfaceoxidation_O  = .true.

  ! With or without fragmentation
!  call parser_read('Soot fragmentation',use_fragmentation,.true.)
  use_fragmentation = .false.

  ! With or without density coupling
!  call parser_read('Soot-density coupling',soot_rho,.true.)
  soot_rho = .false.

  ! With or without soot radiation
!  call parser_read('Soot radiation',soot_rad,.false.)
  soot_rad = .false.

  !select case(trim(chemistry))
  !case('Multistep')

  ! Initialize the specific soot model
  if (use_mom==1) then
     nSootGas = nscalar_ ! + 7 + 1 ! PAH + temperature 
     call soot_momic_init
  else if (use_mom==2) then
     nSootGas = nscalar_ !4 + 7 + 1 ! PAH + temperature 
     call soot_hmom_init(nscalar_,Y_,rho_)
  else if (use_mom==3) then
     nSootGas =  1 + 1
     call soot_dqmom_init
  else if (use_mom==4) then
!     call parser_read('Number of vol. nodes', n_nodeV, 2)
!     call parser_read('Number of area nodes', n_nodeS, 1 )
     n_nodeV = 2
     n_nodeS = 1
     nSootGas =  1 + 1
     call soot_cqmom_init
  else
     nSootGas = 7 +1   ! c2h2,h2,o2,co,OH,O,H + temperature
     call soot_semi_init(nscalar_,Y_)
  end if

  return
end subroutine soot_init



! ======================================================================== !
! Part II: Soot source Routine                                             !
! -> Compute the source terms for soot-related gas-composition and moments !
! ======================================================================== !
subroutine soot_source(timeStep,composition,moment,density,wmixture,srcSoot)
  use soot
  implicit none
  real(WP), intent(in) :: timeStep
  real(WP), intent(in), dimension(nComposition) :: composition
  real(WP), intent(in), dimension(nMoments) :: moment
  real(WP), intent(in) :: density,wmixture  ! unit ???
  real(WP), intent(out), dimension(nSootGas+nMoments) :: srcSoot

  srcSoot = 0.0_WP

  if (use_mom==1) then
     call soot_momic_source(timeStep,composition,moment,density,wmixture,srcSoot)
  else if (use_mom ==2) then
     call soot_hmom_source(timeStep,composition,moment,density,wmixture,srcSoot)
  else if (use_mom ==3) then
     call soot_dqmom_source(timeStep,composition,moment,density,wmixture,srcSoot)
  else if (use_mom ==4) then
     call soot_cqmom_source(timeStep,composition,moment,density,wmixture,srcSoot)
  else 
     call soot_semi_source(timeStep,composition,moment,density,wmixture,srcSoot)
  end if

  return
end subroutine soot_source



! ======================================================================== !
! Part I: Initialize the soot module                                       !
! ======================================================================== !
subroutine soot_rhodot(rhodot,Y_,temp_,soot_,dens_)
  use soot
  implicit none

  real(WP), intent(out) :: rhodot
  real(WP), dimension(nComposition-1), intent(in) :: Y_
  real(WP), dimension(nMoments), intent(in) :: soot_
  real(WP), intent(in) :: temp_,dens_

  rhodot = 0.0D0

  ! soot-related gas-composition
  if (use_mom==1) then
     call soot_momic_rhodot(rhodot,Y_,temp_,soot_,dens_)
  else if (use_mom ==2) then
     call soot_hmom_rhodot(rhodot,Y_,temp_,soot_,dens_)
  else if (use_mom ==3) then
     call soot_dqmom_rhodot(rhodot,Y_,temp_,soot_,dens_)
  else if (use_mom ==4) then
     call soot_cqmom_rhodot(rhodot,Y_,temp_,soot_,dens_)
  else
     call soot_semi_rhodot(rhodot,Y_,temp_,soot_,dens_)
  end if

  return
end subroutine soot_rhodot



! ======================================================================== !
! Part I: Initialize the soot module                                       !
! ======================================================================== !
subroutine soot_poststep(output,composition,soot_,dens_,wmixture)
  use soot
  implicit none

  real(WP), dimension(4+10), intent(out) :: output ! numdens,fv,partdiam,partaggr and source terms
  real(WP), dimension(nMoments), intent(in) :: soot_
  real(WP), intent(in), dimension(nComposition) :: composition
  real(WP), intent(in) :: dens_,wmixture

  output = 0.0_WP
  if (use_mom==1) then
     call soot_momic_poststep(output,composition,soot_,dens_,wmixture)
  else if (use_mom ==2) then
     call soot_hmom_poststep(output,composition,soot_,dens_,wmixture)
  else if (use_mom ==3) then
     call soot_dqmom_poststep(output,composition,soot_,dens_,wmixture)
  else if (use_mom ==4) then
     call soot_cqmom_poststep(output,composition,soot_,dens_,wmixture)
  else
     call soot_semi_poststep(output,composition,soot_,dens_,wmixture)
  end if

  return
end subroutine soot_poststep

subroutine soot_destroy
  use soot
  implicit none

  if (use_mom==1) then
     call soot_momic_destroy
  else if (use_mom ==2) then
     call soot_hmom_destroy
  else if (use_mom ==3) then
     call soot_dqmom_destroy
  else if (use_mom ==4) then
     call soot_cqmom_destroy
  else
     call soot_semi_destroy
  end if

  return

end subroutine soot_destroy

subroutine soot_other(other_,otherArr_)
  use soot
  implicit none

  integer, intent(in) :: other_
  double precision, dimension(other_), intent(out) :: otherArr_

  if (use_mom==1) then
     call soot_momic_other(other_,otherArr_)
  else if (use_mom ==2) then 
     call soot_hmom_other(other_,otherArr_)
  else

  endif

end subroutine soot_other
