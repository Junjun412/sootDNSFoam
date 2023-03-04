
!subroutine semi_driver(timeStep_,nscalar_,Y_,rho_,wmix_,nmom_,sootMom_,sootSrc_,nStats_,sootStats_)
subroutine mom_driver(timeStep_,nscalar_,Y_,rho_,wmix_,nmom_,sootMom_,sootSrc_,nStats_,sootStats_,other_,otherArr_)
  implicit none

  integer, intent(in) :: nscalar_, nmom_, nStats_
  integer, intent(in) :: other_
  double precision, intent(in) :: timeStep_
  double precision, dimension(nscalar_) :: Y_
  double precision, intent(in) :: rho_, wmix_
  double precision, dimension(nmom_), intent(in)    :: sootMom_
  double precision, dimension(nscalar_+nmom_), intent(out)   :: sootSrc_
  double precision, dimension(nStats_), intent(out) :: sootStats_
  double precision, dimension(other_), intent(out) :: otherArr_

!  write(*,*) 'Start soot_init'
  call soot_init(nscalar_,Y_,nmom_,rho_)
!  write(*,*) 'After soot_init'
  call soot_source(timeStep_,Y_,sootMom_,rho_,wmix_,sootSrc_)
!  if(sootMom_(3) > 10000.0) write(*,*) "sootMom_(3)=", sootMom_(3)
  call soot_other(other_, otherArr_)
!  write(*,*) 'After soot_source'
  call soot_poststep(sootStats_,Y_,sootMom_,rho_,wmix_)
!  write(*,*) 'After soot_poststep'
  call soot_chem_destroy
!  write(*,*) 'After soot_chem_destroy'
  call soot_destroy
!  write(*,*) 'After soot_semi_destroy'

  return
end subroutine mom_driver
