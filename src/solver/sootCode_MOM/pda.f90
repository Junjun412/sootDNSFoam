!-------------------------------------------------------------
!                   PD Algorithm
! Dr. Fox's Matlab code converted into fortran subroutine
!-------------------------------------------------------------
! 	  = -1 DSTEQR failure
! istatus =  0 success
! 	  =  1 success but with one or more negative abscissas
!	  =  2 failure (q2)
!	  =  4 failure (q4)
!-------------------------------------------------------------- 
subroutine PDAlgorithm(n, mom, weights, abscissas,istatus)
!% n/2-node quadrature using PD algorithm with n/2 = 1,..,5
!% n is the number of moments [mom(1), ..., mom(n)] 
!% a are abscissas
!% w are weights
!%
  use precision
  implicit none

  integer, intent(in) :: n ! no of moments
  real(WP), dimension(n), intent(in) :: mom
  real(WP), dimension(n/2), intent(out) :: weights, abscissas

  real(WP), dimension(n/2,n/2) :: e_vector
  real(WP) :: m0, m1, m2, m3, m4, m5, m6, m7, m8, m9
  real(WP) :: q1, q2, q3, q4, q5, q6, q7, q8, q9
  real(WP), dimension(n/2) :: a  ! Diagonal of Jacobi matrix => eigenvalues => replace abscissas
  real(WP), dimension(n/2) :: b  ! Off diagonals of the Jacobi matrix
  real(WP), dimension(n/2, n/2) :: v  ! eigenvectors
  real(WP), dimension(n-2) :: tmp  ! temporary array used in dstegr routine
  integer :: info   ! Error flag for the dstegr routine
  integer :: i
  integer :: istatus

  real(WP), dimension(n/2) :: dft_ab

  istatus = 1

  ! Initializing the values
  m0=1.0_WP
  m1=0.0_WP
  m2=0.0_WP
  m3=0.0_WP
  m4=0.0_WP
  m5=0.0_WP
  m6=0.0_WP
  m7=0.0_WP
  m8=0.0_WP
  m9=0.0_WP

  ! V0+-eps in nanopart.f90 !HK! moved up
  if(n.eq.6) then
    dft_ab(1) = 1.16265211584720E-28_WP
    dft_ab(2) = 1.660931594067433E-28_WP
    dft_ab(3) = 2.15921107228766E-28_WP
  end if

  do i = 1, n/2
     weights(i) = mom(1)/real(n/2,WP)
     abscissas(i) = mom(2)/mom(1) 
  end do

  if (n/2.ge.1) then
     m0 = mom(1) 
     m1 = mom(2) 
  end if

  if (n/2.ge.2) then
     m2 = mom(3) 
     m3 = mom(4) 
  end if

  if (n/2.ge.3) then
     m4 = mom(5) 
     m5 = mom(6) 
  end if

  if (n/2.ge.4) then
     m6 = mom(7) 
     m7 = mom(8) 
  end if

  if (n/2.eq.5) then
     m8 = mom(9) 
     m9 = mom(10)
  end if


  ! Constructing the matrix
  e_vector(:,:) = 0.0_WP

  if (n/2.ge.1) then
      q1 = m1/m0 
  end if

  if (n/2.ge.2) then
      q2 = m0*m2 - (m1**2.0_WP) 
      q3 = (m0**2.0_WP)*m3 + (-2.0_WP*m0*m2 + (m1**2.0_WP))*m1 
  end if

  if (n/2.ge.3) then
      if (m2.ge.1.0E+25_WP.or.m5.ge.1.0E+50_WP) then !HK!
      !if (m2.ge.1.0E+40_WP.or.m5.ge.1.0E+150_WP) then !HK!
          istatus = 5
          !print*,'m2 m5',m2,m5
          abscissas(1:3) = dft_ab(1:3)
          return  
      end if
      q4 = -m0*(m3**2.0_WP) + (m0*m4 - (m2**2.0_WP))*m2 + (2.0_WP*m3*m2 - m4*m1)*m1    
      q7 = (m0**2.0_WP)*m3**3.0_WP+(-2.0_WP*(m0**2.0_WP)*m4*m3+(m0**2.0_WP)*m5*m2)*m2+((-2.0_WP*m0*(m3**2.0_WP)+(2.0_WP*m0*m4-(m2**2.0_WP))*m2)*m2 &
    + (2.0_WP*m3*m0*m4 + (-2.0_WP*m5*m0 + 3.0_WP*m3*m2)*m2 + (-(m3**2.0_WP) - 2.0_WP*m4*m2 + m5*m1)*m1)*m1)*m1 
  end if

  if (n/2.ge.4) then
      q5 = -m4**3.0_WP*m0+(2.0_WP*m0*m5*m4+(-m0*m6+m3**2.0_WP)*m3)*m3+(m4*m0*m6-m0*m5**2.0_WP-3.0_WP*m4*m3**2.0_WP &
      +(2.0_WP*m5*m3+m4**2.0_WP-m6*m2)*m2)*m2+((2.0_WP*m4**2.0_WP-2.0_WP*m5*m3)*m3+(-2.0_WP*m5*m4+2.0_WP*m6*m3)*m2 &
      +(m5**2.0_WP-m4*m6)*m1)*m1 
    
      q8 = (-m0**2.0_WP*m4**4.0_WP+(3.0_WP*m0**2.0_WP*m5*m4**2.0_WP+(-m0**2.0_WP*m5**2.0_WP-2.0_WP*m0**2.0_WP*m6*m4+m0**2.0_WP*m7*m3)*m3)*m3)*m3 &
      +(((-2.0_WP*m0**2.0_WP*m5**2.0_WP+2.0_WP*m0**2.0_WP*m6*m4)*m4+(2.0_WP*m0**2.0_WP*m6*m5-2.0_WP*m0**2.0_WP*m4*m7)*m3)*m3 &
      +(m0**2.0_WP*m5**3.0_WP+(-2.0_WP*m0**2.0_WP*m6*m5+m0**2.0_WP*m4*m7)*m4+(-2.0_WP*m4**3.0_WP*m0+(4.0_WP*m0*m5*m4 &
      +(-2.0_WP*m0*m6+m3**2.0_WP)*m3)*m3)*m3+(2.0_WP*m0*m5*m4**2.0_WP+(-4.0_WP*m0*m5**2.0_WP+(2.0_WP*m7*m0-4.0_WP*m4*m3)*m3)*m3 &
      +(2.0_WP*m0*m5*m6-2.0_WP*m0*m4*m7+(3.0_WP*m4**2.0_WP+3.0_WP*m5*m3)*m3+(-2.0_WP*m6*m3-2.0_WP*m5*m4+m7*m2)*m2)*m2)*m2)*m2)*m2 &
      +((2.0_WP*m4**3.0_WP*m0+(-4.0_WP*m0*m5*m4+(2.0_WP*m0*m6-m3**2.0_WP)*m3)*m3)*m3**2.0_WP+(2.0_WP*m0*m4**4.0_WP+(-8.0_WP*m0*m5*m4**2.0_WP &
      +(4.0_WP*m4*m0*m6+6.0_WP*m0*m5**2.0_WP+(-4.0_WP*m7*m0+4.0_WP*m4*m3)*m3)*m3)*m3+((2.0_WP*m0*m5**2.0_WP-2.0_WP*m4*m0*m6)*m4 &
      +(-4.0_WP*m0*m5*m6+4.0_WP*m0*m4*m7-6.0_WP*m5*m3**2.0_WP)*m3+(-2.0_WP*m4**3.0_WP+6.0_WP*m3**2.0_WP*m6+(2.0_WP*m4*m6+m5**2.0_WP-4.0_WP*m7*m3) &
      *m2)*m2)*m2)*m2+(((2.0_WP*m0*m5**2.0_WP-2.0_WP*m4*m0*m6)*m4+(2.0_WP*m0*m4*m7-2.0_WP*m0*m5*m6+(-3.0_WP*m4**2.0_WP+3.0_WP*m5*m3) &
      *m3)*m3)*m3+(-2.0_WP*m0*m5**3.0_WP+(4.0_WP*m0*m5*m6-2.0_WP*m0*m4*m7)*m4+(2.0_WP*m4**3.0_WP+(2.0_WP*m5*m4-4.0_WP*m6*m3)*m3)*m3 &
      +(m5*m4**2.0_WP+(-6.0_WP*m4*m6+m5**2.0_WP+4.0_WP*m7*m3)*m3+(2.0_WP*m7*m4-2.0_WP*m6*m5)*m2)*m2)*m2+(-m4**4.0_WP+(2.0_WP*m5*m4**2.0_WP &
      +(2.0_WP*m4*m6-3.0_WP*m5**2.0_WP)*m3)*m3+((-2.0_WP*m5**2.0_WP+2.0_WP*m4*m6)*m4+(4.0_WP*m6*m5-4.0_WP*m7*m4)*m3)*m2+(m5**3.0_WP &
      +(-2.0_WP*m6*m5+m7*m4)*m4)*m1)*m1)*m1)*m1 
  end if

  if (n/2.eq.5) then
      q6 = m5**4*m0+(-3*m6*m5**2*m0+(m6**2*m0+2*m0*m7*m5+(-m8*m0+m4**2)*m4)*m4)*m4 &
      +((2*m6**2*m0-2*m0*m7*m5)*m5+(2*m8*m0*m5-2*m6*m7*m0-4*m5*m4**2)*m4 &
      +(-m6*m8*m0+m7**2*m0+(3*m5**2+3*m4*m6)*m4+(-2*m7*m4-2*m6*m5+m8*m3)*m3)*m3)*m3 &
      +(-m6**3*m0+(2*m6*m7*m0-m8*m0*m5)*m5+(-m7**2*m0+m6*m8*m0+(3*m5**2-3*m4*m6)*m4)*m4 &
      +(-2*m5**3+(-2*m6*m5+4*m7*m4)*m4+(-3*m8*m4+2*m7*m5+m6**2)*m3)*m3 &
      +(m6*m5**2+(-4*m7*m5+2*m6**2+m8*m4)*m4+(2*m8*m5-2*m6*m7)*m3+(-m6*m8+m7**2)*m2)*m2)*m2 &
      +((-2*m5**3+(4*m6*m5-2*m7*m4)*m4)*m4+(2*m6*m5**2+(-4*m6**2+2*m8*m4)*m4 &
      +(2*m6*m7-2*m8*m5)*m3)*m3+((-2*m6**2+2*m7*m5)*m5+(2*m6*m7-2*m8*m5)*m4 &
      +(2*m6*m8-2*m7**2)*m3)*m2+(m6**3+(-2*m6*m7+m8*m5)*m5+(-m6*m8+m7**2)*m4)*m1)*m1 
    
      q9 = (m0**2*m5**5+(-4*m0**2*m6*m5**3+((3*m0**2*m6**2+3*m0**2*m7*m5)*m5+(-2*m8*m0**2*m5 &
      -2*m0**2*m6*m7+m9*m0**2*m4)*m4)*m4)*m4)*m4**2+(-m0**2*m5**6+(4*m0**2*m6*m5**4 &
      +(-6*m0**2*m7*m5**3+(-2*m0**2*m6**3+6*m8*m0**2*m5**2+(-4*m9*m0**2*m5+m0**2*m7**2+2*m8*m0**2*m6) &
      *m4)*m4)*m4)*m4+((-3*m0**2*m6**2+3*m0**2*m7*m5)*m5**3+((2*m0**2*m6**3+(2*m0**2*m6*m7-4*m8*m0**2*m5) &
      *m5)*m5+(m0**2*m6**2*m7+(-6*m8*m0**2*m6+m0**2*m7**2+4*m9*m0**2*m5)*m5+(-2*m8*m0**2*m7+2*m9*m0**2*m6) &
      *m4)*m4)*m4+(-m0**2*m6**4+(2*m0**2*m6**2*m7+(-3*m0**2*m7**2+2*m8*m0**2*m6)*m5)*m5 &
      +((-2*m0**2*m7**2+2*m8*m0**2*m6)*m6+(4*m8*m0**2*m7-4*m9*m0**2*m6-2*m0*m5**3)*m5 &
      +(6*m6*m5**2*m0+(-2*m6**2*m0-4*m0*m7*m5+(2*m8*m0-m4**2)*m4)*m4)*m4)*m4+(m0**2*m7**3 &
      +(-2*m8*m0**2*m7+m9*m0**2*m6)*m6+2*m0*m6*m5**3+((-8*m6**2*m0+2*m0*m7*m5)*m5 &
      +(-2*m8*m0*m5+8*m6*m7*m0+(-2*m9*m0+5*m5*m4)*m4)*m4)*m4+(-2*m0*m8*m5**2+2*m6**3*m0 &
      +(4*m9*m0*m5-4*m7**2*m0+(-6*m5**2-4*m4*m6)*m4)*m4+(-2*m9*m0*m6+2*m8*m0*m7+m5**3 &
      +(6*m6*m5+3*m7*m4)*m4+(-2*m8*m4-m6**2-2*m7*m5+m9*m3)*m3)*m3)*m3)*m3)*m3)*m3)*m3 &
      +(((-2*m0**2*m6**3+(4*m0**2*m6*m7-2*m8*m0**2*m5)*m5)*m5+(2*m0**2*m6**2*m7+(-4*m0**2*m7**2 &
      +2*m9*m0**2*m5)*m5+(-2*m9*m0**2*m6+2*m8*m0**2*m7)*m4)*m4)*m4**2+((2*m0**2*m6**3 &
      +(-4*m0**2*m6*m7+2*m8*m0**2*m5)*m5)*m5**2+(2*m0**2*m6**4+(-8*m0**2*m6**2*m7+(6*m0**2*m7**2 &
      +4*m8*m0**2*m6-4*m9*m0**2*m5)*m5)*m5+((2*m0**2*m7**2-2*m8*m0**2*m6)*m6+(4*m9*m0**2*m6 &
      -4*m8*m0**2*m7+4*m0*m5**3)*m5+(-12*m6*m5**2*m0+(8*m0*m7*m5+4*m6**2*m0+(-4*m8*m0+2*m4**2) &
      *m4)*m4)*m4)*m4)*m4+(((2*m0**2*m7**2-2*m8*m0**2*m6)*m6+(-2*m8*m0**2*m7+2*m9*m0**2*m6+2*m0*m5**3) &
      *m5)*m5+(-2*m0**2*m7**3+(4*m8*m0**2*m7-2*m9*m0**2*m6)*m6-12*m0*m6*m5**3+((22*m6**2*m0+2*m0*m7*m5) &
      *m5+(-20*m6*m7*m0+(6*m9*m0-10*m5*m4)*m4)*m4)*m4)*m4+((4*m6**2*m0-4*m0*m7*m5)*m5**2 &
      +(-4*m6**3*m0+(-8*m6*m7*m0+12*m8*m0*m5)*m5+(-12*m9*m0*m5+12*m7**2*m0+(8*m5**2+12*m4*m6) &
      *m4)*m4)*m4+(-2*m0*m6**2*m7+(4*m7**2*m0-2*m9*m0*m5)*m5+(8*m9*m0*m6+4*m5**3-8*m8*m0*m7 &
      +(-12*m6*m5-12*m7*m4)*m4)*m4+(-6*m6*m5**2+(4*m7*m5+2*m6**2+10*m8*m4)*m4+(2*m8*m5+2*m6*m7 &
      -6*m9*m4)*m3)*m3)*m3)*m3)*m3)*m3+((-m0**2*m6**4+(3*m0**2*m6**2*m7+(-m0**2*m7**2-2*m8*m0**2*m6 &
      +m9*m0**2*m5)*m5)*m5)*m5+(((-2*m0**2*m7**2+2*m8*m0**2*m6)*m6+(-2*m9*m0**2*m6+2*m8*m0**2*m7 &
      -4*m0*m5**3)*m5)*m5+(m0**2*m7**3+(-2*m8*m0**2*m7+m9*m0**2*m6)*m6+14*m0*m6*m5**3 &
      +((-8*m6**2*m0-10*m0*m7*m5)*m5+(4*m6*m7*m0+6*m8*m0*m5+(-2*m9*m0-m5*m4)*m4)*m4) &
      *m4)*m4)*m4+(2*m0*m6*m5**4+((-10*m6**2*m0+4*m0*m7*m5)*m5**2+(-4*m6**3*m0+(20*m6*m7*m0 &
      -10*m8*m0*m5)*m5+(4*m6*m8*m0-6*m7**2*m0+(13*m5**2-8*m4*m6)*m4)*m4)*m4)*m4 &
      +((-2*m6**3*m0+(8*m6*m7*m0-6*m8*m0*m5)*m5)*m5+(10*m0*m6**2*m7+(-8*m6*m8*m0 &
      -16*m7**2*m0+14*m9*m0*m5)*m5+(-8*m9*m0*m6-18*m5**3+8*m8*m0*m7+(-4*m6*m5+12*m7*m4)*m4) &
      *m4)*m4+((-2*m7**2*m0+2*m6*m8*m0)*m6+(4*m8*m0*m7-4*m9*m0*m6)*m5+(18*m6*m5**2 &
      +(10*m7*m5-4*m6**2-14*m8*m4)*m4)*m4+((3*m6**2-3*m7*m5)*m5+(-2*m6*m7-14*m8*m5 &
      +11*m9*m4)*m4+(4*m9*m5-2*m6*m8-m7**2)*m3)*m3)*m3)*m3)*m3+((-2*m6**2*m0+2*m0*m7*m5) &
      *m5**3+((8*m6**3*m0+(-12*m6*m7*m0+4*m8*m0*m5)*m5)*m5+(-4*m0*m6**2*m7+(10*m7**2*m0 &
      -4*m6*m8*m0-2*m9*m0*m5)*m5+(4*m9*m0*m6-4*m5**3-4*m8*m0*m7+4*m6*m5*m4)*m4)*m4)*m4 &
      +((-4*m0*m6**2*m7+(8*m6*m8*m0-4*m9*m0*m5)*m5)*m5+(8*m5**4+(-20*m7*m5+8*m6**2+4*m8*m4) &
      *m4**2)*m4+(2*m0*m7**3+(-4*m8*m0*m7+2*m9*m0*m6)*m6-6*m6*m5**3+((-8*m6**2+2*m7*m5)*m5 &
      +(-4*m6*m7+22*m8*m5-6*m9*m4)*m4)*m4+(4*m8*m5**2+(-12*m9*m5+8*m6*m8)*m4+(-2*m9*m6 &
      +2*m8*m7)*m3)*m3)*m3)*m3+(((2*m7**2*m0-2*m6*m8*m0)*m6+(-2*m8*m0*m7+2*m9*m0*m6)*m5) &
      *m5+(-2*m0*m7**3+(4*m8*m0*m7-2*m9*m0*m6)*m6-4*m6*m5**3+((-m6**2+10*m7*m5)*m5 &
      +(-4*m8*m5-2*m6*m7+m9*m4)*m4)*m4)*m4+((7*m6**2-4*m7*m5)*m5**2+(-2*m6**3+(4*m6*m7-8*m8*m5) &
      *m5+(-6*m6*m8+4*m9*m5+5*m7**2)*m4)*m4+(m6**2*m7+(-6*m6*m8+m7**2+4*m9*m5)*m5 &
      +(-6*m8*m7+6*m9*m6)*m4)*m3)*m3+((-2*m6**3+2*m6*m7*m5)*m5+(2*m6**2*m7+(4*m6*m8-6*m7**2)*m5 &
      +(-2*m9*m6+2*m8*m7)*m4)*m4+((2*m6*m8-2*m7**2)*m6+(4*m8*m7-4*m9*m6)*m5)*m3 &
      +(m7**3+(-2*m8*m7+m9*m6)*m6)*m2)*m2)*m2)*m2)*m2)*m2+((-2*m5**4*m0+(6*m6*m5**2*m0 &
      +(-2*m6**2*m0-4*m0*m7*m5+(2*m8*m0-m4**2)*m4)*m4)*m4)*m4**3+((4*m0*m6*m5**3+(-12*m6**2*m0*m5 &
      +(12*m6*m7*m0+(-4*m9*m0+6*m5*m4)*m4)*m4)*m4)*m4**2+(-2*m0*m6*m5**4+((2*m6**2*m0+4*m0*m7*m5) &
      *m5**2+(6*m6**3*m0+(-4*m6*m7*m0-8*m8*m0*m5)*m5+(-6*m7**2*m0+12*m9*m0*m5-4*m6*m8*m0 &
      +(-11*m5**2-4*m4*m6)*m4)*m4)*m4)*m4+((-4*m6*m7*m0+4*m8*m0*m5)*m5**2+(-4*m0*m6**2*m7 &
      +(4*m7**2*m0+8*m6*m8*m0-8*m9*m0*m5)*m5+(10*m5**3-4*m9*m0*m6+4*m8*m0*m7+(4*m6*m5+6*m7*m4)*m4) &
      *m4)*m4+((2*m7**2*m0-2*m6*m8*m0)*m6+(4*m9*m0*m6-4*m8*m0*m7-3*m5**3)*m5+(-6*m6*m5**2 &
      +(6*m6**2-6*m7*m5-6*m8*m4)*m4)*m4+(6*m7*m5**2+(4*m8*m5-8*m6*m7+4*m9*m4)*m4 &
      +(-4*m9*m5+2*m6*m8+m7**2)*m3)*m3)*m3)*m3)*m3)*m3+(2*m0*m5**6+(-8*m0*m6*m5**4 &
      +((8*m0*m7*m5+4*m6**2*m0)*m5**2+(4*m6**3*m0+(-8*m6*m7*m0-4*m8*m0*m5)*m5+(-4*m6*m8*m0 &
      +4*m9*m0*m5+2*m7**2*m0+(-4*m5**2+4*m4*m6)*m4)*m4)*m4)*m4)*m4+((8*m6**2*m0-8*m0*m7*m5)*m5**3 &
      +((-8*m6**3*m0+8*m0*m8*m5**2)*m5+(-4*m0*m6**2*m7+(16*m6*m8*m0-12*m9*m0*m5)*m5+(4*m5**3 &
      +(8*m6*m5-12*m7*m4)*m4)*m4)*m4)*m4+(2*m0*m6**4+(6*m7**2*m0-12*m6*m8*m0+4*m9*m0*m5)*m5**2 &
      +((4*m7**2*m0-4*m6*m8*m0)*m6+(8*m9*m0*m6-8*m8*m0*m7-2*m5**3)*m5+(-6*m6*m5**2 &
      +(-14*m6**2+8*m7*m5+14*m8*m4)*m4)*m4)*m4+(-4*m0*m7**3+(8*m8*m0*m7-4*m9*m0*m6)*m6 &
      +8*m6*m5**3+((4*m6**2-16*m7*m5)*m5+(16*m6*m7-12*m9*m4)*m4)*m4+(-2*m6**3-2*m8*m5**2 &
      +(4*m7**2-8*m6*m8+8*m9*m5)*m4+(-4*m8*m7+4*m9*m6)*m3)*m3)*m3)*m3)*m3+((-2*m6**3*m0 &
      +(4*m6*m7*m0-2*m8*m0*m5)*m5)*m5**2+(-2*m0*m6**4+(8*m0*m6**2*m7+(-6*m7**2*m0-4*m6*m8*m0 &
      +4*m9*m0*m5)*m5)*m5+((-2*m7**2*m0+2*m6*m8*m0)*m6+(4*m8*m0*m7-4*m9*m0*m6+2*m5**3)*m5 &
      +(-6*m6*m5**2+(10*m7*m5-4*m6**2-2*m8*m4)*m4)*m4)*m4)*m4+(((-4*m7**2*m0+4*m6*m8*m0)*m6 &
      +(4*m8*m0*m7-4*m9*m0*m6-4*m5**3)*m5)*m5+(4*m0*m7**3+(-8*m8*m0*m7+4*m9*m0*m6)*m6+8*m6*m5**3 &
      +((2*m6**2-2*m7*m5)*m5+(-20*m8*m5+12*m6*m7+4*m9*m4)*m4)*m4)*m4+((-14*m6**2+14*m7*m5)*m5**2 &
      +(10*m6**3+(-8*m6*m7-2*m8*m5)*m5+(16*m9*m5-16*m7**2)*m4)*m4+(-2*m6**2*m7+(12*m6*m8 &
      -2*m7**2-8*m9*m5)*m5+(-12*m9*m6+12*m8*m7)*m4)*m3)*m3)*m3+(2*m6*m5**4+((2*m6**2-8*m7*m5) &
      *m5**2+((-4*m6*m7+10*m8*m5)*m5+(-4*m9*m5+4*m6*m8-2*m7**2)*m4)*m4)*m4+((4*m6**3 &
      +(4*m8*m5-8*m6*m7)*m5)*m5+(-8*m6**2*m7+(16*m7**2-8*m9*m5)*m5)*m4+((6*m7**2-6*m6*m8)*m6 &
      +(-12*m8*m7+12*m9*m6)*m5)*m3)*m3+(m6**4+(-2*m6**2*m7+(-2*m6*m8+3*m7**2)*m5)*m5 &
      +((-2*m6*m8+2*m7**2)*m6+(-4*m8*m7+4*m9*m6)*m5)*m4+(-4*m7**3+(8*m8*m7-4*m9*m6)*m6)*m3) &
      *m2)*m2)*m2)*m2+(((2*m6**3*m0+(-4*m6*m7*m0+2*m8*m0*m5)*m5)*m5+(-2*m0*m6**2*m7 &
      +(4*m7**2*m0-2*m9*m0*m5)*m5+(2*m9*m0*m6-2*m8*m0*m7+3*m5**3+(-6*m6*m5+3*m7*m4)*m4) &
      *m4)*m4)*m4**2+((-2*m6**3*m0+(4*m6*m7*m0-2*m8*m0*m5)*m5)*m5**2+(-2*m0*m6**4+(8*m0*m6**2*m7 &
      +(-6*m7**2*m0-4*m6*m8*m0+4*m9*m0*m5)*m5)*m5+((-2*m7**2*m0+2*m6*m8*m0)*m6+(4*m8*m0*m7 &
      -4*m9*m0*m6-5*m5**3)*m5+(6*m6*m5**2+(7*m6**2-4*m7*m5-4*m8*m4)*m4)*m4)*m4)*m4+(((-2*m7**2*m0 &
      +2*m6*m8*m0)*m6+(-2*m9*m0*m6+2*m8*m0*m7+3*m5**3)*m5)*m5+(2*m0*m7**3+(-4*m8*m0*m7 &
      +2*m9*m0*m6)*m6-4*m6*m5**3+((-5*m6**2+8*m7*m5)*m5+(-10*m6*m7+4*m8*m5+4*m9*m4)*m4)*m4)*m4 &
      +((3*m6**2-6*m7*m5)*m5**2+(-4*m6**3+(12*m6*m7-2*m8*m5)*m5+(-m7**2-8*m9*m5+6*m6*m8)*m4)*m4 &
      +(3*m6**2*m7+(-6*m6*m8-3*m7**2+6*m9*m5)*m5+(-2*m9*m6+2*m8*m7)*m4)*m3)*m3)*m3)*m3 &
      +((2*m0*m6**4+(-6*m0*m6**2*m7+(2*m7**2*m0+4*m6*m8*m0-2*m9*m0*m5)*m5)*m5)*m5 &
      +(((4*m7**2*m0-4*m6*m8*m0)*m6+(4*m9*m0*m6-4*m8*m0*m7+2*m5**3)*m5)*m5+(-2*m0*m7**3 &
      +(4*m8*m0*m7-2*m9*m0*m6)*m6-6*m6*m5**3+((8*m6**2-2*m7*m5)*m5+(-6*m6*m7+4*m8*m5)*m4) &
      *m4)*m4)*m4+(-2*m6*m5**4+(6*m6**2*m5**2+(-2*m6**3+(-12*m6*m7+8*m8*m5)*m5+(14*m7**2-4*m6*m8 &
      -8*m9*m5)*m4)*m4)*m4+((2*m6*m7-2*m8*m5)*m5**2+((-4*m6*m8+2*m7**2+2*m9*m5)*m5 &
      +(-14*m8*m7+14*m9*m6)*m4)*m4+((-4*m7**2+4*m6*m8)*m6+(8*m8*m7-8*m9*m6)*m5)*m3)*m3)*m3 &
      +((-m6**2+m7*m5)*m5**3+((-6*m6**3+(14*m6*m7-8*m8*m5)*m5)*m5+(5*m6**2*m7+(-2*m6*m8-9*m7**2 &
      +6*m9*m5)*m5+(-2*m9*m6+2*m8*m7)*m4)*m4)*m4+(-3*m6**4+(10*m6**2*m7+(-2*m6*m8-9*m7**2 &
      +4*m9*m5)*m5)*m5+((-6*m7**2+6*m6*m8)*m6+(-12*m9*m6+12*m8*m7)*m5)*m4+(4*m7**3 &
      +(-8*m8*m7+4*m9*m6)*m6)*m3)*m3+(((2*m6*m8-2*m7**2)*m6+(-2*m9*m6+2*m8*m7)*m5)*m5 &
      +(2*m7**3+(-4*m8*m7+2*m9*m6)*m6)*m4)*m2)*m2)*m2+(-m5**6+(4*m6*m5**4+((-4*m6**2-2*m7*m5)*m5**2 &
      +(-2*m6**3+(8*m6*m7-2*m8*m5)*m5+(2*m6*m8-3*m7**2)*m4)*m4)*m4)*m4+((-2*m6**2+2*m7*m5)*m5**3 &
      +((4*m6**3-4*m6*m7*m5)*m5+(2*m6**2*m7+(-4*m6*m8-2*m7**2+4*m9*m5)*m5+(4*m8*m7-4*m9*m6) &
      *m4)*m4)*m4+(m6**4+(-6*m6**2*m7+(6*m6*m8+3*m7**2-4*m9*m5)*m5)*m5+((-2*m6*m8+2*m7**2)*m6 &
      +(-4*m8*m7+4*m9*m6)*m5)*m4)*m3)*m3+((2*m6**3+(-4*m6*m7+2*m8*m5)*m5)*m5**2+(2*m6**4 &
      +(-8*m6**2*m7+(6*m7**2+4*m6*m8-4*m9*m5)*m5)*m5+((-2*m6*m8+2*m7**2)*m6+(-4*m8*m7+4*m9*m6) &
      *m5)*m4)*m4+(((4*m7**2-4*m6*m8)*m6+(-4*m8*m7+4*m9*m6)*m5)*m5+(-4*m7**3+(8*m8*m7 &
      -4*m9*m6)*m6)*m4)*m3)*m2+((-m6**4+(3*m6**2*m7+(-2*m6*m8-m7**2+m9*m5)*m5)*m5)*m5 &
      +(((2*m6*m8-2*m7**2)*m6+(-2*m9*m6+2*m8*m7)*m5)*m5+(m7**3+(-2*m8*m7+m9*m6)*m6)*m4)*m4) &
      *m1)*m1)*m1)*m1 
  end if

  if (n/2.ge.1) then
      e_vector(1,1) = q1 
   end if

  if (n/2.ge.2) then
    !q2check = q2
    !if (q2.le.1.0E-10_WP) then
      if (q2.le.1.0E-30_WP) then
        !print*,'2nd-order canonical moment (q2) is invalid in PD algoritm'
        !print*, q2
          istatus = 2
          return 
       end if
       e_vector(2,2) = (q3/q2)/m0 
       e_vector(1,2) = sqrt(q2)/m0 
       e_vector(2,1) = e_vector(1,2)  
  end if

  if (n/2.ge.3) then
      !q4check = q4
      !if (q4.le.1.0E-10_WP) then
      if (q4.le.1.0E-30_WP) then
       ! print*, '4th-order canonical moment (q4) is invalid in PD algoritm'
       ! print*, q4
          istatus = 4
          return  
       end if
       e_vector(3,3) = (q7/q4)/q2  
       e_vector(2,3) = sqrt(m0*q4)/q2  
       e_vector(3,2) = e_vector(2,3)  
  end if

  if (n/2.ge.4) then
      if (q5.le.1e-30_WP) then
      !   print*, '6th-order canonical moment (q5) is invalid in PD algoritm'
      !   print*, q5
          return
      end if
      e_vector(4,4) = (q8/q5)/q4  
      e_vector(3,4) = sqrt(q2*q5)/q4  
      e_vector(4,3) = e_vector(3,4)  
  end if

  if (n/2.eq.5) then
      if (q6.le.0.0_WP) then
      !  print*, '8th-order canonical moment (q6) is invalid in PD algoritm'
      !  print*, q6
          return
      end if
      e_vector(5,5) = (q9/q6)/q5  
      e_vector(4,5) = sqrt(q4*q6)/q5  
      e_vector(5,4) = e_vector(4,5)  
  end if

  ! To Calculate the Weights & Abscissas
  ! e_vector(n/2,n/2) are the elements of the trigonal matrix
  ! a(i) are the diagonal & b(i) are the codiagonal elements
  do i = 1,n/2-1
     a(i) = e_vector(i,i) 
     b(i) = e_vector(i,i+1) 
  end do
  a(n/2) = e_vector(n/2,n/2) 
 
  v = 0.0_WP
  do i = 1,n/2
     v(i,i) = 1.0_WP
  end do
 
  ! For the LAPack tridiagonal eigenvalue solver:
  call dsteqr('I', n/2, a, b, v, n/2, tmp, info)
  !  ! For the EisPack tridiagonal eigenvalue solver:
  !  v = 0.0_WP
  !  do i = 1,n/2
  !     v(i,i) = 1.0_WP
  !  end do
  !  do i = n/2,2,-1
  !     b(i) = b(i-1)
  !  end do
  !  call imtql2(n/2, n/2, a, b, v, info)

  do i = 1,n/2
     abscissas(i) = a(i)
     weights(i) = m0*v(1,i)**2.0_WP
  end do

  istatus = 0

  ! V0+-eps in nanopart.f90
  !if(n.eq.6) then
  !  dft_ab(1) = 1.16265211584720E-28_WP
  !  dft_ab(2) = 1.660931594067433E-28_WP
  !  dft_ab(3) = 2.15921107228766E-28_WP
  !end if

  do i = 1, n/2
     if(abscissas(i).le.0.0_WP) then
        ! print *,' reached in low abscissa limit',abscissas
        if(n.eq.6) then
          abscissas(i) = dft_ab(i) ! V0 in nanopart.f90
        else
          abscissas(i) = 1E-30_WP
        end if
        istatus = 1
     end if
  end do
  
  !if(info.ne.0) then
  !  write(*,'(A,I3,6ES12.5)') ' dsteqr error:  =',info,mom
  !  istatus = -1
  !end if

  return
end subroutine PDAlgorithm
