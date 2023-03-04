!---------------------------------------------------------!
! Modified according to NGAsoot code for fdfsolver        !
! Wang Han Jul. 30 2016, advisored by Venkat Raman        !
!---------------------------------------------------------!
module soot_hmom_mod ! harding code for Adelaide sooting flame
  use soot_mod
  use soot_chem_mod
  implicit none

contains

  ! Part I: Initialize the soot module - soot_init
  ! Part II:Soot source Routine - soot_source

  ! Compute the fractional moment for cogulation, condension, surfacereaction, surfaceOxidiation, and fragmention
  ! -----------------------------
  function soot_fracmom(i,j)
    implicit none

    real(WP) :: soot_fracmom
    real(WP) :: m00, m10, m01, m20, m11, m02
    real(WP) :: f1,f2,f3,f4,f5,f6
    real(WP), intent(in) :: i,j

    m00 = mom(1) - mom(nMoments)
    m10 = mom(2) -NUCL_NBRC*mom(nMoments)
    m01 = mom(3) -NUCL_SURF*mom(nMoments)

    if (nMoments.eq.4) then
       if (m00.lt.1e-25_WP .or. m10.lt.1e-25_WP .or. m01.lt.1e-25_WP) then
          soot_fracmom = mom(1) * (mom(2)/mom(1))**i * (mom(3)/mom(1))**j
       else
          soot_fracmom = mom(nMoments)*NUCL_NBRC**i*NUCL_SURF**j + m00*(m10/m00)**i*(m01/m00)**j
       end if
    elseif (nMoments.eq.7) then
       m20 = mom(4) - (NUCL_NBRC*NUCL_NBRC)*mom(nMoments)
       m11 = mom(5) - (NUCL_NBRC*NUCL_SURF)*mom(nMoments)
       m02 = mom(6) - (NUCL_SURF*NUCL_SURF)*mom(nMoments)

       if (m00.lt.1e-25_WP .or. m10.lt.1e-25_WP .or. m01.lt.1e-25_WP .or. &
            m20.lt.1e-25_WP .or. m11.lt.1e-25_WP .or. m02.lt.1e-25) then
          f1 = mom(1)
          f2 = (mom(2)**(2.0_WP)*mom(1)**(-1.5_WP)*mom(4)**(-0.5_WP))**i
          f3 = (mom(3)**(2.0_WP)*mom(1)**(-1.5_WP)*mom(6)**(-0.5_WP))**j
          f4 = (mom(4)**(0.5_WP)*mom(1)**(0.5_WP)*mom(2)**(-1.0_WP))**(i*i)
          f5 = (mom(5)*mom(1)/mom(2)/mom(3))**(i*j)
          f6 = (mom(6)**(0.5_WP)*mom(1)**(0.5_WP)*mom(3)**(-1.0_WP))**(j*j)
          soot_fracmom = f1*f2*f3*f4*f5*f6
       else
          f1 = m00
          f2 = (m10**(2.0_WP)*m00**(-1.5_WP)*m20**(-0.5_WP))**i
          f3 = (m01**(2.0_WP)*m00**(-1.5_WP)*m02**(-0.5_WP))**j
          f4 = (m20**(0.5_WP)*m00**(0.5_WP)*m10**(-1.0_WP))**(i*i)
          f5 = (m11*m00/m10/m01)**(i*j)
          f6 = (m02**(0.5_WP)*m00**(0.5_WP)*m01**(-1.0_WP))**(j*j)
          soot_fracmom = mom(nMoments)*NUCL_NBRC**i*NUCL_SURF**j + f1*f2*f3*f4*f5*f6
       end if
    end if

    return
  end function soot_fracmom




  ! Compute the fractional moment of the second mode
  ! ------------------------------------------------
  function soot_fracmomlarge(i,j)
    implicit none
    real(WP) :: soot_fracmomlarge
    real(WP) :: m00, m10, m01, m20, m11, m02
    real(WP), intent(in) :: i,j
    integer :: n

    m00 = mom(1) - mom(nMoments)
    m10 = mom(2) -(NUCL_NBRC)*mom(nMoments)
    m01 = mom(3) -(NUCL_SURF)*mom(nMoments)

    if (nMoments.eq.4) then
       if (m00.lt.1e-25_WP .or. m10.lt.1e-25_WP .or. m01.lt.1e-25_WP) then
          soot_fracmomlarge = 1.0e-60_WP * NUCL_NBRC**i * NUCL_SURF**j
       else
          soot_fracmomlarge = soot_fracmom(i,j) - mom(nMoments)*NUCL_NBRC**i*NUCL_SURF**j
       end if
    elseif (nMoments.eq.7) then
       m20 = mom(4) - (NUCL_NBRC*NUCL_NBRC)*mom(nMoments)
       m11 = mom(5) - (NUCL_NBRC*NUCL_SURF)*mom(nMoments)
       m02 = mom(6) - (NUCL_SURF*NUCL_SURF)*mom(nMoments)

       if (m00.lt.1e-25_WP .or. m10.lt.1e-25_WP .or. m01.lt.1e-25_WP .or. &
            m20.lt.1e-25_WP .or. m11.lt.1e-25_WP .or. m02.lt.1e-25) then
          soot_fracmomlarge = 1.0e-60_WP * NUCL_NBRC**i * NUCL_SURF**j
       else
          soot_fracmomlarge = soot_fracmom(i,j) - mom(nMoments)*NUCL_NBRC**i*NUCL_SURF**j
       end if
    end if

    return
  end function soot_fracmomlarge



  ! Compute the collision coefficients
  ! ----------------------------------
  function soot_beta_nucl()
    implicit none
    real(WP) :: soot_beta_nucl

    ! 2.2_WP * 4.0_WP * sqrt(2.0_WP) = 12.45
    soot_beta_nucl = 12.45_WP * Cfm * sqrtT * DIMER_NBRC**(1.0_WP/6.0_WP)

    return
  end function soot_beta_nucl



  function soot_beta_cond()
    implicit none
    real(WP) :: soot_beta_cond
    real(WP) :: dV,S0v

    dV = DIMER_NBRC

    S0v =          soot_fracmom(2.0_WP*av       , 2.0_WP*as       ) * dV**(-3.0_WP/6.0_WP) + &
          2.0_WP * soot_fracmom(       av       ,        as       ) * dV**(-1.0_WP/6.0_WP) + &
                   soot_fracmom(          0.0_WP,           0.0_WP) * dV**(1.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(2.0_WP*av-1.0_WP, 2.0_WP*as       ) * dV**(3.0_WP/6.0_WP) + &
                   soot_fracmom(       av-1.0_WP,        as       ) * dV**(5.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(         -1.0_WP,           0.0_WP) * dV**(7.0_WP/6.0_WP)

    soot_beta_cond = Cfm * sqrtT * S0v

    return
  end function soot_beta_cond



  ! Create the moments
  ! ------------------
  subroutine soot_define_moments
    implicit none

    allocate(moments(nMoments,dim))

    if (nMoments.ge.4) then
       ! Total number
       moments(1,1) = 0.0_WP
       moments(1,2) = 0.0_WP
       ! Total volume
       moments(2,1) = 1.0_WP
       moments(2,2) = 0.0_WP
       ! Total surface
       moments(3,1) = 0.0_WP
       moments(3,2) = 1.0_WP
    end if
    if (nMoments.ge.7) then
       ! Total volume*volume
       moments(4,1) = 2.0_WP
       moments(4,2) = 0.0_WP
       ! Total volume*surface
       moments(5,1) = 1.0_WP
       moments(5,2) = 1.0_WP
       ! Total surface*surface
       moments(6,1) = 0.0_WP
       moments(6,2) = 2.0_WP
    end if

    return
  end subroutine soot_define_moments



  ! 2.2) =============== soot_hmom_compute_src ====================!
  ! Compute the source terms
  ! ------------------------
  subroutine soot_hmom_compute_src
    implicit none

    ! Compute PAH molecules
    call soot_pah_molecules
!    if (betaC_out .gt. 2.5e3) then
!            write(*,*) betaC_out
!            write(*,*) mom(1)
!            write(*,*) mom(2)
!            write(*,*) mom(3)
!            write(*,*) mom(4)
!    end if

    ! Compute source terms
    call soot_hmom_src_mom

    return
  end subroutine soot_hmom_compute_src



  ! 2.2.1) Compute PAH concentration to get DIMER_conc for nucleation and condension
  ! -------------------------
  subroutine soot_pah_molecules
    implicit none
    real(WP) :: prodRate_local
    real(WP) :: betaN,betaC,delta
    integer :: n

    prodRate_local = prodRate

    ! Nucleation collision coeff
    betaN = soot_beta_nucl()

    ! Condensation collision coeff
    betaC = soot_beta_cond()

    ! Solve quadratic equation
    delta = betaC**2 + 4.0_WP*betaN*prodRate_local
    DIMER_conc = (sqrt(delta)-betaC)/(2.0_WP*betaN)

    betaN_out = betaN
    betaC_out = betaC
    prodRate_out = prodRate_local
    return
  end subroutine soot_pah_molecules




  ! 2.2.2) Compute soot formation and growth process for source terms of moments
  ! ----------------------------------------------------------------------------
  subroutine soot_hmom_src_mom
    implicit none
    integer :: i
    real(WP) :: src_tmp

    ! Source terms for moments -- M00 M10 M01
    do i=1,nMoments-1
       kmom = moments(i,:)
       src_mom(i) = 0.0_WP
       ! Nucleation
       if (use_nucleation) then
           call soot_nucleation(kmom,src_tmp)
           src_mom(i) = src_mom(i) + src_tmp
       end if
       ! Coagulation
       if (use_coagulation) then
           call soot_coagulation(kmom,src_tmp)
           src_mom(i) = src_mom(i) + src_tmp
       end if
       ! Condensation
       if (use_condensation) then
           call soot_condensation(kmom,src_tmp)
           src_mom(i) = src_mom(i) + src_tmp
       end if
       ! Surface Reaction
       if (use_surfacereaction) then
           call soot_surfacereaction(kmom,src_tmp)
           src_mom(i) = src_mom(i) + src_tmp
       end if
       ! Surface Oxidation
       if (use_surfaceoxidation) then
           call soot_surfaceoxidation(kmom,src_tmp)
           src_mom(i) = src_mom(i) + src_tmp
       end if
       ! Fragmentation
       if (use_fragmentation) then
           call soot_fragmentation(kmom,src_tmp)
           src_mom(i) = src_mom(i) + src_tmp
       end if
    end do

    ! Source term for weight of first peak -- N00
    src_mom(nMoments) = 0.0_WP
    ! Nucleation
    if (use_nucleation) then
        call soot_nucleationsmall(src_tmp)
        src_mom(nMoments) = src_mom(nMoments) + src_tmp
    end if
    ! Coagulation
    if (use_coagulation) then
        call soot_coagulationsmall(src_tmp)
        src_mom(nMoments) = src_mom(nMoments) + src_tmp
    end if
    ! Condensation
    if (use_condensation) then
        call soot_condensationsmall(src_tmp)
        src_mom(nMoments) = src_mom(nMoments) + src_tmp
    end if
    ! Surface Reaction
    if (use_surfacereaction) then
        call soot_surfacereactionsmall(src_tmp)
        src_mom(nMoments) = src_mom(nMoments) + src_tmp
    end if
    ! Surface Oxidation
    if (use_surfaceoxidation) then
        call soot_surfaceoxidationsmall(src_tmp)
        src_mom(nMoments) = src_mom(nMoments) + src_tmp
    end if
    ! Fragmentation
    if (use_fragmentation) then
        call soot_fragmentationsmall(src_tmp)
        src_mom(nMoments) = src_mom(nMoments) + src_tmp
    end if

    return
  end subroutine soot_hmom_src_mom



  !=== Begine HMOM Model=====================================================================================!
  !=================!
  ! Nucleation      !
  !=================!
  ! Compute Nucleation source term for moment k
  ! -------------------------------------------
  subroutine soot_nucleation(k,src)
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src

    src = 0.5_WP * soot_beta_nucl() * DIMER_conc**2 * &
          NUCL_NBRC**k(1) * NUCL_SURF**k(2)

    return
  end subroutine soot_nucleation




  ! Compute Nucleation source term for weight of first mode
  ! -------------------------------------------------------
  subroutine soot_nucleationsmall(src)
    implicit none

    real(WP), intent(out) :: src

    src = 0.5_WP * soot_beta_nucl() * DIMER_conc**2

    return
  end subroutine soot_nucleationsmall




  !=================!
  ! Coagulation     !
  !=================!

  ! Compute Coagulation source term for moment k
  ! --------------------------------------------
  subroutine soot_coagulation(k,src)
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: V0,S0,Vi,Si
    real(WP) :: ss_fm,ss_cn,ss
    real(WP) :: sl_fm,sl_cn,sl
    real(WP) :: ll_fm,ll_cn,ll
    real(WP) :: slip

    V0 = NUCL_NBRC
    S0 = NUCL_SURF
    Vi = soot_fracmomlarge(1.0_WP,0.0_WP)/soot_fracmomlarge(0.0_WP,0.0_WP)
    Si = soot_fracmomlarge(0.0_WP,1.0_WP)/soot_fracmomlarge(0.0_WP,0.0_WP)

    slip = lambda * MUsqrtW_RHOsqrtT

    if ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
       ss = 0.0_WP
    else
       ss_fm = 2.2_WP * Cfm * sqrtT * 2.0_WP**2.5_WP * (2.0_WP**(k(1)+2.0_WP/3.0_WP*k(2)-1.0_WP) - 1.0_WP) * &
               V0**(k(1)+2.0_WP/3.0_WP*k(2)+1.0_WP/6.0_WP) * mom(nMoments) * mom(nMoments)
       ss_cn = 4.0_WP * Ccn * T_MU * (1.0_WP + 1.257_WP*slip*V0**(-1.0_WP/3.0_WP)) * (2.0_WP**(k(1)+2.0_WP/3.0_WP*k(2)-1.0_WP) - 1.0_WP) * &
               V0**(k(1)+2.0_WP/3.0_WP*k(2)) * mom(nMoments) * mom(nMoments)
       ss = (ss_fm * ss_cn) / (ss_fm + ss_cn)
    end if

    if ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
       sl = 0.0_WP
    else
       if ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
          sl_fm = -soot_psi_sl(0.0_WP,0.0_WP,0.0_WP,0.0_WP)
          sl_cn = -Ccn*T_MU * (2.000_WP                                 *mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                             NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                             NUCL_NBRC**(-1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                               1.257_WP*slip*NUCL_NBRC**(-1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                               1.257_WP*slip                            *mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                               1.257_WP*slip*NUCL_NBRC**(-2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                               1.257_WP*slip*NUCL_NBRC**(1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as))

       elseif ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.1.01_WP .and. k(2).gt.0.99_WP)) then
          sl_fm = FitC*soot_psi_sl(-2.0_WP*FitExp-1.0_WP,3.0_WP*FitExp+1.0_WP,1.0_WP,0.0_WP) - &
                       soot_psi_sl(               0.0_WP, 0.0_WP,0.0_WP,1.0_WP)
          sl_cn = Ccn*T_MU * ( FitC * ( &
                                       2.000_WP*     NUCL_NBRC*mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP, 3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**(4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-as+3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**(2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**(2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC*mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**(1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**(4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-2.0_WP*FitExp-1.0_WP,-2.0_WP*as+3.0_WP*FitExp+1.0_WP)&
                             ) - ( &
                                       2.000_WP*     NUCL_NBRC**(2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                                     NUCL_NBRC*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                                     NUCL_NBRC**(1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                       1.257_WP*slip*NUCL_NBRC**(1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                       1.257_WP*slip*NUCL_NBRC**(2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                       1.257_WP*slip*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                       1.257_WP*slip*NUCL_NBRC*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as)))
       elseif ((k(1).lt.2.01_WP .and. k(1).gt.1.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
          sl_fm = 2.0_WP*soot_psi_sl(1.0_WP,0.0_WP,1.0_WP,0.0_WP)
          sl_cn = 2.0_WP*Ccn*T_MU * (2.000_WP*     NUCL_NBRC                 *mom(nMoments)*soot_fracmomlarge(           1.0_WP, 0.0_WP   ) + &
                                                   NUCL_NBRC**(4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av+1.0_WP,-       as) + &
                                                   NUCL_NBRC**(2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av+1.0_WP,        as) + &
                                     1.257_WP*slip*NUCL_NBRC**(2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          +1.0_WP, 0.0_WP   ) + &
                                     1.257_WP*slip*NUCL_NBRC*mom(nMoments)*soot_fracmomlarge(-       av+1.0_WP,-       as) + &
                                     1.257_WP*slip*NUCL_NBRC**(1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av+1.0_WP,        as) + &
                                     1.257_WP*slip*NUCL_NBRC**(4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av+1.0_WP,-2.0_WP*as))
       elseif ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.1.01_WP .and. k(2).gt.0.99_WP)) then
          sl_fm = FitC*soot_psi_sl(-2.0_WP*FitExp       ,3.0_WP*FitExp+1.0_WP,1.0_WP,0.0_WP) + &
                       soot_psi_sl(               0.0_WP,              1.0_WP,1.0_WP,0.0_WP) + &
                  FitC*soot_psi_sl(-2.0_WP*FitExp-1.0_WP,3.0_WP*FitExp+1.0_WP,2.0_WP,0.0_WP) - &
                       soot_psi_sl(               0.0_WP,              0.0_WP,1.0_WP,1.0_WP)
          sl_cn = Ccn*T_MU * ( FitC * ( &
                                       2.000_WP*     NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp,           3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp,-       as+3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp,           3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp,-       as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-2.0_WP*FitExp,-2.0_WP*as+3.0_WP*FitExp+1.0_WP) &
                             ) + ( &
                                       2.000_WP*     NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                                                     NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                                                     NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as+1.0_WP) &
                             ) + FitC * ( &
                                       2.000_WP*     NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 7.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 7.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-2.0_WP*FitExp-1.0_WP,-2.0_WP*as+3.0_WP*FitExp+1.0_WP) &
                             ) - ( &
                                       2.000_WP*     NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                                     NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                                     NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                       1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                       1.257_WP*slip*NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                       1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                       1.257_WP*slip*NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as)))
       elseif ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.2.01_WP .and. k(2).gt.1.99_WP)) then
          sl_fm = 2.0_WP*FitC*     soot_psi_sl(-2.0_WP*FitExp-1.0_WP,3.0_WP*FitExp+2.0_WP,1.0_WP,0.0_WP) + &
                         FitC*FitC*soot_psi_sl(-4.0_WP*FitExp-2.0_WP,6.0_WP*FitExp+2.0_WP,2.0_WP,0.0_WP) - &
                                   soot_psi_sl(               0.0_WP,              0.0_WP,0.0_WP,2.0_WP)
          sl_cn = 2.0_WP*Ccn*T_MU * ( 2.0_WP*FitC * ( &
                                              2.000_WP*     NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+2.0_WP) + &
                                                            NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+2.0_WP) + &
                                                            NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-2.0_WP*FitExp-1.0_WP,-2.0_WP*as+3.0_WP*FitExp+2.0_WP) &
                                    ) + FitC*FitC * ( &
                                              2.000_WP*     NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(          -4.0_WP*FitExp-2.0_WP,           6.0_WP*FitExp+2.0_WP) + &
                                                            NUCL_NBRC**( 7.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-4.0_WP*FitExp-2.0_WP,-       as+6.0_WP*FitExp+2.0_WP) + &
                                                            NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-4.0_WP*FitExp-2.0_WP,        as+6.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -4.0_WP*FitExp-2.0_WP,           6.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(-       av-4.0_WP*FitExp-2.0_WP,-       as+6.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-4.0_WP*FitExp-2.0_WP,        as+6.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 7.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-4.0_WP*FitExp-2.0_WP,-2.0_WP*as+6.0_WP*FitExp+2.0_WP) &
                                    ) - ( &
                                              2.000_WP*     NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                                            NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                                            NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                              1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                              1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                              1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                              1.257_WP*slip*NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as)))
       endif
       sl = (sl_fm * sl_cn) / (sl_fm + sl_cn)
    end if

    if ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
       ll = 0.0_WP
    elseif ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.1.01_WP .and. k(2).gt.0.99_WP)) then
       ll = 0.0_WP
    else
       if ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
          ll_fm = -0.5_WP*soot_psi_ll(0.0_WP,0.0_WP,0.0_WP,0.0_WP)
          ll_cn = -0.5*Ccn*T_MU * (2.000_WP*     soot_fracmomlarge( 0.0_WP   , 0.0_WP   )*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                                 soot_fracmomlarge(        av,        as)*soot_fracmomlarge(-       av,-       as) + &
                                                 soot_fracmomlarge(-       av,-       as)*soot_fracmomlarge(        av,        as) + &
                                   1.257_WP*slip*soot_fracmomlarge(-       av,-       as)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                   1.257_WP*slip*soot_fracmomlarge( 0.0_WP   , 0.0_WP   )*soot_fracmomlarge(-       av,-       as) + &
                                   1.257_WP*slip*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as)*soot_fracmomlarge(        av,        as) + &
                                   1.257_WP*slip*soot_fracmomlarge(        av,        as)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as))
       elseif ((k(1).lt.2.01_WP .and. k(1).gt.1.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
          ll_fm = soot_psi_ll(1.0_WP,0.0_WP,1.0_WP,0.0_WP)
          ll_cn = Ccn*T_MU * (2.000_WP*     soot_fracmomlarge(           1.0_WP, 0.0_WP   )*soot_fracmomlarge(           1.0_WP, 0.0_WP   ) + &
                                            soot_fracmomlarge(        av+1.0_WP,        as)*soot_fracmomlarge(-       av+1.0_WP,-       as) + &
                                            soot_fracmomlarge(-       av+1.0_WP,-       as)*soot_fracmomlarge(        av+1.0_WP,        as) + &
                              1.257_WP*slip*soot_fracmomlarge(-       av+1.0_WP,-       as)*soot_fracmomlarge(           1.0_WP, 0.0_WP   ) + &
                              1.257_WP*slip*soot_fracmomlarge(           1.0_WP, 0.0_WP   )*soot_fracmomlarge(-       av+1.0_WP,-       as) + &
                              1.257_WP*slip*soot_fracmomlarge(-2.0_WP*av+1.0_WP,-2.0_WP*as)*soot_fracmomlarge(        av+1.0_WP,        as) + &
                              1.257_WP*slip*soot_fracmomlarge(        av+1.0_WP,        as)*soot_fracmomlarge(-2.0_WP*av+1.0_WP,-2.0_WP*as))
       elseif ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.1.01_WP .and. k(2).gt.0.99_WP)) then
          ll_fm = soot_psi_ll(1.0_WP,0.0_WP,0.0_WP,1.0_WP)
          ll_cn = Ccn*T_MU * (2.000_WP*     soot_fracmomlarge(           1.0_WP, 0.0_WP   )*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                                            soot_fracmomlarge(        av+1.0_WP,        as)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                                            soot_fracmomlarge(-       av+1.0_WP,-       as)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(-       av+1.0_WP,-       as)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(           1.0_WP, 0.0_WP   )*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(-2.0_WP*av+1.0_WP,-2.0_WP*as)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(        av+1.0_WP,        as)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as+1.0_WP))
       elseif ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.2.01_WP .and. k(2).gt.1.99_WP)) then
          ll_fm = soot_psi_ll(0.0_WP,1.0_WP,0.0_WP,1.0_WP)
          ll_cn = Ccn*T_MU * (2.000_WP*     soot_fracmomlarge( 0.0_WP   ,           1.0_WP)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                                            soot_fracmomlarge(        av,        as+1.0_WP)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                                            soot_fracmomlarge(-       av,-       as+1.0_WP)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(-       av,-       as+1.0_WP)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge( 0.0_WP   ,           1.0_WP)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as+1.0_WP)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(        av,        as+1.0_WP)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as+1.0_WP))
       endif
       ll = (ll_fm * ll_cn) / (ll_fm + ll_cn)
    end if

    src = ss + sl + ll

    return
  end subroutine soot_coagulation




  ! Compute Coagulation source term for weight of first mode
  ! --------------------------------------------------------
  subroutine soot_coagulationsmall(src)
    !use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: V0,Vi,Si
    real(WP) :: ss_fm,ss_cn,ss
    real(WP) :: sl_fm,sl_cn,sl
    real(WP) :: slip

    V0 = NUCL_NBRC
    Vi = soot_fracmomlarge(1.0_WP,0.0_WP)/soot_fracmomlarge(0.0_WP,0.0_WP)
    Si = soot_fracmomlarge(0.0_WP,1.0_WP)/soot_fracmomlarge(0.0_WP,0.0_WP)

    slip = lambda * MUsqrtW_RHOsqrtT

    ss_fm = -2.2_WP * Cfm * sqrtT * 2.0_WP**2.5_WP * V0**(1.0_WP/6.0_WP) * mom(nMoments) * mom(nMoments)
    ss_cn = -4.0_WP * Ccn * T_MU * (1.0_WP + 1.257_WP*slip*V0**(-1.0_WP/3.0_WP)) * mom(nMoments) * mom(nMoments)
    ss = (ss_fm * ss_cn) / (ss_fm + ss_cn)

    sl_fm = -soot_psi_sl(0.0_WP,0.0_WP,0.0_WP,0.0_WP)
    sl_cn = -Ccn * T_MU * (2.000_WP*                                 mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                         NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                         NUCL_NBRC**(-1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                           1.257_WP*slip*NUCL_NBRC**(-1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                           1.257_WP*slip*                            mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                           1.257_WP*slip*NUCL_NBRC**(-2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                           1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as))
    sl = (sl_fm * sl_cn) / (sl_fm + sl_cn)

    src = ss + sl

    return
  end subroutine soot_coagulationsmall




  ! Compute PsiSL for free molecular coagulation source term
  ! --------------------------------------------------------
  function soot_psi_sl(x,y,a,b)
    !use soot
    implicit none

    real(WP), intent(in) :: x,y,a,b
    real(WP) :: soot_psi_sl
    real(WP) :: psi1,psi2,psi3,psi4

    psi1 =          NUCL_NBRC**( 1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(         -0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**(-1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av-0.5_WP+x,       as+y) + &
                    NUCL_NBRC**(-3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y)

    psi2 =          NUCL_NBRC**( 7.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(         -0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**( 5.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av-0.5_WP+x,       as+y) + &
                    NUCL_NBRC**( 3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) + &
                    NUCL_NBRC**( 1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(          0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**(-1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av+0.5_WP+x,       as+y) + &
                    NUCL_NBRC**(-3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av+0.5_WP+x,2.0_WP*as+y)

    psi3 =          NUCL_NBRC**( 13.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(         -0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**( 11.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av-0.5_WP+x,       as+y) + &
                    NUCL_NBRC**(  9.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) + &
                    NUCL_NBRC**(  7.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(          0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**(  5.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av+0.5_WP+x,       as+y) + &
                    NUCL_NBRC**(  3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av+0.5_WP+x,2.0_WP*as+y) + &
                    NUCL_NBRC**(  1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(          1.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**( -1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av+1.5_WP+x,       as+y) + &
                    NUCL_NBRC**( -3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av+1.5_WP+x,2.0_WP*as+y)
    if (nMoments.gt.4) then
    psi4 =          NUCL_NBRC**( 19.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(         -0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**( 17.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av-0.5_WP+x,       as+y) + &
                    NUCL_NBRC**( 15.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) + &
                    NUCL_NBRC**( 13.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(          0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**( 11.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av+0.5_WP+x,       as+y) + &
                    NUCL_NBRC**(  9.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av+0.5_WP+x,2.0_WP*as+y) + &
                    NUCL_NBRC**(  7.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(          1.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**(  5.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av+1.5_WP+x,       as+y) + &
                    NUCL_NBRC**(  3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av+1.5_WP+x,2.0_WP*as+y) + &
                    NUCL_NBRC**(  1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(          2.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**( -1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av+2.5_WP+x,       as+y) + &
                    NUCL_NBRC**( -3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av+2.5_WP+x,2.0_WP*as+y)
    end if

    !if (nMoments==4) then
    !   soot_psi_sl = 2.2_WP * Cfm * sqrtT * psi1**(3.0_WP/8.0_WP)*psi2**(3.0_WP/4.0_WP) * psi3**(-1.0_WP/8.0_WP)
    !else if (nMoments==7) then
    !   soot_psi_sl = 2.2_WP * Cfm * sqrtT * psi1**(5.0_WP/16.0_WP)*psi2**(15_WP/16_WP) * psi3**(-5.0_WP/16.0_WP) * psi4**(1.0_WP/16.0_WP)
    !end if
    soot_psi_sl = 2.2_WP * Cfm * sqrtT * sqrt(psi1*psi2)

    return
  end function soot_psi_sl




  ! Compute PsiLL for free molecular coagulation source term
  ! --------------------------------------------------------
  function soot_psi_ll(x,y,a,b)
    !use soot
    implicit none

    real(WP), intent(in) :: x,y,a,b
    real(WP) :: soot_psi_ll
    real(WP) :: psi1,psi2,psi3,psi4

    psi1 =          soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(         -0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av-0.5_WP+x,       as+y) * soot_fracmomlarge(       av-0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         -0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av-0.5_WP+a,2.0_WP*as+b)

    psi2 =          soot_fracmomlarge(2.0_WP*av+0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(         -0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av+0.5_WP+x,       as+y) * soot_fracmomlarge(       av-0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         +0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av-0.5_WP+a,2.0_WP*as+b) + &
                    soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(          0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av-0.5_WP+x,       as+y) * soot_fracmomlarge(       av+0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         -0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av+0.5_WP+a,2.0_WP*as+b)

    psi3 =          soot_fracmomlarge(2.0_WP*av+1.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(         -0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av+1.5_WP+x,       as+y) * soot_fracmomlarge(       av-0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         +1.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av-0.5_WP+a,2.0_WP*as+b) + &
                    soot_fracmomlarge(2.0_WP*av+0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(          0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av+0.5_WP+x,       as+y) * soot_fracmomlarge(       av+0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(          0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av+0.5_WP+a,2.0_WP*as+b) + &
                    soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(          1.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av-0.5_WP+x,       as+y) * soot_fracmomlarge(       av+1.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         -0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av+1.5_WP+a,2.0_WP*as+b)
    if (nMoments.gt.4) then
    psi4 =          soot_fracmomlarge(2.0_WP*av+2.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(         -0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av+2.5_WP+x,       as+y) * soot_fracmomlarge(       av-0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         +2.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av-0.5_WP+a,2.0_WP*as+b) + &
                    soot_fracmomlarge(2.0_WP*av+1.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(          0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av+1.5_WP+x,       as+y) * soot_fracmomlarge(       av+0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         +1.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av+0.5_WP+a,2.0_WP*as+b) + &
                    soot_fracmomlarge(2.0_WP*av+0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(          1.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av+0.5_WP+x,       as+y) * soot_fracmomlarge(       av+1.5_WP+a,       as+b) + &
                    soot_fracmomlarge(          0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av+1.5_WP+a,2.0_WP*as+b) + &
                    soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(          2.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av-0.5_WP+x,       as+y) * soot_fracmomlarge(       av+2.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         -0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av+2.5_WP+a,2.0_WP*as+b)
    end if

    !if (nMoments==4) then
    !   soot_psi_ll = 2.2_WP * Cfm * sqrtT * psi1**(3.0_WP/8.0_WP)*psi2**(3.0_WP/4.0_WP) * psi3**(-1.0_WP/8.0_WP)
    !else if (nMoments==7) then
    !   soot_psi_ll = 2.2_WP * Cfm * sqrtT * psi1**(5.0_WP/16.0_WP)*psi2**(15_WP/16_WP) * psi3**(-5.0_WP/16.0_WP) * psi4**(1.0_WP/16.0_WP)
    !end if

    soot_psi_ll = 2.2_WP * Cfm * sqrtT * sqrt(psi1*psi2)

    return
  end function soot_psi_ll




  !=================!
  ! Condensation    !
  !=================!
  ! Compute Condensation source term for moment k
  ! ---------------------------------------------
  subroutine soot_condensation(k,src)
    !use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: dV
    real(WP) :: S0v,S0s

    dV = DIMER_NBRC

    S0v =          soot_fracmom(k(1)+2.0_WP*av-1.0_WP,k(2)+2.0_WP*as) * dV**( 3.0_WP/6.0_WP) + &
          2.0_WP * soot_fracmom(k(1)+       av-1.0_WP,k(2)+       as) * dV**( 5.0_WP/6.0_WP) + &
                   soot_fracmom(k(1)          -1.0_WP,k(2)          ) * dV**( 7.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(k(1)+2.0_WP*av-2.0_WP,k(2)+2.0_WP*as) * dV**( 9.0_WP/6.0_WP) + &
                   soot_fracmom(k(1)+       av-2.0_WP,k(2)+       as) * dV**(11.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(k(1)          -2.0_WP,k(2)          ) * dV**(13.0_WP/6.0_WP)

    S0s =          soot_fracmom(k(1)-2.0_WP*FitExp+2.0_WP*av-1.0_WP, k(2)+3.0_WP*FitExp+2.0_WP*as) * dV**( 3.0_WP/6.0_WP) + &
          2.0_WP * soot_fracmom(k(1)-2.0_WP*FitExp+       av-1.0_WP, k(2)+3.0_WP*FitExp+       as) * dV**( 5.0_WP/6.0_WP) + &
                   soot_fracmom(k(1)-2.0_WP*FitExp          -1.0_WP, k(2)+3.0_WP*FitExp          ) * dV**( 7.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(k(1)-2.0_WP*FitExp+2.0_WP*av-2.0_WP, k(2)+3.0_WP*FitExp+2.0_WP*as) * dV**( 9.0_WP/6.0_WP) + &
                   soot_fracmom(k(1)-2.0_WP*FitExp+       av-2.0_WP, k(2)+3.0_WP*FitExp+       as) * dV**(11.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(k(1)-2.0_WP*FitExp          -2.0_WP, k(2)+3.0_WP*FitExp          ) * dV**(13.0_WP/6.0_WP)

    src =  Cfm * sqrtT * (k(1) * S0v + FitC * k(2) * S0s) * DIMER_conc

    return
  end subroutine soot_condensation




  ! Compute Condensation source term for weight of first mode
  ! ---------------------------------------------------------
  subroutine soot_condensationsmall(src)
    !use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: dV,V0

    dV = DIMER_NBRC
    V0 = NUCL_NBRC

    src = -Cfm * sqrtT * sqrt(1.0_WP/dV + 1.0_WP/V0) * (dV**(1.0_WP/3.0_WP) + V0**(1.0_WP/3.0_WP))**2.0_WP * DIMER_conc * mom(nMoments)

    return
  end subroutine soot_condensationsmall




  !=================!
  ! SurfaceReaction !
  !=================!
  ! Compute Surface Reaction source term for moment k
  ! -------------------------------------------------
  subroutine soot_surfacereaction(k,src)
    !use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: dV,chicarb,V0
    
!    dV = 1.0_WP
    dV =  NUCL_NBRC**(2.0_WP/3.0_WP) / min_nbrC
    chicarb = chisoot * (36.0_WP*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    src = wCoeff * chicarb * dV * (k(1) * soot_fracmom(k(1)-1.0_WP, k(2)+1.0_WP) + &
          k(2) * FitC * soot_fracmom(k(1)-1.0_WP-2.0_WP*FitExp,k(2)+1.0_WP+3.0_WP*FitExp))
   
    return
  end subroutine soot_surfacereaction




  ! Compute Surface Reaction source term for moment k
  ! -------------------------------------------------
  subroutine soot_surfacereactionsmall(src)
   ! use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: dV,V0,chicarb

!    dV = 1.0_WP 
    dV = NUCL_NBRC**(2.0_WP/3.0_WP) / min_nbrC
    V0 = NUCL_NBRC
    chicarb = chisoot * (36.0*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    src = -wCoeff * chicarb  * V0**(2.0_WP/3.0_WP) * mom(nMoments)

    return
  end subroutine soot_surfacereactionsmall




  !=================!
  ! SurfaceOxidiation !
  !=================!

  ! Compute Surface Oxidation source term for moment k
  ! --------------------------------------------------
  subroutine soot_surfaceoxidation(k,src)
    !use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: dV,V0,chicarb,small,large

    dV = 1.0_WP
!    dV = NUCL_NBRC**(2.0_WP/3.0_WP) / min_nbrC
    V0 = NUCL_NBRC
    chicarb = chisoot * (36.0_WP*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    small = -oxCoeff * chicarb * dV * V0**(k(1)-1.0_WP+(2.0_WP/3.0_WP)*(k(2)+1.0_WP)) * mom(nMoments)
    large = -oxCoeff * chicarb * dV * (k(1)+(2.0_WP/3.0_WP)*k(2)) * soot_fracmomlarge(k(1)-1.0_WP,k(2)+1.0_WP)
   
    src = small + large
!    src = large 
!    src = -oxCoeff * chicarb * dV * k(1) * soot_fracmom(k(1)-1.0_WP,k(2)+1.0_WP)

    return
  end subroutine soot_surfaceoxidation




  ! Compute Surface Reaction source term for moment k
  ! -------------------------------------------------
  subroutine soot_surfaceoxidationsmall(src)
    !use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: dV,V0,chicarb,small,inter,large

    dV = 1.0_WP
    V0 = NUCL_NBRC
    chicarb = chisoot * (36.0*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    ! Oxidation of the small particles
    small = -oxCoeff * chicarb * dV * V0**(-1.0_WP/3.0_WP) * mom(nMoments)

    ! Oxidation of the larger particles as they become small
    ! Proportionality constant is the ratio between the volumes of the two modes
    inter = V0 / (soot_fracmomlarge(1.0_WP,0.0_WP) / soot_fracmomlarge(0.0_WP,0.0_WP))
    large = inter * oxCoeff * chicarb * dV * soot_fracmomlarge(-1.0_WP,1.0_WP)

    src = small + large
!     src = small
    return
  end subroutine soot_surfaceoxidationsmall




  !=================!
  ! Fragmentation   !
  !=================!
  ! Compute Fragmentation source term for moment k
  ! ----------------------------------------------
  subroutine soot_fragmentation(k,src)
    !use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: chicarb

    chicarb = chisoot * (36.0*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    src = (2.0_WP**(1.0_WP-k(1)-k(2))-1.0_WP) * 2.0_WP*o2Coeff*chicarb*soot_fracmomlarge(k(1)-1.0_WP,k(2)+1.0_WP)

    return
  end subroutine soot_fragmentation




  ! Compute Fragmentation source term for delta function
  ! ----------------------------------------------------
  subroutine soot_fragmentationsmall(src)
    !use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: V0,chicarb,inter

    V0 = NUCL_NBRC
    chicarb = chisoot * (36.0*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    inter = V0 / (soot_fracmomlarge(1.0_WP,0.0_WP) / soot_fracmomlarge(0.0_WP,0.0_WP))

    src = inter * 2.0_WP*o2Coeff*chicarb*soot_fracmomlarge(-1.0_WP,1.0_WP)

    return
  end subroutine soot_fragmentationsmall

  !=== End HMOM Model=====================================================================================!

end module soot_hmom_mod




! ======================================================================== !
! Part I: Initialize the soot module                                       !
! ======================================================================== !
subroutine soot_hmom_init(nscalar_,Y_,rho_)
  use soot_hmom_mod
  implicit none
  integer, intent(in) :: nscalar_
  double precision, intent(in) :: rho_
  double precision, dimension(nscalar_) :: Y_
!  real(WP) :: min_nbrC, min_nbrH

  ! Initialize the finite rate chemistry portion of soot, if applicable
  call soot_chem_init
  call soot_chem_C_dimer(Y_,rho_,C_dimer)
  
  ! Define some parameters
  CarbonToDiam = (6.0_WP*MolarMassSoot/(Pi*SootDensity*Avogadro))**(1.0_WP/3.0_WP)
  Cfm = Avogadro * (8.0_WP*Pi*Rgas/MolarMassSoot)**(0.5_WP) &
       * (0.75_WP*MolarMassSoot/(Avogadro*Pi*SootDensity))**(2.0_WP/3.0_WP)
  ! Need square root temperature for free-molecular regime constant
  Ccn = 8.0_WP * Rgas / 3.0_WP
  ! Need temperature / viscosity for continuum regime constant
  lambda = 3.0_WP * sqrt(Pi/(8.0_WP*Rgas)) / (6.0_WP*MolarMassSoot/(Pi*sootDensity*Avogadro))**(1.0_WP/3.0_WP)
  ! Need viscosity * sqrt molar mass / (density * sqrt tempeature) for slip factor constant
  av = 1.0_WP - (2.0_WP / Df)
  as = (3.0_WP / Df) - 1.0_WP
  FitC = 2.0_WP/3.0_WP
  FitExp = -0.2043_WP

  ! Allocate arrays for HMOM
  allocate(frac(dim))
  allocate(kmom(dim))
  allocate(mom(nMoments))
  allocate(old_mom(nMoments))
  allocate(src_mom(nMoments))

  ! Define the momenst used for HMOM
  call soot_define_moments

  ! number(nbr) of C in DIMER (TWO PAH MOLECULRES)
  min_nbrC = max(1E-3,C_dimer) !20.0_WP
  min_nbrH = 16.0_WP
!  write(*,*) "C_dimer = ", C_dimer
  ! Set dimer size
  DIMER_NBRC = min_nbrC
  DIMER_NBRH = min_nbrH

  ! Set nucleated particle size
  NUCL_NBRC = 2.0_WP*min_nbrC !*W_C/SootDensity
  NUCL_SURF = NUCL_NBRC**(2.0_WP/3.0_WP)
  NUCL_NBRH = 2.0_WP*min_nbrH !*W_C/SootDensity

  return
end subroutine soot_hmom_init



! ======================================================================== !
! Part II: Soot source Routine                                             !
! -> Compute the source terms for the pressure and finite rate chemistry   !
! ======================================================================== !
subroutine soot_hmom_source(timeStep,composition,moment,density,wmixture,srcSoot)
  use soot_hmom_mod
  implicit none
  real(WP), intent(in) :: timeStep
  real(WP), intent(in), dimension(nComposition) :: composition
  real(WP), intent(in), dimension(nMoments) :: moment
  real(WP), intent(in) :: density,wmixture  ! unit ???
  real(WP), intent(out), dimension(nSootGas+nMoments) :: srcSoot

  ! Update the chemical species/temperature source terms for finite rate chemistry
  ! only function of gas species
  call soot_chem_gas_source(composition,moment,density,srcSoot(1:nSootGas)) 


  ! Updata source terms for Moment equations
  ! function of gas species and momnets
  call soot_hmom_moment_source(timeStep,composition,moment,density,wmixture,srcSoot(nSootGas+1:nSootGas+nMoments))

  return
end subroutine soot_hmom_source




subroutine soot_hmom_moment_source(dt,SC_,mom_,dens_,Wmix,srcMoments)
  use soot_hmom_mod
  implicit none
  real(WP), intent(in) :: dt
  real(WP), intent(in), dimension(nComposition) :: SC_
  real(WP), intent(in), dimension(nMoments) :: mom_
  real(WP), intent(in) :: dens_, Wmix  !unit ???
  real(WP), intent(out), dimension(nMoments) :: srcMoments
 ! For viscosity by empirical formula
  real(WP),parameter :: T0 = 294_WP
  real(WP),parameter :: a  = 0.000038275_WP
  real(WP),parameter :: b  = 1.67_WP

  integer  :: i,j,k,n
  real(WP) :: t_start,rate

  real(WP) :: dt_
  integer  :: step

  ! Evaluate all quantities at half timestep for second order
  ! Get Some quantites for global variables
  dt_ = dt
  dens = dens_
  temp = SC_(nComposition)            ! Temperature
  MassFracs = SC_(1:nComposition-1)   ! Mass fractions of species
  sqrtT = sqrt(temp)                  ! sqrt of temperature
  VISC = a * (temp/T0)**b *dens_      ! viscosity 
  T_MU = temp/ VISC
  MUsqrtW_RHOsqrtT = VISC/dens_*sqrt(Wmix/temp)

  !print *,dens,temp,MassFracs(1:45)

  !== 2.1) =============================================================!
  ! Compute quantities which are only functions of the combustion scalars
  call soot_chem_ksg(MassFracs,temp,dens,wCoeff)  ! surface growth by reaction
  call soot_chem_kox(MassFracs,temp,dens,oxCoeff) ! surface oxidation
  call soot_chem_ko2(MassFracs,temp,dens,o2Coeff) ! fracm
  !print *, 'wCoeff',wcoeff,oxcoeff,o2coeff
  call soot_chem_dimerprodrate(MassFracs,temp,dens,prodRate)
  local_nbrC = DIMER_NBRC
  ! Rescale dimer production rate for HMOM
  prodRate = local_nbrC*prodRate / DIMER_NBRC

  !== 2.2) =============================================================!
  ! Get the values
  mom = mom_*dens  
  ! Save old values/soot
  old_mom = mom

  ! clip mom initial values before compute src !!!! important for the first step
  ! mom =0
  !======================!
  mom = mom / dens
  call soot_hmom_clip
  mom = mom * dens
  !======================!

!  call soot_hmom_compute_src

!  DO n=1,nMoments
!     rate = -dt*src_mom(n)/max(1.E-60_WP,mom(n))
!     src_mom(n) = src_mom(n)/max(1.0_WP, rate) !mom - old_mom
!  end do
!  srcMoments = src_mom

!============================
  t_start = 0.0_WP
  step = 0
  do while (t_start.lt.dt .and. step.lt.10000)
     call soot_hmom_compute_src
     dt_ = dt - t_start
     rate = 3.0_WP*maxval(abs(-dt_*src_mom)/mom)
     if (rate.gt.1.0_WP) then
        dt_ = dt_/rate
     end if
     mom = mom + dt_*src_mom
     mom = mom / dens
     call soot_hmom_clip
     mom = mom * dens
     t_start = t_start+dt_
     step = step+1
  end do
  if (t_start.lt.dt) then
     call soot_hmom_compute_src
     mom = mom + (dt-t_start)*src_mom
  end if


  ! Copy back the source terms
  mom = mom / dens
  call soot_hmom_clip
  mom = mom * dens

  srcMoments = (mom - old_mom) / dt

!      if (betaC_out .gt. 2.5e3) then
!           write(*,*) betaC_out
 !          write(*,*) old_mom(1),mom(1)
!           write(*,*) old_mom(2),mom(2)
!           write(*,*) old_mom(3),mom(3)
!           write(*,*) old_mom(4),mom(4)
!           write(*,*) dens,temp,Wmix
!           do i=1,25
!              write(*,*) MassFracs(i)
!           enddo
!      end if


  return
end subroutine soot_hmom_moment_source




! ======================== !
! Clip the soot quantities !
! ======================== !
subroutine soot_hmom_clip !for soot_hmom_source
  use soot_hmom_mod
  implicit none
  integer :: n

  ! Check for globally small moments
  if ( mom(1).lt.SMALLWEIGHT              .or. &
       mom(2).lt.SMALLWEIGHT*NUCL_NBRC    .or. &
!       mom(3).lt.SMALLWEIGHT*NUCL_NBRC**2 .or. &
!       mom(4).lt.SMALLWEIGHT*NUCL_NBRC**3 ) then
       mom(3).lt.SMALLWEIGHT*NUCL_SURF     ) then
!       mom(6).lt.SMALLWEIGHT*NUCL_NBRC*NUCL_NBRC ) then
     mom(2) = max(mom(2), SMALLWEIGHT*NUCL_NBRC)
!     mom(1) = mom(2) / NUCL_NBRC
     mom(1) = max(mom(1),SMALLWEIGHT)
     mom(3) = max(mom(3),SMALLWEIGHT*NUCL_SURF) !    NUCL_SURF * mom(1)
!     mom(4) = max(mom(4),SMALLWEIGHT*NUCL_NBRC**3)
    ! mom(5) = mom(1)*NUCL_SURF
    ! mom(6) = max(mom(6),SMALLWEIGHT*NUCL_NBRC*NUCL_NBRC)
!     mom(nMoments) = mom(1)
  end if

  ! Check for size of second mode
!  if (mom(2).lt.NUCL_NBRC*mom(1) .or. &
!      mom(3).lt.NUCL_SURF*mom(1)) then
!     mom(1) = mom(2) / NUCL_NBRC
!     mom(3) = mom(1) * NUCL_SURF
!  end if
  
  if (mom(2).lt.NUCL_NBRC*mom(1)) then
     mom(1) = mom(2) / NUCL_NBRC
  end if
  if (mom(3).lt.NUCL_SURF*mom(1)) then
     mom(3) = mom(1) * NUCL_SURF
  end if

  ! Check for (co)variance of second mode
!  if (nMoments.ge.7) then
!     mom(3) = max(mom(3),mom(2)*mom(2)/mom(1))
!     mom(4) = max(mom(4),mom(2)*mom(2)*mom(2)/mom(1))
!     mom(6) = max(mom(6),mom(2)*mom(5)/mom(1))
!  end if

  ! Check for small weight of first mode
  if (mom(nMoments).lt.SMALLWEIGHT) then
     do n=1,nMoments-1
        mom(n) = mom(n) + (SMALLWEIGHT-mom(nMoments))*NUCL_NBRC**moments(n,1)*NUCL_SURF**moments(n,2)
     end do
     mom(nMoments) = SMALLWEIGHT
  end if
  if (mom(nMoments).gt.mom(1)) then
      mom(nMoments) = 0.999_WP * mom(1)
  end if

  ! Clip PAH mass fraction
  !if (use_pah) then
  !   if (mom(nEquations).lt.SMALLY) then
  !      mom(nEquations) = SMALLY
  !   end if
  !end if

  return
end subroutine soot_hmom_clip




subroutine soot_hmom_rhodot(rhodot,Y_,temp_,soot_,dens_)
  use soot_hmom_mod
  implicit none

  real(WP), intent(out) :: rhodot
  real(WP), dimension(nComposition-1), intent(in) :: Y_
  real(WP), dimension(nMoments), intent(in) :: soot_
  real(WP), intent(in) :: temp_,dens_

  rhodot = 0.0_WP
  call soot_chem_rhodot(rhodot,Y_,temp_,soot_,dens_)

  return
end subroutine soot_hmom_rhodot




subroutine soot_hmom_poststep(output,SC_,mom_,dens_,Wmix)
  use soot_hmom_mod
  implicit none

  real(WP), dimension(4+10), intent(out) :: output ! numdens,fv,partdiam,partaggr
  real(WP), dimension(nMoments), intent(in) :: mom_
  real(WP), intent(in), dimension(nComposition) :: SC_
  real(WP), intent(in) :: dens_, Wmix  !unit ???
 ! For viscosity by empirical formula
  real(WP),parameter :: T0 = 294_WP
  real(WP),parameter :: a  = 0.000038275_WP
  real(WP),parameter :: b  = 1.67_WP
  real(WP):: src_tmp

  output = 0.0_WP

 ! Evaluate all quantities at half timestep for second order
  ! Get Some quantites for global variables
  dens = dens_
  temp = SC_(nComposition)            ! Temperature
  MassFracs = SC_(1:nComposition-1)   ! Mass fractions of species
  sqrtT = sqrt(temp)                  ! sqrt of temperature
  VISC = a * (temp/T0)**b *dens_      ! viscosity 
  T_MU = temp/ VISC
  MUsqrtW_RHOsqrtT = VISC/dens_*sqrt(Wmix/temp)

  !print *,dens,temp,MassFracs(1:45)

  !== 2.1) =============================================================!
  ! Compute quantities which are only functions of the combustion scalars
  call soot_chem_ksg(MassFracs,temp,dens,wCoeff)  ! surface growth by reaction
  call soot_chem_kox(MassFracs,temp,dens,oxCoeff) ! surface oxidation
  call soot_chem_ko2(MassFracs,temp,dens,o2Coeff) ! fracm
  !print *, 'wCoeff',wcoeff,oxcoeff,o2coeff
  call soot_chem_dimerprodrate(MassFracs,temp,dens,prodRate)
  local_nbrC = DIMER_NBRC
  ! Rescale dimer production rate for HMOM
  prodRate = local_nbrC*prodRate / DIMER_NBRC

  mom = mom_
  call soot_hmom_clip
  mom = mom * dens_

  output(1) = soot_fracmom(0.0_WP,0.0_WP) * Avogadro*1.0e-6_WP
  output(2) = soot_fracmom(1.0_WP,0.0_WP) * MolarMassSoot/SootDensity
!  output(3) = soot_fracmom(1.0_WP,-1.0_WP) / soot_fracmom(0.0_WP,0.0_WP) * CarbonToDiam*1.0E9_WP
  output(3) = (soot_fracmom(2.0_WP,-2.0_WP)/soot_fracmom(1.0_WP,-1.0_WP)) * CarbonToDiam*1.0E9_WP
  output(4) = soot_fracmom(-2.0_WP,3.0_WP) / soot_fracmom(0.0_WP,0.0_WP)

  if (output(1)<0.0_WP) output(1) = 0.0_WP
  if (output(2)<0.0_WP) output(2) = 0.0_WP
  if (output(3)<0.0_WP.or.output(3).ge.1000.0_WP) output(3) = 0.0_WP
  if (output(4)<0.0_WP) output(4) = 0.0_WP


  ! Prepare for compute source terms
  call soot_pah_molecules

  ! Number density source terms
  kmom = moments(1,:)
  call soot_nucleation(kmom,src_tmp)
  output(5) = src_tmp !* Avogadro*1.0e-6_WP
  call soot_coagulation(kmom,src_tmp)
  output(6) = src_tmp !* Avogadro*1.0e-6_WP
  call soot_condensation(kmom,src_tmp)
  output(7) = src_tmp !* Avogadro*1.0e-6_WP
  call soot_surfacereaction(kmom,src_tmp)
  output(8) = src_tmp !* Avogadro*1.0e-6_WP
  call soot_surfaceoxidation(kmom,src_tmp)
  output(9) = src_tmp !* Avogadro*1.0e-6_WP
  
  ! Volume fraction source terms
  kmom = moments(2,:)
  call soot_nucleation(kmom,src_tmp)
!  call soot_nucleationsmall(src_tmp)
  output(10) = src_tmp !* MolarMassSoot !/SootDensity
  call soot_coagulation(kmom,src_tmp)
!  call soot_coagulationsmall(src_tmp)
  output(11) = src_tmp  !* MolarMassSoot !/SootDensity
  call soot_condensation(kmom,src_tmp)
!  call soot_condensationsmall(src_tmp)
  output(12) = src_tmp  !* MolarMassSoot !/SootDensity
  call soot_surfacereaction(kmom,src_tmp)
!  call soot_surfacereactionsmall(src_tmp)
  output(13) = src_tmp  !* MolarMassSoot !/SootDensity
  call soot_surfaceoxidation(kmom,src_tmp)
!  call soot_surfaceoxidationsmall(src_tmp)
  output(14) = src_tmp  !* MolarMassSoot !/SootDensity

  return
end subroutine soot_hmom_poststep

subroutine soot_hmom_destroy
   use soot_hmom_mod
   implicit none

   deallocate(moments)
   deallocate(frac)
   deallocate(kmom)
   deallocate(mom)
   deallocate(old_mom)
   deallocate(src_mom)

   return
end subroutine soot_hmom_destroy

subroutine soot_hmom_other(other_,otherArr_)
  use soot_hmom_mod
  implicit none
  integer, intent(in) :: other_
  double precision, dimension(other_), intent(out) :: otherArr_
  real(WP) :: C_SootStar

  call soot_chem_ksg(MassFracs,temp,dens,wCoeff)  ! surface growth by reaction
  call soot_chem_kox(MassFracs,temp,dens,oxCoeff) ! surface oxidation
  call soot_chem_sootstar(MassFracs,temp,dens,C_SootStar)

  otherArr_(1) = wCoeff
  otherArr_(2) = oxCoeff
  otherArr_(3) = DIMER_conc !(sqrt(soot_beta_cond()**2 + 4.0_WP*soot_beta_nucl()*prodRate)-soot_beta_cond())/(2.0_WP*soot_beta_nucl()) ! DIMER_conc
  otherArr_(4) = C_SootStar
  otherArr_(5) = betaN_out !soot_beta_nucl()
  otherArr_(6) = betaC_out !soot_beta_cond()
  otherArr_(7) = prodRate_out
  otherArr_(8) = aromaticConc
  otherArr_(9) = C_dimer

end subroutine soot_hmom_other



