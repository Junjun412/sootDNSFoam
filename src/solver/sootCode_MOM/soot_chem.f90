!---------------------------------------------------------!
! Modified according to NGAsoot code for fdfsolver        !
! Wang Han Jul. 30 2016, advisored by Venkat Raman        !
!---------------------------------------------------------!
module soot_chem_mod ! harding code for reduced Ethylene Mech
  use soot_mod
  implicit none

contains

  ! call by HMOM and SEMI
  ! subroutine soot_chem_init (use_hmom==1?)
  ! subroutine soot_chem_gas_source
  !             |- call subroutine soot_chem_hmom_gas_source
  !             |- call subroutine soot_chem_semi_gas_source
  ! subroutine soot_chem_rhodot (use_hmom==1)

  ! only call by HMOM in fdf_soot_hmom.f90
  ! subroutine soot_chem_ksg
  !             |-call soot_sootstar
  ! subroutine soot_chem_kox
  !             |-call soot_sootstar
  ! subroutine soot_chem_ko2
  !             |-call soot_sootstar
  ! subroutine soot_chem_dimerprodrate

  ! only call by SEMI in fdf_soot_semi.f90
  ! subroutine soot_chem_cC2H2
  ! subroutine soot_chem_cO2

  !--------------------- SEMI -------------------------!
  subroutine soot_chem_semi_gas_source(Y_,temp_,moment_,dens_,src_) !subroutine soot_finitechem_pahprodrates
    implicit none

    real(WP), intent(in), dimension(nComposition-1) :: Y_
    real(WP), intent(in), dimension(nMoments) :: moment_
    real(WP), intent(in) :: temp_, dens_
    real(WP), intent(out), dimension(nSootGas-1) :: src_
    real(WP) :: A_s,ka,kb,kz,kT,x

!    ka = 20.0_WP * exp(-15098.0_WP/temp_)
!    kb = 4.46e-3 * exp(-7650.0_WP/temp_)
!    kz = 21.3_WP * exp(2063.0_WP/temp_)
!    kT = 1.51e+5 * exp(-48817.0_WP/temp_)
!    x  = 1.0_WP/(1+kT/(kb*(dens_*Y_(i_O2)/W_GAS(i_O2))))
    
    src_ = 0.0_WP
    !c2h2,h2,o2,co
    ! Nucleation
    src_(1) = src_(1) - W_GAS(i_C2H2) * (A_NUCL*exp(-E_NUCL/temp_)) * (dens_*Y_(i_C2H2)/W_GAS(i_C2H2))
    src_(2) = src_(2) + W_GAS(i_H2) * (A_NUCL*exp(-E_NUCL/temp_)) * (dens_*Y_(i_C2H2)/W_GAS(i_C2H2))

    ! Surface Growth
    A_s = Pi*(6.0_WP*dens_*moment_(2)/Pi/SootDensity)**(2.0_WP/3.0_WP)*(dens_*moment_(1))**(1.0_WP/3.0_WP)
    src_(1) = src_(1) - W_GAS(i_C2H2) * (A_GROWTH*exp(-E_GROWTH/temp_)) * sqrt(A_s) * (dens_*Y_(i_C2H2)/W_GAS(i_C2H2))
    src_(2) = src_(2) + W_GAS(i_H2)   * (A_GROWTH*exp(-E_GROWTH/temp_)) * sqrt(A_s) * (dens_*Y_(i_C2H2)/W_GAS(i_C2H2))
    src_(3) = src_(3) - 0.5_WP*W_GAS(i_O2) * (A_OXI*sqrt(temp_)*exp(-E_OXI/temp_)) * A_s  *  (dens_*Y_(i_O2)/W_GAS(i_O2))
    src_(4) = src_(4) + W_GAS(i_CO) * (A_OXI*sqrt(temp_)*exp(-E_OXI/temp_)) * A_s * (dens_*Y_(i_O2)/W_GAS(i_O2))
!    src_(3) = src_(3) - 0.5_WP*W_GAS(i_O2) * 120.0_WP * (ka*(dens_*Y_(i_O2)/W_GAS(i_O2))*x/(1+kz*(dens_*Y_(i_O2)/    &
!              W_GAS(i_O2))) + kb*(dens_*Y_(i_O2)/W_GAS(i_O2))*(1-x))
!    src_(4) = src_(4) + W_GAS(i_CO) * 120.0_WP * (ka*(dens_*Y_(i_O2)/W_GAS(i_O2))*x/(1+kz*(dens_*Y_(i_O2)/W_GAS(i_O2)))   &
!              + kb*(dens_*Y_(i_O2)/W_GAS(i_O2))*(1-x))
    src_(4) = src_(4) + W_GAS(i_CO) * phi_OH * K_OXI_OH / sqrt(temp_) * A_s * (dens_*Y_(i_OH)/W_GAS(i_OH))
    src_(4) = src_(4) + W_GAS(i_CO) * phi_O  * K_OXI_O  / sqrt(temp_) * A_s * (dens_*Y_(i_O)/W_GAS(i_O))
    src_(5) = src_(5) - W_GAS(i_OH) * phi_OH * K_OXI_OH / sqrt(temp_) * A_s * (dens_*Y_(i_OH)/W_GAS(i_OH)) !i_OH
    src_(6) = src_(6) - W_GAS(i_O)  * phi_O  * K_OXI_O  / sqrt(temp_) * A_s * (dens_*Y_(i_O)/W_GAS(i_O)) !i_O
    src_(7) = src_(7) + W_GAS(i_H)  * phi_OH * K_OXI_OH / sqrt(temp_) * A_s * (dens_*Y_(i_OH)/W_GAS(i_OH)) !i_H
    return
  end subroutine soot_chem_semi_gas_source


  subroutine soot_chem_cC2H2(Y_,dens_,C_C2H2_)
    implicit none

    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(in) :: dens_
    real(WP), intent(out) :: C_C2H2_

    C_C2H2_ = dens_ * Y_(i_C2H2) / W_GAS(i_C2H2)

    return
  end subroutine soot_chem_cC2H2



  subroutine soot_chem_cO2(Y_,dens_,C_O2_)
    implicit none

    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(in) :: dens_
    real(WP), intent(out) :: C_O2_

    C_O2_ = dens_ * Y_(i_O2) / W_GAS(i_O2)

    return
  end subroutine soot_chem_cO2

  subroutine soot_chem_cOH(Y_,dens_,C_OH_)
    implicit none

    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(in) :: dens_
    real(WP), intent(out) :: C_OH_

    C_OH_ = dens_ * Y_(i_OH) / W_GAS(i_OH)

    return
  end subroutine soot_chem_cOH

  subroutine soot_chem_cO(Y_,dens_,C_O_)
    implicit none

    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(in) :: dens_
    real(WP), intent(out) :: C_O_

    C_O_ = dens_ * Y_(i_O) / W_GAS(i_O)

    return
  end subroutine soot_chem_cO

  !--------------------- MOM -------------------------!
  subroutine soot_chem_mom_gas_source(SC_,dens_,srcPAH)
    implicit none

    real(WP), intent(in), dimension(nComposition) :: SC_
    real(WP), intent(in) :: dens_  !unit ???
    real(WP), intent(out), dimension(nSootGas) :: srcPAH
    integer :: i,j,k,n
    real(WP) :: chicarb,reacRate,oxiRate_O2,oxiRate_OH,dV

    srcPAH = 0.0D0
!    dV = 1.0_WP
    dV =  NUCL_NBRC**(2.0_WP/3.0_WP) / min_nbrC

    call soot_pah_prodrates(SC_(1:nComposition-1),SC_(nComposition),dens_,pahprodrate)
    call soot_chem_ksg(SC_(1:nComposition-1),SC_(nComposition),dens_,reacRate)  ! surface growth by reaction
    call soot_chem_kox1(SC_(1:nComposition-1),SC_(nComposition),dens_,oxiRate_O2)
    call soot_chem_kox2(SC_(1:nComposition-1),SC_(nComposition),dens_,oxiRate_OH)
    chicarb = chisoot * (36.0_WP*Pi)**(1.0_WP/3.0_WP)*(MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    srcPAH(nPAH+1)=min(0.0,srcPAH(nPAH+1)-W_GAS(i_C2H2)*dV*reacRate*chicarb*mom(3)*MolarMassSoot) !*MolarMassSoot/SootDensity
    srcPAH(nPAH+2)=min(0.0,srcPAH(nPAH+2)-W_GAS(i_O2)*dV*oxiRate_O2*chicarb*mom(3)*MolarMassSoot)
    srcPAH(nPAH+3)=min(0.0,srcPAH(nPAH+3)-W_GAS(i_OH)*dV*oxiRate_OH*chicarb*mom(3)*MolarMassSoot)
    srcPAH(nPAH+4)=max(0.0,srcPAH(nPAH+4)+2.0_WP*W_GAS(i_CO)*dV*oxiRate_O2*chicarb*mom(3)*MolarMassSoot)
    srcPAH(nPAH+4)=max(0.0,srcPAH(nPAH+4)+W_GAS(i_CO)*dV*oxiRate_OH*chicarb*mom(3)*MolarMassSoot)
    srcPAH(nPAH+5)=max(0.0,srcPAH(nPAH+5)+W_GAS(i_H) *dV*oxiRate_OH*chicarb*mom(3)*MolarMassSoot)
    srcPAH(nPAH+5)=max(0.0,srcPAH(nPAH+5)+W_GAS(i_H) *dV*reacRate*chicarb*mom(3)*MolarMassSoot)

    do n=1,nPAH
       if (i_PAH(n).ne.-1) then
       srcPAH(n) = srcPAH(n) + pahprodrate(n)
       srcPAH(nPAH+6) = srcPAH(nPAH+6) + pahprodrate(n)
       end if
    end do

!     srcPAH(nPAH+6) = srcPAH(nPAH+6)+srcPAH(nPAH+1)+srcPAH(nPAH+2)+srcPAH(nPAH+3)+srcPAH(nPAH+4)+srcPAH(nPAH+5)

    return
  end subroutine soot_chem_mom_gas_source

  subroutine soot_chem_C_dimer(SC_,dens_,C_dimer_)
    implicit none

    double precision, intent(in) :: dens_
    double precision, intent(in), dimension(nComposition) :: SC_
    real(WP), intent(out) :: C_dimer_
    real(WP) :: C_pah_total, C_total
    integer :: n

    call soot_pah_prodrates(SC_(1:nComposition-1),SC_(nComposition),dens_,pahprodrate)
    C_pah_total = -1.0E-10
    C_total = -1.0E-10
    do n=1,nPAH
        C_pah_total = C_pah_total + pahprodrate(n)*size_PAH(n)
        C_total = C_total + pahprodrate(n)
    end do
    C_dimer_ = 2.0 * C_pah_total / C_total

    end subroutine soot_chem_C_dimer

  subroutine soot_pah_prodrates(Y_,temp_,dens_,pahprodrate_) !subroutine soot_finitechem_pahprodrates
    implicit none

    real(WP), intent(in), dimension(nComposition-1) :: Y_
    real(WP), intent(in) :: temp_, dens_
    real(WP), intent(out), dimension(nPAH) :: pahprodrate_
    real(WP) :: aromconc_i,aromconc_j,betadimer_ij
    real(WP) :: vol_i, vol_j, diam_i, diam_j, mu_ij
    integer :: i,j,j1
    real(WP) :: DA, sd_lin, KB, pwd_lin, sticking_coeff_T_lin, dPAH_i, dPAH_j

    DA = 1.395_WP*sqrt(3.0)
    sd_lin=(2.07e-13*temp_-2.14e-11)*1e10
    KB = Rgas/AVOGADRO

    do i=1,nPAH
       if (i_PAH(i).ne.-1) then
       pahprodrate_(i) = 0.0
       do j=1,nPAH
          aromconc_i = dens_ * Y_(i_PAH(i)) /(W_GAS(i_PAH(i)))
          aromconc_j = dens_ * Y_(i_PAH(j)) /(W_GAS(i_PAH(j)))
          vol_i  = MolarMassSoot * size_PAH(i) / (AVOGADRO*SootDensity)
          vol_j  = MolarMassSoot * size_PAH(j) / (AVOGADRO*SootDensity)
          diam_i = 2.0_WP * (3.0_WP*vol_i/(4.0_WP*Pi))**(1.0_WP/3.0_WP)
          diam_j = 2.0_WP * (3.0_WP*vol_j/(4.0_WP*Pi))**(1.0_WP/3.0_WP)
          mu_ij  = vol_i * vol_j * sootDensity / (vol_i + vol_j)
          betadimer_ij = sqrt(8.0_WP*Pi*Rgas*temp_/(mu_ij*AVOGADRO))* &
                         0.25_WP * (diam_i+diam_j)**2.0_WP * AVOGADRO

          dPAH_i=DA*sqrt(2.0*size_PAH(i)/3.0)
          dPAH_j=DA*sqrt(2.0*size_PAH(j)/3.0)
          pwd_lin=5.0e-20/6.0*(1.0/2.0*dPAH_i*dPAH_j/((dPAH_i+dPAH_j+sd_lin) &
                  *sd_lin)+1.0/2.0*dPAH_i*dPAH_j/((dPAH_i+sd_lin)*(dPAH_j+sd_lin)) &
                  +log((dPAH_i+dPAH_j+sd_lin)*sd_lin/(dPAH_i+sd_lin)/(dPAH_j+sd_lin)))
          sticking_coeff_T_lin=1.0-((1.0+pwd_lin/KB/temp_)*exp(-pwd_lin/KB/temp_))

!          sticking_coeff_T_lin = 0.01 * sticking_coeff_T_lin ! give a factor

          if (i==j) then
             pahprodrate_(i) = pahprodrate_(i) - (W_GAS(i_PAH(i))) *0.5_WP  &
                               * 2.0_WP*sticking_coeff_T_lin *    &
                               betadimer_ij * aromconc_i * aromconc_j
          else
             pahprodrate_(i) = pahprodrate_(i) - (W_GAS(i_PAH(i))) * 0.5_WP &
                               * sticking_coeff_T_lin * betadimer_ij  &
                               * aromconc_i * aromconc_j
          end if
      end do
      end if
    end do
    return
  end subroutine soot_pah_prodrates


!! wCoeff                            
  subroutine soot_chem_ksg(Y_,temp_,dens_,ksg)
    implicit none

    real(WP), intent(in) :: temp_,dens_
    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(out) :: ksg
    real(WP) :: C_C2H2_, C_SootStar

    C_C2H2_ = dens_ * Y_(i_C2H2) / W_GAS(i_C2H2)

    call soot_chem_sootstar(Y_,temp_,dens_,C_SootStar)

    ksg = A4_s * temp_**n4_s * exp(-E4_s/temp_) * C_C2H2_ * C_SootStar

    return
  end subroutine soot_chem_ksg



  subroutine soot_chem_kox(Y_,temp_,dens_,kox)
    implicit none

    real(WP), intent(in) :: temp_,dens_
    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(out) :: kox
    real(WP) :: k6, C_OH_, C_O2_, C_SootStar

    C_OH_ = dens_ * Y_(i_OH) / W_GAS(i_OH)
    C_O2_ = dens_ * Y_(i_O2) / W_GAS(i_O2)

    call soot_chem_sootstar(Y_,temp_,dens_,C_SootStar)

    k6 = 8.94_WP * eff6 * sqrt(temp_) * AVOGADRO

    kox = A5_s * temp_**n5_s * exp(-E5_s/temp_) * C_O2_ * C_SootStar + &
         (0.5_WP/chisoot) * k6 * C_OH_

!    kox = A5_s * temp_**n5_s * exp(-E5_s/temp_) * C_O2_ * 2.0_WP * C_SootStar + &
!         (0.5_WP/chisoot) * k6 * C_OH_
!    kox = A5_s * temp_**n5_s * exp(-E5_s/temp_) * C_O2_ * C_SootStar + &

    return
  end subroutine soot_chem_kox



  subroutine soot_chem_kox1(Y_,temp_,dens_,kox)
    implicit none

    real(WP), intent(in) :: temp_,dens_
    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(out) :: kox
    real(WP) :: k6, C_OH_, C_O2_, C_SootStar

!    C_OH_ = dens_ * Y_(i_OH) / W_GAS(i_OH)
    C_O2_ = dens_ * Y_(i_O2) / W_GAS(i_O2)

    call soot_chem_sootstar(Y_,temp_,dens_,C_SootStar)

!    k6 = 8.94_WP * eff6 * sqrt(temp_) * AVOGADRO
    kox = A5_s * temp_**n5_s * exp(-E5_s/temp_) * C_O2_ * C_SootStar 

    return
  end subroutine soot_chem_kox1



  subroutine soot_chem_kox2(Y_,temp_,dens_,kox)
    implicit none

    real(WP), intent(in) :: temp_,dens_
    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(out) :: kox
    real(WP) :: k6, C_OH_, C_O2_, C_SootStar

    C_OH_ = dens_ * Y_(i_OH) / W_GAS(i_OH)
!    C_O2_ = dens_ * Y_(i_O2) / W_GAS(i_O2)

!    call soot_chem_sootstar(Y_,temp_,dens_,C_SootStar)

    k6 = 8.94_WP * eff6 * sqrt(temp_) * AVOGADRO
    kox = (0.5_WP/chisoot) * k6 * C_OH_

    return
  end subroutine soot_chem_kox2



  subroutine soot_chem_ko2(Y_,temp_,dens_,ko2)
    implicit none

    real(WP), intent(in) :: temp_,dens_
    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(out) :: ko2
    real(WP) :: C_O2_, C_SootStar

    C_O2_ = dens_ * Y_(i_O2) / W_GAS(i_O2)

!    call soot_chem_sootstar(Y_,temp_,dens_,C_SootStar)

!    ko2 = A5_s * temp_**n5_s * exp(-E5_s/temp_) * C_O2_ * C_SootStar

    ko2 = A5_s * temp_**n5_s * exp(-E5_s/temp_) * C_O2_

    return
  end subroutine soot_chem_ko2



  subroutine soot_chem_dimerprodrate(Y_,temp_,dens_,prodRate_)
    implicit none
    real(WP), dimension(nComposition-1), intent(in) :: Y_ ! N_tot
    real(WP), intent(in) :: temp_,dens_
    real(WP), intent(out) :: prodRate_
    integer :: i, j, j1
    real(WP) :: aromconc_i,aromconc_j,betadimer_ij
    real(WP) :: vol_i, vol_j, diam_i, diam_j, mu_ij
    real(WP) :: DA, sd_lin, KB, pwd_lin, sticking_coeff_T_lin, dPAH_i, dPAH_j

    prodRate_ = 0.0_WP
    aromaticConc = 0.0_WP
    DA = 1.395_WP*sqrt(3.0)
    sd_lin=(2.07e-13*temp_-2.14e-11)*1e10
    KB = Rgas/AVOGADRO

    do i=1,nPAH
       if (i_PAH(i).ne.-1) then
          aromaticConc = aromaticConc + dens_ * Y_(i_PAH(i)) / W_GAS(i_PAH(i))
          j1 = i
          do j=j1,nPAH
             aromconc_i = dens_ * Y_(i_PAH(i)) / W_GAS(i_PAH(i))
             aromconc_j = dens_ * Y_(i_PAH(j)) / W_GAS(i_PAH(j))
             vol_i  = MolarMassSoot * size_PAH(i) / (AVOGADRO*SootDensity)
             vol_j  = MolarMassSoot * size_PAH(j) / (AVOGADRO*SootDensity)
             diam_i = 2.0_WP * (3.0_WP*vol_i/(4.0_WP*Pi))**(1.0_WP/3.0_WP)
             diam_j = 2.0_WP * (3.0_WP*vol_j/(4.0_WP*Pi))**(1.0_WP/3.0_WP)
             mu_ij  = vol_i * vol_j * sootDensity / (vol_i + vol_j)
             betadimer_ij = sqrt(8.0_WP*Pi*Rgas*temp_/(mu_ij*AVOGADRO))*  &
                            0.25_WP * (diam_i+diam_j)**2.0_WP * AVOGADRO

             dPAH_i=DA*sqrt(2.0*size_PAH(i)/3.0)
             dPAH_j=DA*sqrt(2.0*size_PAH(j)/3.0)
             pwd_lin=5.0e-20/6.0*(1.0/2.0*dPAH_i*dPAH_j/((dPAH_i+dPAH_j+sd_lin) &
                     *sd_lin)+1.0/2.0*dPAH_i*dPAH_j/((dPAH_i+sd_lin)*(dPAH_j+sd_lin)) &
                     +log((dPAH_i+dPAH_j+sd_lin)*sd_lin/(dPAH_i+sd_lin)/(dPAH_j+sd_lin)))
             sticking_coeff_T_lin=1.0-((1.0+pwd_lin/KB/temp_)*exp(-pwd_lin/KB/temp_))
 !            sticking_coeff_T_lin = 0.01 * sticking_coeff_T_lin ! give a factor

             prodRate_ = prodRate_ + 0.5_WP * sticking_coeff_T_lin*   &
                         betadimer_ij * aromconc_i * aromconc_j
          end do
       end if
    end do
    return
  end subroutine soot_chem_dimerprodrate



  !--- For ksg, kox, ko2
  subroutine soot_chem_sootstar(Y_,temp_,dens_,C_SootStar)
    implicit none

    real(WP), dimension(nComposition-1), intent(in) :: Y_
    real(WP), intent(in) :: temp_,dens_
    real(WP), intent(out) :: C_SootStar
    real(WP) :: C_OH_,C_H_,C_H2O_,C_H2_,C_C2H2_
    real(WP) :: k1f,k2f,k3f,k1b,k2b,k3b,k4

    C_OH_ = dens_ * Y_(i_OH) / W_GAS(i_OH)  !mol/m3
    C_H_ = dens_ * Y_(i_H) / W_GAS(i_H)
    C_H2O_ = dens_ * Y_(i_H2O) / W_GAS(i_H2O)
    C_H2_ = dens_ * Y_(i_H2) / W_GAS(i_H2)
    C_C2H2_ = dens_ * Y_(i_C2H2) / W_GAS(i_C2H2)

    k1f = A1f_s * temp_**n1f_s * exp(-E1f_s/temp_)
    k2f = A2f_s * temp_**n2f_s * exp(-E2f_s/temp_)
    k3f = A3f_s * temp_**n3f_s * exp(-E3f_s/temp_)
    k1b = A1b_s * temp_**n1b_s * exp(-E1b_s/temp_)
    k2b = A2b_s * temp_**n2b_s * exp(-E2b_s/temp_)
    k3b = A3b_s * temp_**n3b_s * exp(-E3b_s/temp_)
    k4 = A4_s * temp_**n4_s * exp(-E4_s/temp_)

    if ((C_H2O_+C_H2_+C_H_+C_C2H2_).lt.1e-12_WP) then  !!!!!!!!!!!!!!!!!!!
       C_SootStar = 0.0_WP    ! Avoid infinity or -infinity
    else
       C_SootStar = (k1f*C_OH_ + k2f*C_H_ + k3f) / (k1b*C_H2O_ + k2b*C_H2_ + k3b*C_H_ + k4*C_C2H2_)
       C_SootStar = C_SootStar / (1.0_WP + C_SootStar)
    end if

    return
  end subroutine soot_chem_sootstar

end module soot_chem_mod




!------------------- MOM/SEMI -----------------------!
subroutine soot_chem_init(nscalar_,Y_)  ! harding code for Adelaide sooting flame
!   use finitechem
   use soot_chem_mod
   implicit none
   integer, intent(in) :: nscalar_
   double precision, dimension(nscalar_) :: Y_
   real(WP), dimension(:), pointer :: SC
   character(len=str_short), dimension(:), pointer :: SC_name  !str_short=8 is too small for PAH name

   integer :: i
   
   allocate(SC(nscalar_))
   allocate(SC_name(nscalar_))

   SC = Y_

   if (use_mom/=0) then
      i_PAH = -1
      !Find the PAH indices and surf/ox indices
               i_PAH(1) = 1
               i_PAH(2) = 2
               i_PAH(3) = 3
               i_PAH(4) = 4
!               i_PAH(5) = 5

               i_OH =  nPAH+1
               i_H2O = nPAH+2
               i_H =   nPAH+3
               i_H2 =  nPAH+4
               i_C2H2 =nPAH+5
               i_O2 =  nPAH+6
               i_CO =  nPAH+7
!!        end select
!      end do
      ! PAH carbon number
      size_PAH(1) = 10.0_WP
      size_PAH(2) = 12.0_WP
      size_PAH(3) = 14.0_WP
      size_PAH(4) = 16.0_WP
!      size_PAH(5) = 16.0_WP

      ! PAH collision efficiencies used in Fuel paper
!      gamma_PAH(1) = 1.260e-3_WP
!      gamma_PAH(2) = 2.186e-3_WP
!      gamma_PAH(3) = 5.443e-3_WP
!      gamma_PAH(4) = 7.801e-3_WP

      ! PAH collision efficiencies from Raj's paper
      gamma_PAH(1) = 2.013e-3_WP
      gamma_PAH(2) = 3.018e-3_WP
      gamma_PAH(3) = 5.856e-3_WP
      gamma_PAH(4) = 8.077e-3_WP
!      gamma_PAH(5) = 8.077e-3_WP
   else   ! For semi
     ! Find the indices for C2H2/O2
 !     do i=1,nscalar_-1
 !        select case(trim(SC_name(i)))
 !           case('Y_C2H2')
               i_C2H2 = 1 !C2H2
 !           case('Y_H2')
               i_H2 = 2   ! H2
 !           case('Y_O2')
               i_O2 = 3   ! O2
 !           case('Y_CO')
               i_CO = 4   ! CO
               i_OH = 5   ! OH
               i_O  = 6   ! O
               i_H  = 7   ! H
 !        end select
 !     end do
   end if

!   allocate(W_GAS(nComposition-1))
   allocate(MassFracs(nComposition-1))

   W_GAS = 0.0_WP 
   if (use_mom/=0) then
       W_GAS(i_PAH(1)) = 10*W_C+8*W_H
       if (i_PAH(2)/=-1) then
          W_GAS(i_PAH(2))  = 12*W_C+8*W_H
          W_GAS(i_PAH(3))  = 14*W_C+10*W_H
          W_GAS(i_PAH(4))  = 16*W_C+10*W_H
!          W_GAS(i_PAH(5))  = 16*W_C+10*W_H
       end if
       W_GAS(i_H2O) = 2*W_H+W_O !H2O
   end if
   W_GAS(i_H2) = 2*W_H     !H2
   W_GAS(i_C2H2) = 2*W_C+2*W_H !C2H2 
   W_GAS(i_O2) = 2*W_O       !O2
   W_GAS(i_CO) = W_C+W_O     !CO
   W_GAS(i_OH) = W_O+W_H !OH
   W_GAS(i_O)  = W_O
   W_GAS(i_H)  = W_H

   return
end subroutine soot_chem_init



subroutine soot_chem_gas_source(SC_,moment_,dens_,srcGAS)
  use soot_chem_mod
  implicit none

  real(WP), intent(in), dimension(nComposition) :: SC_
  real(WP), intent(in), dimension(nMoments) :: moment_
  real(WP), intent(in) :: dens_ 
  real(WP), intent(out), dimension(nSootGas) :: srcGAS
  real(WP) :: rhodot

  srcGAS = 0.0D0

  ! soot-related gas-composition
  if (use_mom/=0) then
     call soot_chem_mom_gas_source(SC_(1:nComposition),dens_,srcGAS(1:nSootGas-1))
  else
     call soot_chem_semi_gas_source(SC_(1:nComposition-1),SC_(nComposition),moment_,dens_,srcGAS(1:nSootGas-1))
  end if

  
  ! soot-related temperature
  call soot_chem_rhodot(rhodot,SC_(1:nComposition-1),SC_(nComposition),moment_,dens_)
  srcGAS(nSootGas) = srcGAS(nSootGas) + SC_(nComposition) * rhodot

  return
end subroutine soot_chem_gas_source



subroutine soot_chem_rhodot(rhodot,Y_,temp_,moment_,dens_)
  use soot_chem_mod
  implicit none

  real(WP), intent(out) :: rhodot
  real(WP), dimension(nComposition-1), intent(in) :: Y_
  real(WP), dimension(nMoments), intent(in) :: moment_
  real(WP), intent(in) :: temp_,dens_
  integer :: i
  real(WP) :: aromconc,betadimer,A_s,ka,kb,kz,kT,x

!    ka = 20.0_WP * exp(-15098.0_WP/temp_)
!    kb = 4.46e-3 * exp(-7650.0_WP/temp_)
!    kz = 21.3_WP * exp(2063.0_WP/temp_)
!    kT = 1.51e+5 * exp(-48817.0_WP/temp_)
!    x  = 1.0_WP/(1+kT/(kb*(dens_*Y_(i_O2)/W_GAS(i_O2))))

  rhodot = 0.0_WP
  if (.not.soot_rho) return
   if (use_mom/=0) then
      do i=1,nPAH
         if (i_PAH(i).ne.-1) then
         aromconc  = dens_ * Y_(i_PAH(i)) / W_GAS(i_PAH(i))
         betadimer = 4.0 * sqrt(2.0) * size_PAH(i)**(1.0_WP/6.0_WP) * &
                    (6.0_WP/Pi)**(2.0_WP/3.0_WP) * AVOGADRO * &
                     sqrt(Pi*Rgas*temp_/(2.0_WP*AVOGADRO*SootDensity)) * &
                    (MolarMassSoot/(AVOGADRO*SootDensity))**(1.0_WP/6.0_WP)
         rhodot    = rhodot - W_GAS(i_PAH(i)) * gamma_PAH(i) * betadimer * aromconc * aromconc
         end if
      end do
   else
      ! Nucleation
      rhodot = rhodot - W_GAS(i_C2H2) * (A_NUCL*exp(-E_NUCL/temp_)) * (dens_*Y_(i_C2H2)/W_GAS(i_C2H2))
      rhodot = rhodot + W_GAS(i_H2) * (A_NUCL*exp(-E_NUCL/temp_)) * (dens_*Y_(i_C2H2)/W_GAS(i_C2H2))

      ! Surface Growth
      A_s = Pi*(6.0_WP*dens_*moment_(2)/Pi/SootDensity)**(2.0_WP/3.0_WP)*(dens_*moment_(1))**(1.0_WP/3.0_WP)
      rhodot = rhodot - W_GAS(i_C2H2) * (A_GROWTH*exp(-E_GROWTH/temp_)) * sqrt(A_s) * (dens_*Y_(i_C2H2)/W_GAS(i_C2H2))
      rhodot = rhodot + W_GAS(i_H2) * (A_GROWTH*exp(-E_GROWTH/temp_)) * sqrt(A_s) * (dens_*Y_(i_C2H2)/W_GAS(i_C2H2))
      ! Oxidation by O2
      rhodot = rhodot - 0.5_WP*W_GAS(i_O2) * (A_OXI*sqrt(temp_)*exp(-E_OXI/temp_)) * A_s * (dens_*Y_(i_O2)/W_GAS(i_O2))
      rhodot = rhodot + W_GAS(i_CO) * (A_OXI*sqrt(temp_)*exp(-E_OXI/temp_)) * A_s * (dens_*Y_(i_O2)/W_GAS(i_O2))
!      rhodot = rhodot - 0.5_WP*W_GAS(i_O2) * 120.0_WP * (ka*(dens_*Y_(i_O2)/W_GAS(i_O2))*x/(1+kz*(dens_*Y_(i_O2)   &
!               /W_GAS(i_O2))) + kb*(dens_*Y_(i_O2)/W_GAS(i_O2))*(1-x))
!      rhodot = rhodot + W_GAS(i_CO) * 120.0_WP * (ka*(dens_*Y_(i_O2)/W_GAS(i_O2))*x/(1+kz*(dens_*Y_(i_O2)/         &
!               W_GAS(i_O2))) + kb*(dens_*Y_(i_O2)/W_GAS(i_O2))*(1-x))
      ! Oxidation by OH
      rhodot = rhodot - W_GAS(i_OH) * phi_OH * K_OXI_OH / sqrt(temp_) * A_s * (dens_*Y_(i_OH)/W_GAS(i_OH))
      rhodot = rhodot + W_GAS(i_CO) * phi_OH * K_OXI_OH / sqrt(temp_) * A_s * (dens_*Y_(i_OH)/W_GAS(i_OH))
      rhodot = rhodot + W_GAS(i_H)  * phi_OH * K_OXI_OH / sqrt(temp_) * A_s * (dens_*Y_(i_OH)/W_GAS(i_OH))
      ! Oxidation by O
      rhodot = rhodot - W_GAS(i_O)  * phi_O * K_OXI_O / sqrt(temp_) * A_s * (dens_*Y_(i_O)/W_GAS(i_O))
      rhodot = rhodot + W_GAS(i_CO) * phi_O * K_OXI_O / sqrt(temp_) * A_s * (dens_*Y_(i_O)/W_GAS(i_O))

   end if

   return
end subroutine soot_chem_rhodot

subroutine soot_chem_destroy
   use soot_mod
   use soot_chem_mod
   implicit none

!   deallocate(SC)
!   deallocate(SC_name)
!   deallocate(W_GAS)
   deallocate(MassFracs)

   return
end subroutine soot_chem_destroy 
