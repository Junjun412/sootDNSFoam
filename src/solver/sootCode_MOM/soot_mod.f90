module soot_mod ! harding code for Adelaide sooting flame
  use string
  use precision
  implicit none

  ! Size of problem
  integer, parameter :: dim = 2
  ! DQMOM
  integer, parameter :: dim_dqmom = 3
  integer, parameter :: nDeltas  = 2

  ! Number of aromatic species
  integer, parameter :: nPAH = 4         ! HMOM: number of PAH
  ! PAH indices
  integer, dimension(nPAH) :: i_PAH

  ! PAH carbon number
  integer, dimension(nPAH) :: size_PAH

  ! PAH collision efficiencies
  real(WP), dimension(nPAH) :: gamma_PAH

  ! PAH production rates
  real(WP), dimension(nPAH) :: pahprodrate

  ! Dimer carbon number
  real(WP) :: C_dimer

  ! Parameters from input file
  integer :: nMoments     ! HMOM/SEMI: =2 for semi; =4 for HMOM
  integer :: nEquations     ! HMOM/SEMI: =2 for semi; =4 for HMOM
  integer :: nSootGas     ! HMOM/SEMI: soot-related gas-composition
  integer :: nComposition !ns+1
  integer :: use_mom 
  logical :: soot_rho
  logical :: soot_rad
  logical :: use_pah
  logical :: use_nucleation
  logical :: use_coagulation
  logical :: use_condensation
  logical :: use_surfacereaction
  logical :: use_surfaceoxidation
  logical :: use_surfaceoxidation_OH
  logical :: use_surfaceoxidation_O
  logical :: use_fragmentation

  ! Trigonometric parameters
  real(WP), parameter :: Pi    = 3.1415926535897932385_WP
  real(WP), parameter :: twoPi = 6.2831853071795864770_WP

  ! Soot relevant constants
  real(WP), parameter :: SootDensity = 1800.0_WP    ! kg/m^3
  real(WP), parameter :: MolarMassSoot = 12.0e-3_WP ! kg/mol
  real(WP), parameter :: Avogadro = 6.022e23_WP     ! 1/mol
  real(WP), parameter :: Rgas = 8.314_WP            ! J/mol/K
  real(WP), parameter :: Df = 1.8_WP                ! 1
  real(WP), parameter :: chisoot = 1.0_WP * 2.3e19_WP        ! 1/m^2
  ! Semi-empirical model parameters
  real(WP), parameter :: A_NUCL   = 10000.0_WP          ! 1/s
  real(WP), parameter :: E_NUCL   = 21100.0_WP        ! K
  real(WP), parameter :: A_GROWTH = 600.0_WP        ! m^(3/2)/(m_S)/s
  real(WP), parameter :: E_GROWTH = 12100.0_WP      ! K
  real(WP), parameter :: A_OXI    = 1.0e4_WP           ! m^3/(m_S)^2/s
  real(WP), parameter :: E_OXI    = 19680.0_WP         ! K
  real(WP), parameter :: phi_OH   = 0.2_WP
  real(WP), parameter :: K_OXI_OH = 1270.0_WP
  real(WP), parameter :: phi_O    = 0.2_WP
  real(WP), parameter :: K_OXI_O  = 665.5_WP
  real(WP), parameter :: C_COAG   = 9.0_WP
  ! Derived constants
  real(WP) :: CarbonToDiam                          ! HMOM/Semi
  real(WP) :: Cfm                                   ! HMOM
  real(WP) :: Ccn                                   ! HMOM
  real(WP) :: lambda                                ! HMOM
  real(WP) :: av,as                                 ! HMOM
  real(WP) :: Cred                                  ! DQMOM
  real(WP) :: expDf                                 ! DQMOM

  ! Approximation of small surface area change
  real(WP) :: FitC                                  ! HMOM
  real(WP) :: FitExp                                ! HMOM

  ! Clipping for HMOM
  real(WP), parameter :: SMALLWEIGHT = 1.0e-20_WP
  real(WP), parameter :: BIG_NP = 1.0e6_WP
  real(WP), parameter :: SMALLSEP    = 1.1_WP

  ! Clipping for PAH equation
  real(WP), parameter :: SMALLY = 1.0e-20_WP

  ! ----------  Local environment -----------
  ! Particle sizes
  real(WP) :: DIMER_NBRC,DIMER_NBRH
  real(WP) :: NUCL_NBRC,NUCL_SURF,NUCL_NBRH
  real(WP) :: min_nbrC, min_nbrH

  ! Relevant quantities to soot formation
  real(WP) :: temp,dens,wmol                             ! HMOM/Semi
  real(WP), dimension(:), pointer :: MassFracs       ! HMOM/Semi
  real(WP) :: DIMER_conc,sqrtT,T_MU,MUsqrtW_RHOsqrtT ! HMOM
  real(WP) :: VISC                                   ! HMOM
  real(WP) :: c_C2H2, c_O2, c_OH, c_O                ! Semi

  ! HMOM/Semi moments
  real(WP), dimension(:), pointer :: mom, old_mom

  ! HMOM/Semi source terms
  real(WP), dimension(:),   pointer :: src_mom

  ! ------------- HMOM -------------
  ! Quantities from the gas phase
  real(WP) :: wCoeff,oxCoeff,o2Coeff,oxCoeff_1,oxCoeff_2
  real(WP) :: prodRate,local_nbrC
  real(WP) :: y_pah_table,src_chem_pos,src_chem_neg,src_dimer
  real(WP) :: prodRate_out,betaN_out,betaC_out
  ! HMOM moment indices
  real(WP), dimension(:),   pointer :: frac
  real(WP), dimension(:),   pointer :: kmom
  real(WP), dimension(:,:), pointer :: moments


  ! DQMOM quantities
  real(WP), dimension(:), pointer :: weight
  real(WP), dimension(:), pointer :: nbrC
  real(WP), dimension(:), pointer :: surf
  real(WP), dimension(:), pointer :: nbrH
  real(WP), dimension(:), pointer :: dp
  real(WP), dimension(:), pointer :: np

  ! DQMOM src terms and matrix
  real(WP), dimension(:,:), pointer :: matrix
  real(WP), dimension(:),   pointer :: src_dqmom
  ! working arrays
  real(WP), dimension(:), pointer :: sol
  real(WP), dimension(:), pointer :: old_sol


  ! ------------- CQMOM -------------
  ! Number of nodes
  integer :: n_nodeV,n_nodeS
  ! Number of moments to determine V and S
  integer :: n_momentV,n_momentS
  ! Number of total weight or abscissa (N_weight = N_abscissa)
  integer :: n_wa
  real(WP), dimension(:), pointer :: wt_V
  real(WP), dimension(:,:), pointer :: wt_S
  real(WP), dimension(:), pointer :: volume
  real(WP), dimension(:,:), pointer :: area
  real(WP), dimension(:,:), pointer :: Ail
  real(WP), dimension(:,:), pointer :: dp_cqmom
  real(WP), dimension(:,:), pointer :: np_cqmom

  ! ------------- Soot_finitechem ---------------
  ! Surf/ox indices (7)
  integer :: i_OH,i_H2O,i_H,i_H2,i_C2H2,i_O2,i_CO,i_O
  integer :: i_CO2,i_CH4

  ! PAH production rates
!  real(WP),allocatable :: W_GAS(:)
  real(WP),dimension(300) :: W_GAS

  ! Surface chemistry rate parameters 
  real(WP), parameter :: A1f_s = 6.72e1_WP / 1e6_WP
  real(WP), parameter :: n1f_s = 3.33_WP
  real(WP), parameter :: E1f_s = 6.09_WP * 1e3_WP / Rgas
  real(WP), parameter :: A1b_s = 6.44e-1_WP / 1e6_WP
  real(WP), parameter :: n1b_s = 3.79_WP
  real(WP), parameter :: E1b_s = 27.96_WP * 1e3_WP / Rgas
  real(WP), parameter :: A2f_s = 1.00e8_WP / 1e6_WP
  real(WP), parameter :: n2f_s = 1.80_WP
  real(WP), parameter :: E2f_s = 68.42_WP * 1e3_WP / Rgas
  real(WP), parameter :: A2b_s = 8.68e4_WP / 1e6_WP
  real(WP), parameter :: n2b_s = 2.36_WP
  real(WP), parameter :: E2b_s = 25.46_WP * 1e3_WP / Rgas
  real(WP), parameter :: A3f_s = 1.13e16_WP
  real(WP), parameter :: n3f_s = -0.06_WP
  real(WP), parameter :: E3f_s = 476.05_WP * 1e3_WP / Rgas
  real(WP), parameter :: A3b_s = 4.17e13_WP / 1e6_WP
  real(WP), parameter :: n3b_s = 0.15_WP
  real(WP), parameter :: E3b_s = 0.00_WP * 1e3_WP / Rgas
  real(WP), parameter :: A4_s = 2.52e9_WP / 1e6_WP
  real(WP), parameter :: n4_s = 1.10_WP
  real(WP), parameter :: E4_s = 17.13_WP * 1e3_WP / Rgas
  real(WP), parameter :: A5_s = 2.20e12_WP / 1e6_WP
  real(WP), parameter :: n5_s = 0.00_WP
  real(WP), parameter :: E5_s = 31.38_WP * 1e3_WP / Rgas
  real(WP), parameter :: eff6 = 0.13_WP
  ! Parameters for mixture fraction source term
  real(WP) :: nu_h,nu_c,nu_o2
  real(WP) :: fuel_c,fuel_h,fuel_o
  real(WP) :: oxi_c,oxi_h,oxi_o
  real(WP), parameter :: W_C = 0.01201_WP           ! kg/mol
  real(WP), parameter :: W_H = 0.001008_WP          ! kg/mol
  real(WP), parameter :: W_O = 0.0160_WP            ! kg/mol

  real(WP) :: aromaticConc
end module soot_mod
