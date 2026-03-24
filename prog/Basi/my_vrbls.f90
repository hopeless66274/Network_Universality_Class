
  !*******************************************************************
  ! Bond Percolation on the 5D Simple-Cubic Lattice

  ! Error bars are calculated using the blocking technique. 
  ! Let ' \chi^2 = \sum_{t=1}^{T} <(O_t -<O>)^2> ' be the squared chi for
  ! 'T' blocks of observable 'O'. Assuming each block of data is independent
  ! of each other, the error bar of 'O' is given by 'Dev = \sqrt{\chi^2/T(T-1)}.

  ! Reliabity of the obtained errors is monitored by t=1 correlation,
  ! for which tolerance is set by variable 'tol' (default: tol=0.20d0).

  ! Composite quantities like Binder ratios are calculated in each block, and
  ! the associated error bars are obtained from their fluctuations.

  ! Results are written into a special file 'dat.***' if the number of
  ! blocks is less than 125 or correlation is too big. Data in each 
  ! block will be also printed out in this case.

  ! Default number of extensive simulation is 'NBlck=1024'.

  ! For test purpose, for which huge amount of information will be 
  ! printed out, 'NBlck' should be set smaller but >2.

  ! Dynamical behavior is not studied.

  ! 'my_vrbls.f90', 'carlo.f90', 'monte.f90', 'measure.f90', 
  ! 'write2file.f90' and etc need to be modified for new projects.

  !  Look for 'PROJECT-DEPENDENT'.

  !*******************************************************************

  ! Look for 'PROJECT-DEPENDENT' for different projects
  !============== variable module ====================================
  MODULE my_vrbls
    IMPLICIT NONE
    !-- common parameters and variables ------------------------------
    ! PROJECT-INDEPENDENT parameters
    double precision, parameter :: pi=4.0d0*atan(1.0d0)
    ! THIS IS ALMOST PROJECT-INDEPENDENT 
    double precision, parameter :: tm32   = 1.d0/(2.d0**32)
    double precision, parameter :: eps    = 1.d-14         ! very small number
    double precision, parameter :: tol    = 1.20d0         ! tolerance for Cor
    logical                     :: prt                     ! flag for write2file
    integer,          parameter :: Mxint  = 2147483647     ! maximum integer
    integer,          parameter :: Mnint  =-2147483647     ! minimum integer

    integer,          parameter :: MxBlck = 2**20          ! maximum number of blocks for statistics
    integer,          parameter :: MnBlck = 2              ! minimum number of blocks

		integer						 :: Nense 													! ensembles 
    integer            :: NBlck                            ! # blocks
    integer            :: Nsamp                            ! # samples in unit 'NBlck'
    integer            :: Totsamp
    integer            :: Ntoss                            ! # samples to be thrown away

    double precision   :: rmuv,radv                        ! auxillary variables to draw a random vertex
    double precision   :: rmub,radb                        ! auxillary variables to draw a neighboring vertex

    character(8)  :: ident     = 'PerER V0'                ! identifier
    character(14) :: datafile  = 'PerERV0.dat'             ! datafile
    character(11) :: datafile1 = 'dat.big_cor'             ! datafile for big correlation
 !   character(14) :: datafile2 = 'Per7DV0N_s.dat'
    !-----------------------------------------------------------------

    !-- parameters and variables -------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, parameter :: D = 2                            ! dimensionality
    integer, parameter :: nnb  = 2*D                       ! # neighboring vertices

    integer            :: L,  Lm1, Lm2, Lh, L4, L4h, Vol, VoE, Volm1
    double precision   :: wv, we                           ! 1/Vol, 1/E
    double precision   :: pb, fpb                          ! bond prob. 
    !-----------------------------------------------------------------

    !-- Lattice and State --------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    !! 3^7=2187;   2^15=32768
    integer, parameter :: MxL  =4096                       ! maximum linear size
    integer(4), parameter :: MxV  = MxL**D                    ! maximum volumn
		integer(4), parameter :: MxE  = MxV*20


		integer(4), dimension(1:MxV)    :: ptr 
!		integer(4), dimension(1:MXV)		:: siteorder
		integer(4), dimension(1:MxE)		:: siteorder
!		integer(4), dimension(0:1000000):: PathStorage
 		integer(4),	dimension(MxV):: V2E, Vbn 
 		integer(4), dimension(MxE):: E2E, E2V, E2b 
    integer(4), dimension(1:MxV):: gaplist

    integer(4), dimension(1000):: HisMDC1,HisMNC1,HismaxDC1, Hismaxnk,HisPmaxDC1, HisPmaxnk

    !-----------------------------------------------------------------

    !-- Histogram ----------------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer(4), parameter :: MxHis  = 1920
    integer(4), parameter :: unHis  = 64
    integer(8)            :: szHis(0:MxHis)
    logical               :: at_cri
    integer(4)            :: nmHis, nmHist
    integer(4)            :: MiMDC1, MxMDC1, MiMNC1, MxMNC1  
    double precision      :: His_MDC1(MxHis), His_MNC1(MxHis)
    integer(4)            :: MimaxDC1, MxmaxDC1, Mimaxnk, Mxmaxnk
    double precision      :: His_maxDC1(MxHis), His_maxnk(MxHis)
    integer(4)            :: MiPmaxDC1, MxPmaxDC1, MiPmaxnk, MxPmaxnk
    double precision      :: His_PmaxDC1(MxHis), His_Pmaxnk(MxHis)

    integer(4):: MiHisDC1, MxHisDC1 
    double precision :: nsHisDC1(MxHis), ns1HisDC1(MxHis)

    !-- Observables --------------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    
    integer, parameter  :: NObs_b = 24                        ! #basic     observables
    integer, parameter  :: NObs_c = 9                         ! #composite observables
    integer, parameter  :: NObs = NObs_b+NObs_c               ! #total     observables
!    double precision, dimension(1:MxL) :: fsin, fcos
    double precision   :: S2,  S4, Cs1, Sh1, Ux1, Cs ,Csm, tot
    double precision   :: nk,nkL, ns0, ns1, ns2, nsh, nc0, nc1
    double precision   :: C1,C2,C3 
		double precision 	 :: Nor
		double precision  :: MaxDC1, Pmaxdc1, Maxnk, Pmaxnk, MaxS2, PmaxS2, MDC1, MNC1, MSC1
		double precision :: Pc

    double precision, dimension(1:nnb*MxL)::nkave, nkave2, Scave, dcave, c1ave, Scave2, dcave2, c1ave2 
		double precision   :: Gamma, Prob

    double precision :: nmGap, nmGap1, SGap, S1Gap

    !   1. Rx       2. R2,         3. R3       4. M^2       5. M^4      

    ! Quan.35-50 
    !  36. Qm      37. Qs         38. Q1      39. xim     
    !  40. Qe      41. Ce         42. Re   
    !-----------------------------------------------------------------

    !-- Statistics ---------------------------------------------------
    ! THIS IS PROJECT-DEPENDENT 
!		integer, parameter ::	Obs = 4
!		double precision, dimension(MxV,Obs) :: Quan, Quan2 
!		double precision, dimension(MxV,Obs) :: Ave, Dev
    double precision, dimension(NObs_b)      :: Quan            ! Measured quantities
    double precision, dimension(NObs,MxBlck) :: Obs             ! 1st--#quan.  2nd--#block
    double precision, dimension(NObs)        :: Ave, Dev, Cor   ! average, error bars, and correlation of observables
 !   real(4),dimension(1:Mxs)  :: Obs_ns
!    integer(4) :: smax=0
    !-----------------------------------------------------------------

    !-- Random-number generator---------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    integer, parameter           :: mult=32781
    integer, parameter           :: mod2=2796203, mul2=125
    integer, parameter           :: len1=9689,    ifd1=471
    integer, parameter           :: len2=127,     ifd2=30
    integer, dimension(1:len1)   :: inxt1
    integer, dimension(1:len2)   :: inxt2
    integer, dimension(1:len1)   :: ir1
    integer, dimension(1:len2)   :: ir2
    integer                      :: ipnt1, ipnf1
    integer                      :: ipnt2, ipnf2
    integer, parameter           :: mxrn = 10000
    integer, dimension(1:mxrn)   :: irn(mxrn)

    integer                      :: Seed                   ! random-number seed
    integer                      :: nrannr                 ! random-number counter
    !-----------------------------------------------------------------

    !-- time-checking variables --------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    character( 8)         :: date
    character(10)         :: time
    character(5 )         :: zone
    integer, dimension(8) :: tval
    double precision      :: t_prev, t_curr, t_elap
    integer               :: h_prev, h_curr
    double precision      :: t_init, t_simu, t_meas, t_toss, t_anas
    !-----------------------------------------------------------------
  END MODULE my_vrbls
  !===================================================================
