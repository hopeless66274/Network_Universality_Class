
  !==============Measurement =========================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE measure
    implicit none
    integer :: Wx
    integer(4)::bn, tmp

 !   C1 = C1 *wv; C2 = C2 *wv; C3 = C3 *wv;        nk = nk *wv
 !   Csm= Cs /ns0;       Cs = Cs *wv

    if (at_cri) then 
        nmHis = nmHis + 1
        nmHist = nmHist + 1
        tmp = int(MDC1+1.d-3)
        HisMDC1(nmHis) =  tmp      
      
        tmp = int(MNC1+1.d-3)
        HisMNC1(nmHis) = tmp
      
        tmp =int(maxDC1+1.d-3) 
        HismaxDC1(nmHis) = tmp

        tmp = int(maxnk+1.d-3)
        Hismaxnk(nmHis) = tmp
        
        tmp =int(PmaxDC1+1.d-3)
        HisPmaxDC1(nmHis) = tmp
        
        tmp = int(Pmaxnk+1.d-3)
        HisPmaxnk(nmHis) = tmp
        
    endif

    MDC1 = MDC1*wv; MNC1 = MNC1*wv;  MSC1 = MSC1*wv
		MaxDc1 = MaxDc1*wv; Maxnk = Maxnk*wv;      maxS2 = maxS2*wv
		PmaxDc1 = PmaxDc1*wv;  Pmaxnk = Pmaxnk*wv; PmaxS2 = PmaxS2*wv

    Quan( 1) = MaxDc1  
    Quan( 2) = PmaxDc1 
    Quan( 3) = Maxnk 
    Quan( 4) = Pmaxnk 
    Quan( 5) = MaxS2 
		Quan( 6) = PmaxS2 
    Quan( 7) = abs(PmaxDc1- Pmaxnk)
    Quan( 8) = abs(PmaxS2-Pmaxnk) 
		Quan( 9) = abs(PmaxDc1-PmaxS2)
		Quan(10) = abs(Pmaxnk-Pc)
    Quan(11) = MDC1
    Quan(12) = MNC1
    Quan(13) = MSC1

    Quan(14) = MaxDc1**2.d0
    Quan(15) = PmaxDc1**2.d0
    Quan(16) = Maxnk**2.d0
    Quan(17) = Pmaxnk**2.d0
    Quan(18) = MaxS2**2.d0
		Quan(19) = PmaxS2**2.d0  
    Quan(20) = MDC1**2.d0 
    Quan(21) = MNC1**2.d0

    Quan(22) = nmGap1 
    Quan(23) = nmGap1**2.d0 
    Quan(24) = SGap

    return
  END SUBROUTINE measure
  !======================================================================

  !==============Calculate reweighting ==================================
  !! O=Ave(b1)/Ave(bw)
  !! THIS IS PROJECT-DEPENDENT
  SUBROUTINE cal_rew(jb,b1,bw)
    implicit none
    integer, intent(in) :: jb, b1,bw
    integer             :: k
    double precision    :: tmp

    !-- Average ----------------------------------------------------
      tmp = Ave(bw);          if(dabs(tmp)>eps) tmp = Ave(b1)/tmp
      Ave(jb) = tmp

    !-- Obs(j,k) series --------------------------------------------
    do k = 1, NBlck
      tmp = Obs(bw,k);        if(dabs(tmp)>eps) tmp = Obs(b1,k)/tmp
      Obs(jb,k) = tmp
    enddo
  END SUBROUTINE cal_rew
  !======================================================================


  !==============Calculate Binder ratio 1================================
  !! Q=Ave(b2)/Ave(b1)**2
  !! THIS IS PROJECT-DEPENDENT
  SUBROUTINE cal_Q(jb,b1,b2)
    implicit none
    integer, intent(in) :: jb, b1,b2
    integer             :: k
    double precision    :: tmp

    !-- Average ----------------------------------------------------
    tmp = Ave(b1)**2;    if(dabs(tmp)>eps) tmp = tmp/Ave(b2)
!    tmp = Ave(b1)**2;    if(dabs(tmp)>eps) tmp = Ave(b2)/tmp
      Ave(jb) = tmp

    !-- Obs(j,k) series --------------------------------------------
    do k = 1, NBlck
      tmp = Obs(b1,k)**2;  if(dabs(tmp)>eps) tmp = tmp/Obs(b2,k)
!      tmp = Obs(b1,k)**2;  if(dabs(tmp)>eps) tmp = Obs(b2,k)/tmp
      Obs(jb,k) = tmp
    enddo
  END SUBROUTINE cal_Q
  !======================================================================



  !==============Calculate Binder ratio 1================================
  !! Q=Ave(b2)*Ave(bw)/Ave(b1)**2
  !! THIS IS PROJECT-DEPENDENT
  SUBROUTINE cal_Q_rew(jb,b1,b2,bw)
    implicit none
    integer, intent(in) :: jb, b1,b2,bw
    integer             :: k
    double precision    :: tmp

    !-- Average ----------------------------------------------------
      tmp = Ave(b1)**2;    if(dabs(tmp)>eps) tmp = Ave(b2)  *Ave(bw)/tmp
      Ave(jb) = tmp

    !-- Obs(j,k) series --------------------------------------------
    do k = 1, NBlck   
      tmp = Obs(b1,k)**2;  if(dabs(tmp)>eps) tmp = Obs(b2,k)*Obs(bw,k)/tmp
      Obs(jb,k) = tmp
    enddo
  END SUBROUTINE cal_Q_rew
  !======================================================================


  !==============Calculate specific-heat-like quantity ==================
  !! THIS IS PROJECT-INDEPENDENT
  !! C = (Ave(b2)-Ave(b1)^2)*Vol
  SUBROUTINE cal_sp_heat(jb,b1,b2)
    implicit none
    integer, intent(in) :: jb, b1, b2
    integer             :: k
    double precision    :: tmp

    !-- Average ----------------------------------------------------
      Ave(jb) = Vol*(Ave(b2)-Ave(b1)**2.d0)

    !-- Obs(j,k) series --------------------------------------------
    do k = 1, NBlck
      Obs(jb,k) = Vol*(Obs(b2,k)-Obs(b1,k)**2.d0)
    enddo
  END SUBROUTINE cal_sp_heat
  !======================================================================


  !==============Calculate specific-heat-like quantity ==================
  !! THIS IS PROJECT-INDEPENDENT
  !! C = Ave(b2)/Ave(bw)-(Ave(b1)/Ave(bw))^2
  SUBROUTINE cal_sp_heat_rew(jb,b1,b2,bw)
    implicit none
    integer, intent(in) :: jb, b1, b2,bw
    integer             :: k
    double precision    :: tmp

    !-- Average ----------------------------------------------------
      tmp = Ave(bw);  if(dabs(tmp)>eps)     tmp = 1.d0/tmp
      Ave(jb) = (Ave(b2)*tmp)-(Ave(b1)*tmp)**2.d0
      Ave(jb) = Ave(jb)*Vol

    !-- Obs(j,k) series --------------------------------------------
    do k = 1, NBlck
      tmp = Obs(bw,k);    if(dabs(tmp)>eps) tmp = 1.d0/tmp
      Obs(jb,k) = (Obs(b2,k)*tmp)-(Obs(b1,k)*tmp)**2.d0
      Obs(jb,k) = Obs(jb,k)*Vol
    enddo
  END SUBROUTINE cal_sp_heat_rew
  !======================================================================



  !==============Calculate Convariance ==================================
  !! THIS IS PROJECT-INDEPENDENT 
   !! THIS IS PROJECT-INDEPENDENT
  !! Cov = V(Ave(b3)-Ave(b2)*Ave(b1))
  SUBROUTINE cal_conv(jb,b3,b2,b1)
    implicit none
    integer, intent(in) :: jb, b1, b2, b3
    integer             :: k

    !-- Average ----------------------------------------------------
      Ave(jb) = Vol*(Ave(b3)-Ave(b2)*Ave(b1))

    !-- Obs(j,k) series --------------------------------------------
    do k = 1, NBlck
      Obs(jb,k) = Vol*(Obs(b3,k)-Obs(b2,k)*Obs(b1,k))
    enddo
  END SUBROUTINE cal_conv
  !======================================================================


  !==============Calculate correlation length==================
 !! THIS IS PROJECT-INDEPENDENT
  SUBROUTINE cal_corr_length(jb,b2,b1)
    implicit none
    integer :: k
    integer, intent(in) :: jb, b1, b2
    double precision :: s
    s = 1.d0/(2.d0*L*sin(pi/L))

    !-- Average ----------------------------------------------------
      Ave(jb) = dsqrt(Ave(b1)/Ave(b2)-1)*s
     
    !-- Obs(j,k) series --------------------------------------------
    do k = 1, NBlck
      Obs(jb,k) = dsqrt(Obs(b1,k)/Obs(b2,k)-1)*s
    enddo
      
  END SUBROUTINE cal_corr_length
  !======================================================================



! !==============Calculate energy-like corr. length =====================
!!! THIS IS PROJECT-INDEPENDENT
! SUBROUTINE cal_corr_xie(jb,b3,b2,b1)
!   implicit none
!   integer :: k
!   integer, intent(in) :: jb, b1, b2, b3
!   double precision :: s, s0
!   s = 1.d0/(2.d0*L*sin(pi/L))

!   !-- Average ----------------------------------------------------
!     s0 = (Ave(b2)-Ave(b1)**2.d0)/Ave(b3)
!     s0 = s0-1.d0
!     Ave(jb) = dsqrt(Ave(b1)/Ave(b2)-1)*s
!     Ave(jb) = s0
!    
!   !-- Obs(j,k) series --------------------------------------------
!   do k = 1, NBlck
!     s0 = (Obs(b2,k)-Obs(b1,k)**2.d0)/Obs(b3,k)
!     s0 = s0-1.d0
!     Obs(jb,k) = dsqrt(Obs(b1,k)/Obs(b2,k)-1)*s
!     Obs(jb,k) = s0
!   enddo
      
! END SUBROUTINE cal_corr_xie
! !======================================================================




  !==============Calculate composite observables =====================
  !! THIS IS PROJECT-INDEPENDENT 
  !! call in 'stat_alan'
  SUBROUTINE cal_Obs_comp
    implicit none
    integer    :: jb, b1, b2, bw, i

!    jb = NObs_b+ 1;      call cal_Q(jb,1,2)           ! Q1
!    jb = NObs_b+ 2;      call cal_Q(jb,3,5)           ! Qs
     
     jb = NObs_b + 1;     call cal_sp_heat(jb, 1, 14)
     jb = NObs_b + 2;     call cal_sp_heat(jb, 2, 15)
     jb = NObs_b + 3;     call cal_sp_heat(jb, 3, 16)
     jb = NObs_b + 4;     call cal_sp_heat(jb, 4, 17)
     jb = NObs_b + 5;     call cal_sp_heat(jb, 5, 18)
     jb = NObs_b + 6;     call cal_sp_heat(jb, 6, 19)
     jb = NObs_b + 7;     call cal_sp_heat(jb,11, 20)
     jb = NObs_b + 8;     call cal_sp_heat(jb,12, 21)
     jb = NObs_b + 9;     call cal_sp_heat(jb,22, 23)
    
    return
  END SUBROUTINE cal_Obs_comp
  !===================================================================
