
  !--- PROJECT-INDEPENDENT -------------------------------------------
  !==============Collect data ========================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE coll_data(iblck)
    implicit none
    integer, intent(in) :: iblck
    integer             :: j
    do j = 1, NObs_b
      Obs(j,iblck) = Obs(j,iblck)+ Quan(j)
    enddo
  END SUBROUTINE coll_data 

  !==============Normalize by Nsamp ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE norm_Nsamp(iblck)
    implicit none
    integer, intent(in) :: iblck
    integer             :: j
    double precision    :: nor
    nor = 1.d0/(Nsamp*1.d0)
    do j = 1, NObs_b
      Obs(j,iblck) = Obs(j,iblck)*nor
    enddo
  END SUBROUTINE norm_Nsamp
  !===================================================================

  !==============Statistics ==========================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE stat_analy
    implicit none
    integer          :: j, k, k0
    double precision :: devn, devp, nor

    ! -- calculate average -------------------------------------------
    nor  = 1.d0/(NBlck*1.d0)
    do j = 1, NObs_b
      Ave(j) = nor*Sum(Obs(j,1:NBlck))
    enddo

    Coarsen: do
      ! -- calculate error and t=1 correlation for basics obs.--------
      prt = .true.
      DO j = 1, NObs_b
        devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
        do k = 1,  NBlck
          devn   = Obs(j,k)-Ave(j)
          Dev(j) = Dev(j)+devn*devn
          Cor(j) = Cor(j)+devn*devp
          devp   = devn
        enddo 
        Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
        if(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
        Dev(j)   = dsqrt(Dev(j)/(NBlck-1.d0))
        if(dabs(Cor(j))>tol) prt = .false.
      ENDDO 

      IF(prt)                         EXIT Coarsen 
      IF(NBlck<=2)    THEN
        prt = .false.;                EXIT Coarsen 
      ENDIF

      ! -- coarsen blocking ------------------------------------------
      NBlck = NBlck/2;      nor = nor*2.d0
      DO j = 1, NObs_b
        k0 = 1
        do k   = 1, NBlck
          Obs(j,k) = (Obs(j,k0)+Obs(j,k0+1))*0.5d0
          k0 = k0 +2
        enddo 
      ENDDO 
    enddo Coarsen 

    ! -- define auxillary variables and average of composite obs.-----
    call cal_Obs_comp

    ! -- calculate error and t=1 correlation for composite obs.-----
    do j = 1+NObs_b, NObs
      devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
      DO k = 1,  NBlck
        devn   = Obs(j,k)-Ave(j)
        Dev(j) = Dev(j)+devn*devn
        Cor(j) = Cor(j)+devn*devp
        devp   = devn
      ENDDO
      Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
      IF(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
      Dev(j)   = dsqrt(Dev(j)/(NBlck-1.d0))
    enddo
    return
  END SUBROUTINE stat_analy
  !===================================================================

  !============== Write to files =====================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE write2file 
    implicit none
    integer       :: j, k, Nwri
    double precision :: t_tot

    !-- open file ----------------------------------------------------
    if(prt) then
      open (8,file=datafile,  access='append') 
    else
      open (8,file=datafile1, access='append') 
    endif
    write(8, *);              write(6,*)

    !-- write to data file--------------------------------------------
    write(8,40) ident, Vol, VoE/2, Totsamp, Seed
    40 format(a8,3i12,i8)

    do j = 1, Nobs
      write(8,41) j, Ave(j), Dev(j), Cor(j)
      write(6,41) j, Ave(j), Dev(j), Cor(j)
      41 format(i5,2e22.12,f12.5)
    enddo
    close(8)

    !-- write to output file if #block is too small-------------------
    if(NBlck<=2) then
      write(6,*)
      Nwri = NObs_b;       if(Nwri>5) Nwri = 5
      do k = 1, NBlck
        write(6,42) k,(Obs(j,k),j=1,Nwri)
        42 format(i4,5f16.8) 
      end do
    endif
    write(6,43)   NBlck
    43 format(2x,'NBlck:',i10)

    t_toss = t_toss/60.d0 
    t_simu = t_simu/60.d0
    t_meas = t_meas/60.d0
    t_tot  = t_toss+t_simu+t_meas
    write(6,50) t_toss, (t_toss/60.d0)
    50 format( '  equilibrate time:',f14.7,2x,'m. or ',f12.7,'  h.')
    write(6,51) t_simu, (t_simu/60.d0)
    51 format( '  simulation  time:',f14.7,2x,'m. or ',f12.7,'  h.')
    write(6,52) t_meas, (t_meas/60.d0)
    52 format( '  measure     time:',f14.7,2x,'m. or ',f12.7,'  h.')
    write(6,53) t_tot,  (t_tot/60.d0)
    53 format( '  total CPU   time:',f14.7,2x,'m. or ',f12.7,'  h.')
    return
  END SUBROUTINE write2file
  !===================================================================
