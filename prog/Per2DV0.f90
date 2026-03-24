
  INCLUDE "Basi/my_vrbls.f90"
  !=====Main routine =================================================
  ! Features of this version:
  ! -- Newman-Ziff algorithm 
  ! -- neighbours are not stored, but calculated

  PROGRAM sqa_2bond
    use my_vrbls
    implicit none
    integer :: itoss,isamp,iblck,iense, i,ii
    integer(4):: bk, b1 
    integer(4):: k, j, j1
		character(80):: str0, fname1,fname2
		double precision:: t_tot
    double precision :: Pat, NorSamp, tmp1, tmp2, tmp3, tmp4

    print *, 'L,  Nsamp, Prob, NBlck, Seed'
    read  *,  L,  Nsamp, Prob, NBlck, Seed
    Totsamp = Nsamp*NBlck


		
    !--- Initialization ----------------------------------------------
    call set_time_elapse
    call initialize
    call time_elapse
    t_init = t_elap
    write(6,40) t_init
    40 format(/'        set up time:',f16.7,2x,'s')

    wv = 1.d0/Vol
    L = Floor((Vol*1.d0)**0.5); Lh = L/2 
    nkave(1:nnb*L)= 0.d0;  Scave(1:nnb*L) =0.d0;  dcave(1:nnb*L) = 0.d0;  c1ave(1:nnb*L)= 0.d0 
    nkave2(1:nnb*L)= 0.d0; Scave2(1:nnb*L) =0.d0; dcave2(1:nnb*L) = 0.d0; c1ave2(1:nnb*L)= 0.d0 
    call init_His
!		call def_latt
    at_cri = .true.
    !--- Simulation --------------------------------------------------
		do iblck = 1, NBlck 
      write(*,*) iblck
      nmHis = 0
    	DO isamp = 1, Nsamp
     		call markov;            
				call measure
				call coll_data(iblck)
    	ENDDO
      if(at_cri) call writ_his
			call time_elapse; t_simu = t_simu + t_elap 
			call norm_Nsamp(iblck)	
		enddo 

    if(at_cri) then 
      call writ_his2
      write(str0, '(I8)') Vol 
      fname1 = "ERNss_V"//trim(adjustl(str0)) 
      open(100, file=fname1, access="append") 
      do i =1 ,nnb*L 
        Pat = Lh*i*wv; NorSamp= 1.d0/nmHist; 
        nkave(i) = nkave(i)*NorSamp; Scave(i) = Scave(i)*NorSamp; dcave(i) = dcave(i)*NorSamp; c1ave(i) = c1ave(i)*NorSamp
        nkave2(i) = nkave2(i)*NorSamp; Scave2(i) = Scave2(i)*NorSamp; dcave2(i) = dcave2(i)*NorSamp; c1ave2(i) = c1ave2(i)*NorSamp
        tmp1 = nkave2(i)-nkave(i)**2.d0; tmp1 = ( tmp1/(nmHist-1) )**0.5d0  
        tmp2 = Scave2(i)-Scave(i)**2.d0; tmp2 = ( tmp2/(nmHist-1) )**0.5d0  
        tmp3 = dcave2(i)-dcave(i)**2.d0; tmp3 = ( tmp3/(nmHist-1) )**0.5d0  
        tmp4 = c1ave2(i)-c1ave(i)**2.d0; tmp4 = ( tmp4/(nmHist-1) )**0.5d0  
        write(100, 42) i, Pat, nkave(i), tmp1, Scave(i), tmp2, dcave(i), tmp3, c1ave(i), tmp4  
      enddo 
      write(100,*) 
      close(100)
    42 format(i8, 9e22.14)
    endif

    !--- Statistics --------------------------------------------------
    call stat_analy
    call write2file

    CONTAINS

  INCLUDE "Init/initialize.f90"
  INCLUDE "Simu/markov.f90"
	INCLUDE "Meas/measure.f90" 
	INCLUDE "Meas/statistics.f90"
  INCLUDE "RNG/my_rng.f90"    

END PROGRAM sqa_2bond
!=====================================================================
