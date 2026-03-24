    !--- PROJECT-DEPENDENT ---------------------------------------------
  !==============Initialization ======================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE initialize
    implicit none

    !-- order should be changed --------------------------------------
    call tst_and_prt
    call set_RNG
    call def_latt
!    call def_prob
!   call def_conf
    call def_szHis
    !-- measurement initialization -----------------------------------
    Obs = 0.d0;   Quan = 0.d0
    Ave = 0.d0;   Dev  = 0.d0;     Cor = 0.d0
    return
  END SUBROUTINE initialize
 
  !==============Test and Print ======================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE tst_and_prt
    implicit none

    !-- Test and Print -----------------------------------------------
!    if((NBlck>MxBlck).or.(NBlck<MnBlck)) then
!      write(6,*) 'MnBlck <= NBlck <=MxBlck?';             stop
!    endif

!   if((NBlck>200).and.(NBlck/=MxBlck)) then
!     write(6,*) '"NBlck>200" is supposed for extensive &
!     & simulation. "NBlck=MxBlk" is suggested!';         stop
!   endif

    if(L>MxL) then
      write(6,*) 'L<=MxL?';                               stop
    endif

    write(6,40) Vol    
    40 format(' Percolation on V=',i4,2x,' ER lattice')

    write(6,41) pb    
    41 format(' bond prob:',f12.8)

    write(6,42) Nsamp
    42 format(' Will simulate      ',i10,2x,'steps ')

    write(6,43) NBlck
    43 format(' #Blocks            ',i10)

    write(6,44) Ntoss*NBlck
    44 format(' Throw away         ',i10,2x,'steps')

    write(6,45) Seed
    45 format(' Random-number seed ',i10)

    return
  END SUBROUTINE tst_and_prt

  
  !==============define flipping probability =========================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE def_prob
    implicit none
    integer :: i, mb

    fpb = 1.d0/dlog(1.d0-pb)
    write(6,40) pb,fpb;    40 format(1x,'pb,fpb:',2x,2ES20.8)

    mb  = 1+Floor(dlog(tm32)*fpb)
    write(6,'(2x,"Maximum nxt b:",i20)') mb

    return
  END SUBROUTINE def_prob




  !==============definite Lattice ====================================
  !! THIS IS PROJECT-DEPENDENT 
!           3  2
!           | /
!           |/
!    14 ----o---- 1
!          /|
!         / |
!       13  12
  SUBROUTINE def_latt
    implicit none
		integer(4) :: i, tmp, si  
		double precision:: rd 

!		call GenERGraph
!		call GenBAGraph
		call ReadGraph


    return
  END SUBROUTINE def_latt


! This subroutitne aims to generate the ER lattice with probability Prob
Subroutine GenERGraph
	implicit none 
	integer(4):: bk, k, e, j
	double precision:: rd

! initilaization
    Vol   = L**D;  !VoE = Vol*20   
	Prob  = 2.0/(1.d0*Vol)
!	wv  = 1.d0/Vol;   we = wv/(nnb*0.5d0)
!    Volm1 = Vol-1;   
!    Lm1   = L-1;       Lm2 = L-2;       
	V2E(1:Vol)=0; Vbn(1:Vol)=0; 
	E2E(1:VoE)= 0;E2V(1:VoE)=0;E2b(1:VoE) = 0


	Prob = 2.0/Vol
	rd = rn(); if (rd < tm32) rd = 1.d0; bk =1 + Floor(dlog(rd)/dlog(1-Prob)) 
	e =0; k =1 
	Place_Bond:do  
		if(k>Vol) exit Place_Bond
		if(bk>(Vol-k)) then 
			bk = bk-(Vol-k) 
			k = k + 1
			cycle Place_Bond
		endif 
		e = e+ 1 
		E2E(e) = V2E(k); E2V(e) = k; V2E(k) = e;  E2b(e) = bk; VBn(k) = Vbn(k) + 1 
		e = e+ 1 
		j = k + bk
		E2E(e) =V2E(j);  E2V(e) = j; V2E(j) = e;  E2b(e) = -bk; VBn(j) = Vbn(j) + 1
		rd = rn(); if(rd<tm32) rd = 1.d0; bk = bk + 1 + Floor(dlog(rd)/dlog(1-Prob))
	enddo Place_Bond
	VoE = e/2
	return 
end subroutine 

! This subroutine aims to generate the BA networks 
! m0: initial graph size 
! m =2 : the number of edges added in each step
subroutine GenBAGraph 
	implicit none 
	integer(4):: e, k, j, b, m0, m, TolE, bo
	integer(4):: j0, j0b
	double precision:: rd

! initilaization
	V2E(1:Vol)=0; Vbn(1:Vol)=0; 
	E2E(1:VoE)= 0;E2V(1:VoE)=0;E2b(1:VoE) = 0

	m0 =3; m =2 
	e = 0 

	do k = 1, m0-1 
		do j =k +1, m0 
			e = e+ 1 
			E2E(e) = V2E(k); E2V(e) = k; V2E(k) = e;  E2b(e) = j-k; VBn(k) = Vbn(k) + 1
			e = e+ 1 
			E2E(e) =V2E(j);  E2V(e) = j; V2E(j) = e;  E2b(e) = k-j; Vbn(j) = Vbn(j) + 1 
		enddo 
	enddo 

	! construct the BA network 
	do k = m0+1, Vol  
		j0 =0; j0b = 0 
		TolE = e 
		do b =1,m 
			rd = rn(); bo = Floor((TolE-j0b)*rd) + 1 
			lp_neighbor:do j =1, k-1 
				if(bo .le. Vbn(j)) exit lp_neighbor 
				bo = bo - Vbn(j) 
			end do lp_neighbor 
			j0 = j; j0b = VBn(j0)
			e = e+ 1 
			E2E(e) = V2E(k); E2V(e) = k; V2E(k) = e;  E2b(e) = j-k; VBn(k) = Vbn(k) + 1
			e = e+ 1 
			E2E(e) =V2E(j);  E2V(e) = j; V2E(j) = e;  E2b(e) = k-j; Vbn(j) = Vbn(j) + 1 
		enddo 
	enddo 
	return
end subroutine 

! This subroutine aims to geneate the  scale-free networks
Subroutine GenSFGraph 
 	implicit none 
 	integer(4):: i, TotE, k, Eor, j, e, e1, MaxDeg, Deg, kmin, kmax, InvR1
  	double precision:: rd, InVR, r11, factor, gamma, numtry, cutexpo
  	integer:: maxtry, try, v1, v2
  double precision:: Cum(0:(int(Vol**0.5)+1)) 
  integer, allocatable:: stub(:) 
  integer(4):: pt 

  gamma = 3.50 

  kmin = 2
    
  cutexpo = 1.0/2 ; if(gamma > 3) cutexpo = 1.0/(gamma-1)
  kmax = Floor( (Vol*1.0)**cutexpo + 1.d-3) 

  Cum(kmin-1:kmax) = 0.0 
  do k = kmin, kmax 
      Cum(k) = Cum(k-1) + 1.0/(k*1.0)**gamma 
  enddo 
  Cum =Cum/Cum(kmax) 

  Degree(1:Vol) = 0; TotE=0; MaxDeg = 0 

    lpV: do i =1, Vol 
      rd = rn() 

      lpk:do k = kmin, kmax 
        if(rd < Cum(k) ) then 
           Degree(i) = k
          TotE = TotE + k 
          exit lpk 
        endif 
      enddo lpk 
    enddo lpV 

    if(mod(TotE,2)==1) then 
      Degree(1) = Degree(1)+1  
      ToTE = TotE + 1 
    endif 

    allocate(stub(TotE)) 

    pt = 0
    do i =1, Vol 
      do j = 1, Degree(i) 
        pt = pt + 1
        stub(pt) = i
      enddo 
    enddo 
! shuffle stubs 

  do i = TotE, 2, -1 
    rd = rn() 
    j = int(rd*i) +1  
    
    pt = stub(i); stub(i) = stub(j) ; stub(j) = pt 
  enddo 

 	V2E(1:Vol)=0; Vbn(1:Vol)=0; 
  VoE = TotE
 	E2E(1:VoE)= 0;E2V(1:VoE)=0;E2b(1:VoE) = 0

  maxtry = 1000 
  e  = 0 
  do i =1, TotE, 2 
    v1 = stub(i); v2 = stub(i+1) 

    try =0 
    
    check_edge: do 

      if(v1 == v2) then 
        try = try + 1 
        if(try >maxtry) exit check_edge 

        rd = rn() 
        j = int(rd*(TotE/2))*2 +1 
        v2 = stub(j) 

        cycle 
      endif  

      ! checek multi_edge 

      pt= V2E(v1) 

      do while (pt /= 0) 
        if(v1 + E2b(pt) == v2) then 
          try = try +1 
          if(try > maxtry) exit check_edge 

          rd = rn() 
          j = int(rd*(TotE/2))*2 + 1 

          v2 = stub(j) 
          cycle check_edge 
        endif 

        pt = E2E(pt) 
      enddo 
      exit check_edge 
    enddo check_edge 

   e = e+1; 
   E2E(e)=V2E(v1); E2V(e)=v1; E2b(e)=v2-v1; V2E(v1)=e; Vbn(v1)=Vbn(v1)+1

   e = e+1; 
   E2E(e)=V2E(v2); E2V(e)=v2; E2b(e)=v1-v2; V2E(v2)=e; Vbn(v2)=Vbn(v2)+1

enddo

VoE1 = e/2

write(*,*) "Edges =",VoE1

deallocate(stub)

return
end subroutine  GenSFGraph

!! This subroutitne aims to generate the RR lattice with probability Prob
!Subroutine GenRRGraph
!	implicit none 
!	integer(4):: bk, k, e, j
!	double precision:: rd
!! initilaization
!	V2E(1:Vol)=0; Vbn(1:Vol)=0; 
!	E2E(1:VoE)= 0;E2V(1:VoE)=0;E2b(1:VoE) = 0
!
!	MeanDegree = 3 
!	k = 1
!	Place_Bond:do  k =1, Vol 
!		b = 0 
!		nei:do  
!			if(b > MeanDegree) cycle Place_Bond 
!			rd = rn()
!			j = k + 1 + int(Vol*rd); if(j>Vol) j = j-Vol 
!			ek = V2E(k) 
!			
!		enddo nei  
!		e = e+ 1 
!		E2E(e) = V2E(k); E2V(e) = k; V2E(k) = e;  E2b(e) = bk; VBn(k) = Vbn(k) + 1 
!		e = e+ 1 
!		j = k + bk
!		E2E(e) =V2E(j);  E2V(e) = j; V2E(j) = e;  E2b(e) = -bk; VBn(j) = Vbn(j) + 1
!		rd = rn(); if(rd<tm32) rd = 1.d0; bk = bk + 1 + Floor(dlog(rd)/dlog(1-Prob))
!	enddo Place_Bond
!	return 
!end subroutine GenRRGraph 


! This subroutitne aims to generate the ER lattice with probability Prob
Subroutine ReadGraph
	implicit none 
	integer(4):: bk, k, e, j
	double precision:: rd
	integer(4):: mv1, mv2, me1, stat1
	
	open(10, file="../sc-pkustk11.mtx")
	read(10,*)
	read(10,*) mv1, mv2, me1 

!	Vol = mv1;  VoE = me1*2
!	Pc = 2.0d0*Vol/VoE
! initilaization
	V2E(1:MxV)=0; Vbn(1:MxV)=0; 
	E2E(1:MxV)= 0;E2V(1:MxV)=0;E2b(1:MxV) = 0
	e =0
	Vol = 0
	ReadData:do 
		read(10,*, iostat=stat1) k, j 
		if(stat1 <0) exit ReadData
		if(k*j==0) cycle ReadData
		if(k>Vol) Vol =k 
		if(j>Vol) Vol =j
		if(k==j) cycle ReadData
		e = e+ 1 
		bk = j-k
		E2E(e) = V2E(k); E2V(e) = k; V2E(k) = e;  E2b(e) = bk; VBn(k) = Vbn(k) + 1 
		e = e+ 1 
		j = k + bk
		E2E(e) =V2E(j);  E2V(e) = j; V2E(j) = e;  E2b(e) = -bk; VBn(j) = Vbn(j) + 1
	enddo ReadData 

	
	VoE = e/2;  Pc = Vol/VoE
	if(VoE > MxE .or. Vol > MxV) then 
		write(6,*) "volume or edge exceeed!" 
	endif
	write(6,*) Vol, e/2, e 
!	stop 
!	if(e/2==me1)  write(*,*) "successful!"
!	write(*,*) e, me1	
	close(10)
	return 
end subroutine  ReadGraph


include "Init/def_His.f90"
