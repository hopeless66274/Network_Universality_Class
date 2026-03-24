!==== This program deduced to simulate the percolation with Newman-Ziff algorithm 

Subroutine NewmanZiff 
	implicit none 
	integer(4):: i, si, ii, tmp, b, k, j, j1, e, ppi,pi
	integer(4):: rootpoint, nrp
	integer(4):: cs1, cs2, DC1
	double precision:: rd  
	double precision:: ck, Csc, Csc2 ! cluster size cumulant
	character(80):: str0, fname1
	integer(4):: bn, mil
	integer(4):: kr, jr, j2
	integer(4):: gapcount, G1count

	do i =1,Vol  
		siteorder(i) = i 
	enddo 
	do i =1, Vol 
		rd = rn(); 
		si = i + Floor((Vol-i + 1)*rd) 
		tmp  = siteorder(i) 
		siteorder(i) = siteorder(si) 
		siteorder(si) = tmp 
	enddo 

	ptr(1:Vol) = 0 
	C1 =0.d0; nk = 0.d0; Csc = 0.d0; Csc2 =0.d0

  	DC1 = 0; MaxDC1 = 0.d0; Pmaxdc1 = 0.d0; Maxnk = 0.d0; PMaxnk = 0.d0; PmaxS2 = 0.d0; MaxS2 = 0.d0

	MDC1 = 0.d0; MNC1 = 0.d0; MSC1 =0.d0 
	nmGap1 = 0.d0; SGap = 0.d0
	gapcount = 0; G1count = 0

	LP_lattice: do ii =1,Vol  
		k = siteorder(ii) 
		ptr(k)  = -1 
		nk = nk + 1
		Csc = Csc + 1.0
		Csc2 = Csc2 + 1.d0**2.d0
		if(Vbn(k) < 1) cycle LP_lattice
		e = V2E(k)

		kr = find_root(k); cs1 = -ptr(kr);  rootpoint = kr
		LP_nnb:do  
			if( e == 0 ) exit LP_nnb    
			j = k + E2b(e) 
			e = E2E(e) 
!			write(*,*) ii, k, j, e, ptr(j)
			
			if(ptr(j)==0) cycle LP_nnb 
			kr = find_root(k); cs1 = -ptr(kr); rootpoint= kr
			jr = find_root(j); 
			if(jr == kr) cycle LP_nnb 
			cs2 = -ptr(jr) 
			rootpoint = kr ;  nrp = jr 
			if(cs1 <= cs2) then 
				rootpoint = jr ; nrp = kr 
			endif 

			Csc2 = Csc2 + 2.d0*(cs1*1.d0)*(cs2*1.d0) 
			cs1 =cs1 + cs2 
			ptr(nrp) = rootpoint 
			ptr(rootpoint) = -cs1
			nk = nk -1
		end do LP_nnb 
	
		e = V2E(k) 
		path_comp:do 
			if( e==0) exit path_comp 
			j = k + E2b(e) 
			e = E2E(e)
			if(ptr(j) == 0) cycle path_comp
			j1 = j 
			ptr_cycle1:do 
				if(ptr(j1)<0 .or. ptr(j1)==j1 ) exit ptr_cycle1 
				j2 = ptr(j1);  ptr(j1) = rootpoint 
				j1 = j2  
			enddo ptr_cycle1
			ptr(j1) = rootpoint 
		end do path_comp

		ptr(rootpoint) = -cs1

		if(cs1 > ii) then 
		write(*,*)  cs1, ii 
		stop 
		endif
		ck = cs1*1.d0; 
		if(ck > C1) then
			DC1 = ck -C1 
			C1 = ck
			gapcount = gapcount + 1 
			gaplist(gapcount) = DC1
			if(at_cri) then 
				nmGap = nmGap + 1 
				nmGap1 = nmGap1 + 1 
				bn = func_bnHis(DC1)
				if(bn > MxHisDC1) MxHisDC1 = bn 
				if(bn < MiHisDC1) MiHisDC1 = bn 
				nsHisDC1(bn) = nsHisDC1(bn) + wv
				SGap = SGap + (DC1*wv)**2.d0
			endif

			if(DC1 > MaxDc1) then  
				MaxDc1  = DC1 
				PmaxDC1 = ii 
				MDC1 = C1
				G1count = gapcount
			endif 
		endif 
		if(nk > Maxnk) then 
			Maxnk = nk 
			PMaxnk = ii 
			MNC1 = C1
		endif 

		S2 = Csc2-C1**2.d0 
		if(S2 > MaxS2) then 
			MaxS2 = S2 
			PmaxS2 = ii 
		endif 

	if(mod(ii,Lh)==0) then 
		mil = ii/Lh 
		 nkave(mil) =  nkave(mil) + nk 
		nkave2(mil) = nkave2(mil) + nk**2.d0 
		 Scave(mil) =  Scave(mil) + (csc/nk) 
		Scave2(mil) = Scave2(mil) + (csc/nk)**2.d0 
		 dcave(mil) =  dcave(mil) + DC1 
		dcave2(mil) = dcave2(mil) + DC1**2.d0 
		 c1ave(mil) =  c1ave(mil) + C1
		c1ave2(mil) = c1ave2(mil) + C1**2.d0 
	endif 


	enddo LP_lattice

	if (at_cri) then 
		do i =1, G1count
			DC1 = gaplist(i) 
			bn = func_bnHis(DC1) 
			ns1HisDC1(bn) = ns1HisDC1(bn) + wv
			S1Gap = S1gap + (DC1*wv)**2.d0
		enddo 
	endif 
	
	return
end subroutine 


integer function find_root(k) 
	implicit none 
	integer(4)::k, kr   
	kr = k 
	ft:do 
		if(ptr(kr) < 0) exit ft  
		 kr = ptr(kr) 
	enddo ft 

	find_root = kr
	return 
end function  find_root