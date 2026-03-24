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
	integer(4):: Vertex(2), Rtex(2)


	do i =1,VoE  
		siteorder(i) = i 
	enddo 

	
	ptr(1:Vol) = -1
	C1 =0.d0; nk = 0.d0; Csc = 0.d0; Csc2 =0.d0

  	DC1 = 0; MaxDC1 = 0.d0; Pmaxdc1 = 0.d0; Maxnk = 0.d0; PMaxnk = 0.d0; PmaxS2 = 0.d0; MaxS2 = 0.d0

	MDC1 = 0.d0; MNC1 = 0.d0; MSC1 =0.d0 
	nmGap1 = 0.d0; SGap = 0.d0
	gapcount = 0; G1count = 0


	ii =0 
	LP_lattice: do 
		ii = ii + 1
		if(ii > VoE .or. abs(C1-Vol) < 1.d-2 ) exit LP_lattice

		rd = rn() ; si = ii  + Floor((VoE-ii+1)*rd) 
		tmp=siteorder(ii); siteorder(ii)=siteorder(si); siteorder(si) =tmp  

		e = 2*siteorder(ii)-1
		k = E2V(e);  j = k + E2b(e) 

		kr = find_root(k);  jr = find_root(j) 
		if(kr == jr) cycle LP_lattice
		cs1 = -ptr(kr); cs2 = -ptr(jr)

		if(cs1==1 .and. cs2==1) then 
			nk = nk +1  
		elseif(cs1 >1 .and. cs2 > 1) then 
			nk = nk -1 
		endif 

		if (cs1 .ge. cs2) then 
			Vertex(1) = k ; Vertex(2) = j 
			Rtex(1)  = kr ; Rtex(2) = jr 
		else 
			Vertex(1) = j ; Vertex(2)= k 
			Rtex(1)   = jr; Rtex(2)  = kr 
		endif 

		rootpoint = Rtex(1);  ptr(rootpoint) = ptr(kr) + ptr(jr) 
!		path compression 
		j1 = Vertex(2) 
		ptr_cycle:do 
			if(ptr(j1) <0 ) exit ptr_cycle 
			j2 = ptr(j1); ptr(j1) = rootpoint 
			j1 = j2 
		end do ptr_cycle
		ptr(j1) = rootpoint


		ck = -ptr(rootpoint)*1.d0; 

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