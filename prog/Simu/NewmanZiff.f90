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
		ppi = 0
		LP_nnb:do  
			if( e == 0 ) exit LP_nnb    
			if(ptr(k)<0) then 
				cs1 = -ptr(k)
			else 
				cs1 = -ptr(ptr(k)) 
			endif
			j =  k + E2b(e) 
			if(ptr(j) /= 0) then  
				j1 = j  
				ptr_cycle:do 
					if(ptr(j1)<0) exit ptr_cycle 
					ppi = ppi +1 
					PathStorage(ppi) = j1 
					if(j1==ptr(j1)) stop
					j1 = ptr(j1) 
				enddo ptr_cycle 
				if(j1 /= ptr(k)) then 
					cs2 = -ptr(j1) 
					if(cs1 <=cs2) then 
						rootpoint = j1 
						nrp= k ; if(ptr(k)>0) nrp = ptr(k) 
					else 
						rootpoint = ptr(k) 
						nrp = j1
					endif 
					Csc2 = Csc2 + 2.d0*(cs1*1.d0)*(cs2*1.d0) 
					cs1 =cs1 + cs2 
					ptr(rootpoint) = -cs1
					ptr(nrp) = rootpoint
					ptr(k)   = rootpoint
					nk = nk -1
				endif 
			endif 
			e = E2E(e)
		enddo LP_nnb 

		! path compression 
		if(ppi > 0) then 
			path_comp:do pi =1, ppi 
				tmp = PathStorage(pi) 
				if(tmp==rootpoint) cycle path_comp
				ptr(tmp) = rootpoint 
			enddo path_comp 
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

!			C1 = ck 
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