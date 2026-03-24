  !=======initialize histogram =======================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE init_His
    implicit none
    integer(8) :: i, s, bn, szbn
    
    nmHis = 0 ; nmHist = 0

    MiMDC1 = Vol; MxMDC1 =0; 
    MiMNC1 = Vol; MxMNC1 =0;  
    MimaxDC1=Vol; MxmaxDC1 =0  
    Mimaxnk =Vol; Mxmaxnk = 0
    MiPmaxDC1 =Vol; MxPmaxDC1=0; 
    MiPmaxnk  =Vol; MxPmaxnk =0
    His_PmaxDC1 = 0; His_Pmaxnk=0 
    His_maxDC1 = 0;  His_maxnk =0
    His_MDC1 = 0;    His_MNC1  =0

    HisPmaxDC1(:) = 0; HisPmaxnk(:)=0 
    HismaxDC1 (:) = 0; Hismaxnk (:)=0
    HisMDC1   (:) = 0; HisMNC1  (:)=0
 
    MiHisDC1 = VOl; MxHisDC1 = 0 
    nsHisDC1 = 0
    ns1HisDC1 =0

      return
  END SUBROUTINE init_His
  !===================================================================



  !==========define histogram ========================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE def_szHis
    implicit none
    integer(8) :: i, s, bn, szbn
    double precision :: a
    szHis= 0
    szbn = 1; s = 0
    szHis(0)=0
    do i = 1, MxHis
      s  = s+szbn;    szHis(i) = s
      if(Mod(i,unHis)==0) szbn = szbn*2
!     write(6,'(i5,2x,i12)') i, szHis(i)
    enddo

    a   = 2.d0    
!    Mxm = (L/8.d0)**(-1.d0/3.d0)
!    szm = Mxm/nmBin

 !  write(6,'(2f18.10)') Mxm, szm

 !  do s = 1, 2000
 !  write(6,'(i12,2x,i5)') s, func_bnHis(s)
 !  bn = func_bnHis(s)
 !  enddo

    return
  END SUBROUTINE def_szHis
  !===================================================================





  !============= function sz2bn ======================================
  !! THIS IS PROJECT-DEPENDENT 
  integer(4) FUNCTION func_bnHis(S)
    implicit none
    integer(4)       :: S
    integer(4)       :: i, bn, dbn, szbn
    func_bnHis = MxHis;    IF(S>=szHis(MxHis)) RETURN

    i  = unHis
    LP_I: DO
     IF(i>=MxHis)    EXIT LP_I
     IF(S<=szHis(i)) EXIT LP_I
     i = i+unHis
    ENDDO LP_I
    i    = i-unHis
    bn   = i;        i = i/unHis
    szbn = 2**i

    dbn  = S-szHis(bn)
    dbn  = 1+(dbn-1)/szbn
    func_bnHis = bn+dbn
!   write(6,'(i8,2x,3i8,2x,i8)') S, szbn, bn, dbn, func_bnHis
    return
  END FUNCTION func_bnHis
  !===================================================================





  !============== Write histogram to disk==============================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE writ_his
    implicit none
    character*20     :: str0, str1, str2, str3, str4, str5
    character*120    :: fnam, fnam1
    integer(4)       :: s, i, bn, szbn, j, MxHs,MiHs
    double precision :: pro, nor, nmHs
    double precision :: s0,  a0,  y0, x0, wn
    integer(4):: tmp,tmp1
    character*80:: fstr 
!   write(*,*) Vol
     write(str1, '(I10)') Vol
     write(str3, '(I8)') nmHis
     fstr =  "("//trim(adjustl(str3))//"I12)"
     fnam1  ='OMDC1'//'_V'//trim(adjustl(str1))
     open (21,file=fnam1, access='append')
     write(21, fstr) HisMDC1(1:nmHis) 

     fnam1  ='OMNC1'//'_V'//trim(adjustl(str1))
     open (22,file=fnam1, access= 'append')
     write(22, fstr) HisMNC1(1:nmHis) 

     fnam1  ='OmaxDC1'//'_V'//trim(adjustl(str1))
     open (23,file=fnam1, access= 'append')
     write(23, fstr) HismaxDC1(1:nmHis) 

     fnam1  ='Omaxnk'//'_V'//trim(adjustl(str1))
     open (24,file=fnam1, access= 'append')
     write(24, fstr) Hismaxnk(1:nmHis) 

     fnam1  ='OPmaxDC1'//'_V'//trim(adjustl(str1))
     open (25,file=fnam1, access= 'append')
     write(25, fstr) HisPmaxDC1(1:nmHis) 

     fnam1  ='OPmaxnk'//'_V'//trim(adjustl(str1))
     open (26,file=fnam1, access= 'append')
     write(26, fstr) HisPmaxnk(1:nmHis) 
     close(21); close(22); close(23); close(24); close(25); close(26)

    return
  END SUBROUTINE writ_his
  !===================================================================

  SUBROUTINE writ_his2
    implicit none
    character*20     :: str0, str1, str2, str3, str4, str5
    character*120    :: fnam, fnam1
    integer(4)       :: s, i, bn, szbn, j, MxHs,MiHs
    double precision :: pro, nor, nmHs
    double precision :: s0,  a0,  y0, x0, wn
    integer(4):: tmp,tmp1

     write(str1, '(I10)') Vol
  
     MxHs = MxHisDC1 ; MiHs=MiHisDC1 
     nmHs = nmHist; nor  = 1.d0/nmHs
     fnam  = 'ns.HisDc1'//'_V'//trim(adjustl(str1))
     open (14, file=fnam,access='append')
     write(14,'(i8,i12,i14, 2x, e22.12)') Vol, MxHs-MiHs+1, int(nmHs + 1.d-1), nmGap 
     lp1:do j= MiHs,    MxHs
        bn  = (j-1)/unHis;      szbn=2**bn
        if (nsHisDC1(j) < 0.5*wv) cycle lp1
        pro = nsHisDc1(j)*nor/(szbn*1.d0)
        write(14,42) j,szHis(j),nsHisDc1(j), pro
     enddo lp1    
     write(14,*) 
     close(14)
   42 format(2x,i8,2x,i12,2x,2e20.10)

     MxHs = MxHisDC1 ; MiHs=MiHisDC1 
     nmHs = nmHist; nor  = 1.d0/nmHs
     fnam  = 'ns1.HisDc1'//'_V'//trim(adjustl(str1))
     open (15, file=fnam,access='append')
     write(15,'(i8,i12,i14, 2x, e22.12)') Vol, MxHs-MiHs+1, int(nmHs + 1.d-1), nmGap 
     lp2:do j= MiHs,    MxHs
        bn  = (j-1)/unHis;      szbn=2**bn
        if (ns1HisDC1(j) < 0.5*wv) cycle lp2
        pro = ns1HisDc1(j)*nor/(szbn*1.d0)
        write(15,42) j,szHis(j),ns1HisDc1(j), pro
     enddo lp2    
     write(15,*) 
     close(15)

    return
  END SUBROUTINE writ_his2
