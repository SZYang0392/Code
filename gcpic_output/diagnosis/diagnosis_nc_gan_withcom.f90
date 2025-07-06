!-----------------------------------------------------------------
! for global simulation 
! ic0_save 1; 1:fEflux; 2--fv;  3--gyro phase ---
! 17,181,19,20 --test particles
! 21,22,23   -- Trace particle Orbit
! 26, 27     --- Re-Trace particles

subroutine diagnos_field_nc(box9,iflag)
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use global_case
  use MHDsimulation
  implicit none
  include 'mpif.h'
  type(diag_type)   :: box9
  real,    pointer,dimension(:,:,:) :: plot
  type(bscalar_type), dimension(:), pointer    :: temp,beta_i,beta_e
  type(bvector_type), dimension(:), pointer    :: bfld,efld,ui,ue,delJ,dumm,av1,delB,delE
  type(bvector_type), dimension(:), pointer    :: J_BE,v_S,cv_S,PPPe
  real  :: param(80)
  integer :: ic,i,iplot,ifield,m,ifchecki,ifchecke,n,ifcheck(30),iflag,mn

  character*6 fldname(110)
  character*5 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(4),varid(50),dimids(3)
  save ic
  data ic/0/

  interface 
      subroutine allocate_bscalar(c,m)
        use constants
		use grid_data
	    type(bscalar_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bscalar

      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector

      subroutine dconvert(a,b,c)
		use grid_data
	    type(bscalar_type), dimension(:), pointer    :: a
        real,    pointer,dimension(:,:,:)            :: b
	    type(diag_type)                              :: c
      end subroutine dconvert


      subroutine Adconvert(a,b,c)
		use grid_data
	    type(bscalar_type), dimension(:), pointer    :: a
        real,    pointer,dimension(:,:,:)            :: b
	    type(diag_type)                              :: c
      end subroutine Adconvert


    subroutine dconvert_pqw(a,b)
       use constants
	   use grid_data
	   type(bscalar_type), dimension(:), pointer    :: a
       real,pointer,dimension(:,:,:) :: b
    end subroutine dconvert_pqw

    subroutine cconvert(a,b)
       use constants
      real(p2),pointer,dimension(:,:,:) :: a,b
    end subroutine cconvert
    subroutine cross(a,b,c)
       use constants
       real(p2),pointer,dimension(:,:,:) :: a,b,c
    end subroutine cross
    subroutine getue(a)
       use constants
      real(p2),pointer,dimension(:,:,:) :: a
    end subroutine getue

    subroutine filter_vector(a,m)
       use constants
      real(p2),pointer,dimension(:,:,:) :: a
	  integer :: m
    end subroutine filter_vector

    subroutine becellecenter(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine becellecenter

    subroutine bbcellecenter(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bbcellecenter

    subroutine bconcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bconcarE

    subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bconcarB

    subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bcovcarE

    subroutine gather_scalar(a,b)
       use constants
	   use grid_data
	   type(bscalar_type), dimension(:), pointer    :: a
       real(p2),pointer,dimension(:,:,:) :: b
    end subroutine gather_scalar
    subroutine get_j(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine get_j

    subroutine get_A(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine get_A

    subroutine get_JBE(a)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
    end subroutine get_JBE

    subroutine get_PPe(a)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
    end subroutine get_PPe

    subroutine getJe(m,a)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
	   integer :: m
    end subroutine getJe

    subroutine get_vdevdb(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
	   type(particle_type), dimension(:), pointer    :: b
    end subroutine get_vdevdb
    
	subroutine field3(a)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
    end subroutine field3

    subroutine bconvert_rtheta(a)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
    end subroutine bconvert_rtheta

    subroutine bcovconE(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bcovconE

    subroutine crosscon(a,b,c)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a,b,c
    end subroutine crosscon

    subroutine getUT(m,a)
       use constants
	   use grid_data
	   type(particle_type), dimension(:), pointer    :: a
    end subroutine getUT

	

  end interface

  fldname = (/"Bx    ","By    ","Bz    ","Ni    ","Ne    ","phi   ","A1x   ","A1y   ","A1z   ","Uepar ", & !10
              "Avt   ","Ex    ","Ey    ","Ez    ","Upx   ","Upy   ","Upz   ","Uex   ","Uey   ","Uez   ", & !20
              "Jx    ","Jy    ","Jz    ","B1x   ","B1y   ","B1z   ","E1x   ","E1y   ","E1z   ","EcB1  ", & !30
			  "EcB2  ","EcB3  ","EcB4  ","JcB8  ","JcB9  ","divP1 ","divP2 ","divP3 ","duet1 ","duet2 ", & !40
			  "ui1x  ","ui1y  ","ui1z  ","Ni(1) ","Ti1par","Ti1per","ui2x  ","ui2y  ","ui2z  ","Ni(2) ", & !50
			  "Ti2par","Ti2per","ui3x  ","ui3y  ","ui3z  ","Ni(3) ","T31par","Ti3per","ue1x  ","ue1y  ", & !60
			  "ue1z  ","Ne(1) ","Te1par","Te1per","ue2x  ","ue2y  ","ue2z  ","Ne(2) ","Te2par","Te2per", & !70
			  "betai ","betae ","divB  ","divE  ","ue2x  ","ue2y  ","ue2z  ","Ne(2) ","Te2par","Te2per", & !80
			  "edv11 ","edv12 ","edv21 ","edv22 ","edv31 ","edv32 ","edv11i","edv12i","edv21i","edv22i", & !90
			  "edv31i","edv32i","edv21 ","edv22 ","edv31 ","edv32 ","JB1   ","JE1   ","JB2   ","JE2   ", & !100
			  "PPe11 ","PPe22 ","PPe33 ","PPe12 ","PPe13 ","PPe23 ","JB11  ","JE11  ","JB21  ","JE21  "/)  !110

  ic     = ic +1
  if(iflag ==0 ) then
     idiagf  = idiagf + 1
     write(ct,'(I5.5)')IDIAGf 
  endif

  if(iflag ==1 ) then
     idiagf1 = idiagf1 + 1
     write(ct,'(I5.5)')IDIAGf1
  endif	  

  if(iflag ==2 ) then
     idiagf2 = idiagf2 + 1
     write(ct,'(I5.5)')IDIAGf2
  endif	  

  if(iflag ==3 ) then
     idiagf3 = idiagf3 + 1
     write(ct,'(I5.5)')IDIAGf3
  endif	  


  if( ghybrid) write(ct,'(I5.5)')int(stime + 2.*dt)

  param     = 0.
  param(1)  = stime  ! has been scaled to omega_i
  param(2)  = aelectron
  param(3)  = aion
  param(4)  = br0
  param(5)  = bt0
  param(6)  = bz0
  param(7)  = betae
  param(8)  = deltax !box9%dxd
  param(9)  = deltay !box9%dyd
  param(10) = deltaz !box9%dzd
  param(11) = kinds
! real unit ----
  param(12) = ubfield   !ulength/rhos  !because xdiag is in unit of rhos
  param(13) = uefield   !utime  /alfventime
  param(14) = udensity  !uspeed
  param(15) = uspeed    !uspeed
  param(16) = utime     !uspeed
  param(17) = utemperature
  param(18) = uj0
  param(19) = uphi
  param(20) = ulength*rhos
  param(21) = cva
  param(22) = lamda_i/rhos
  param(23) = lamda_e/rhos
  param(24) = uRe/rhos
  param(25) = lamda_de/rhos
  param(29) = va
  param(30) = alpha
! param(28) = ua0
! param(29) = ua0
! param(30) = uj0
! param(31) = ua0
  if(mod(simtype(2),100) >= 90) param(21:28) = bnvt_sw

! some key parameters
  param(41) = omega_lh/omega_i
  param(75-box9%nfield+1:75)= box9%DiagQuantity(1:box9%nfield) 
  param(77) = simtype(2)
  param(78) = simtype(1)
  param(79) = coordinate
  param(80) = box9%nfield   ! number of field quantity are recorded

  if(mype==0) then
!     if(mod(idiagbox,10) ==1) then

        if(iflag == 0) filename = 'field'//ct//'.nc'
        if(iflag == 1) filename = 'fieldfin'//ct//'.nc'
        if(iflag == 2) filename = 'fieldbox'//ct//'.nc'
        if(iflag == 3) filename = 'fieldinn'//ct//'.nc'
!       Create the file. 
        call check( nf90_create(FILENAME, nf90_clobber, ncid) )
        !创建netcdf文件FILENAME，nf90_clobber表示覆写，文件编号存储为ncid
! define parameters
        call check( nf90_def_dim(ncid, 'param', 80, dimid(1)) )
        !在文件ncid中定义维度，维度名为'param'，维度长度80，维度编号返回为dimid(1)
        !nf90_def_dim的返回值被送入check函数检查是否出错
        call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(1), varid(1)) )
        !在文件ncid中定义变量，变量名为param，数据类型NF90_REAL，变量维度为dimid(1)，变量名返回为varid(1)
        call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )
        !向文件ncid中的变量varid(1)赋予属性UNITS，属性值为'non'
        
        call check( nf90_def_dim(ncid, 'X', box9%nxd, dimid(2)) )
        call check( nf90_def_var(ncid, 'X', NF90_REAL, dimid(2), varid(2)) )
        call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

        call check( nf90_def_dim(ncid, 'Y', box9%nyd, dimid(3)) )
        call check( nf90_def_var(ncid, 'Y', NF90_REAL, dimid(3), varid(3)) )
        call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

        call check( nf90_def_dim(ncid, 'Z', box9%nzd, dimid(4)) )
        call check( nf90_def_var(ncid, 'Z', NF90_REAL, dimid(4), varid(4)) )
        call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

        dimids=(/dimid(2),dimid(3),dimid(4)/)
        do iplot = 1, box9%nfield
           ifield = box9%DiagQuantity(iplot)
           call check( nf90_def_var(ncid, fldname(ifield), NF90_REAL, dimids, varid(5+iplot)) )
           call check( nf90_put_att(ncid, varid(5+iplot), UNITS, 'non') )
        enddo

        call check( nf90_enddef(ncid) )

        call check( nf90_put_var(ncid, varid(1), param) )
        call check( nf90_put_var(ncid, varid(2), box9%xdiag) )
        call check( nf90_put_var(ncid, varid(3), box9%ydiag) )
        call check( nf90_put_var(ncid, varid(4), box9%zdiag) )
!	 else if(mod(idiagbox,10) ==0) then
!        filename = 'field'//ct//'.nc'
!       Create the file. 
!        call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
!        call check( nf90_def_dim(ncid, 'param', 80, dimid(1)) )
!        call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(1), varid(1)) )
!        call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

!        call check( nf90_def_dim(ncid, 'X', cuv%nxd, dimid(2)) )
!        call check( nf90_def_var(ncid, 'X', NF90_REAL, dimid(2), varid(2)) )
!        call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

!        call check( nf90_def_dim(ncid, 'Y', cuv%nyd, dimid(3)) )
!        call check( nf90_def_var(ncid, 'Y', NF90_REAL, dimid(3), varid(3)) )
!        call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

!        call check( nf90_def_dim(ncid, 'Z', cuv%nzd, dimid(4)) )
!        call check( nf90_def_var(ncid, 'Z', NF90_REAL, dimid(4), varid(4)) )
!        call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )
!        dimids=(/dimid(2),dimid(3),dimid(4)/)
!        do iplot = 1, nfield
!           ifield = DiagQuantity(iplot)
!           call check( nf90_def_var(ncid, fldname(ifield), NF90_REAL, dimids, varid(5+iplot)) )
!           call check( nf90_put_att(ncid, varid(5+iplot), UNITS, 'non') )
!        enddo

!        call check( nf90_enddef(ncid) )

!        call check( nf90_put_var(ncid, varid(1), param) )
!        call check( nf90_put_var(ncid, varid(2), cuv%xdiag) )
!        call check( nf90_put_var(ncid, varid(3), cuv%ydiag) )
!        call check( nf90_put_var(ncid, varid(4), cuv%zdiag) )

!	 endif
  endif
  allocate(plot(box9%nxd,box9%nyd,box9%nzd))
! if(mod(idiagbox,10) ==0) allocate(plot(cuv%nxd,cuv%nyd,cuv%nzd))

  call allocate_bscalar(temp,1) 
  call allocate_bscalar(beta_i,1) 
  call allocate_bscalar(beta_e,1) 
  call allocate_bvector(bfld,3) 
  call allocate_bvector(efld,3) 
  call allocate_bvector(ui,  3) 
  call allocate_bvector(ue,  3) 
  call allocate_bvector(delJ,3) 
  call allocate_bvector(dumm,3) 
  call allocate_bvector(av1, 3) 
  call allocate_bvector(delB,3) 
  call allocate_bvector(delE,3) 


! (1) b-field -------
  if(mod(box9%idiagbox,10) ==1 .or. mod(box9%idiagbox,10) ==3) then
     do m=1,mblocks
        dumm(m)%vector = bv1(m)%vector  + bv0(m)%vector 
     enddo
	 call bconcarB(dumm,bfld)
  else
     do m=1,mblocks
	    do n=1,3
		   if(df) then
		      dumm(m)%vector(n,:,:,:) =  bv1(m)%vector(n,:,:,:)  * &
		                                 block(m)%node%h(4,n,:,:,:) 
           else
              dumm(m)%vector(n,:,:,:) = (bv1(m)%vector(n,:,:,:)  + bv0(m)%vector(n,:,:,:))* &
		                                 block(m)%node%h(4,n,:,:,:) 
           endif
        enddo
     enddo
	 call bbcellecenter(dumm,bfld)
  endif
! (2) e-field -------
  if(mod(box9%idiagbox,10) ==1 .or. mod(box9%idiagbox,10) ==3) then
     do m=1,mblocks
        dumm(m)%vector = ev1(m)%vector 
     enddo
	 call bcovcarE(dumm,efld)
  else
     do m=1,mblocks
	    do n=1,3
           dumm(m)%vector(n,:,:,:) =  ev1(m)%vector(n,:,:,:)/block(m)%node%h(3,n,:,:,:)   
        enddo
     enddo
     call becellecenter(dumm,efld)
  endif
! (3) ui-field -------
  if(mod(box9%idiagbox,10) ==1 .or. mod(box9%idiagbox,10) ==3) then
     do m=1,mblocks
        dumm(m)%vector = ji(m)%vector 
     enddo
	 call bconcarE(dumm,ui)
  else
     do m=1,mblocks
	    do n=1,3
           dumm(m)%vector(n,:,:,:) = ji(m)%vector(n,:,:,:)*block(m)%node%h(3,n,:,:,:)
        enddo
     enddo
	 call becellecenter(dumm,ui)
  endif

  do m=1,mblocks
	 do n=1,3
        if(.not. df .and. simtype(1) /= 6 .and. simtype(1) /= 7 ) &
		    ui(m)%vector(n,:,:,:) = ui(m)%vector(n,:,:,:) /max(ni(m)%scalar,1.e-4_p2)
        if(df)ui(m)%vector(n,:,:,:) = ui(m)%vector(n,:,:,:) /max(ni0(m)%scalar,1.e-4_p2)
	 enddo
  enddo
! (4) ue-field --------------
  if(simtype(1) == 6 .or. simtype(1) == 7 .or. simtype(1) == 66 .or. simtype(1) == 67) then
!     do m=1,mblocks
!	    ue(m)%vector(1,:,:,:) = block(m)%ionf(1,:,:,:)%v(1)
!	    ue(m)%vector(2,:,:,:) = block(m)%ionf(1,:,:,:)%v(2)
!	    ue(m)%vector(3,:,:,:) = block(m)%ionf(1,:,:,:)%v(3)
!	 enddo
  else
     if(mod(box9%idiagbox,10) ==1 .or. mod(box9%idiagbox,10) ==3) then
        do m=1,mblocks
           dumm(m)%vector = je(m)%vector/qelectron + jef(m)%vector/qelectron  
        enddo
	    call bconcarE(dumm,ue)
     else
        do m=1,mblocks
	       do n=1,3
              dumm(m)%vector(n,:,:,:) = (je(m)%vector(n,:,:,:)+jef(m)%vector(n,:,:,:))/ &
			                            qelectron *block(m)%node%h(3,n,:,:,:)
           enddo
        enddo
	    call becellecenter(dumm,ue)
     endif

! only to get je --- electro current ---------------------------
     do m=1,mblocks
	    do n=1,3
           if(.not. df) then
		      ue(m)%vector(n,:,:,:) = ue(m)%vector(n,:,:,:) /max(ne(m)%scalar+nef(m)%scalar,1.e-4_p2)
		   endif
           if(df)ue(m)%vector(n,:,:,:) = ue(m)%vector(n,:,:,:) /max(ne0(m)%scalar,1.e-4_p2)
	    enddo
     enddo
  endif

  ifchecki = 0
  ifchecke = 0
  ifcheck  = 0
  do iplot = 1,box9%nfield
     ifield = box9%DiagQuantity(iplot)
	 if((ifield >= 40 .and. ifield <= 58) .or. (ifield >= 87 .and. ifield<=92)) then
	    ifchecki =   ifchecki + 1
	    if(ifchecki ==1) call getUT(0,ions)
     endif
	 if((ifield >= 59 .and. ifield<=70) .or.(ifield >= 81 .and. ifield<=86)  ) then
	    ifchecke =   ifchecke + 1
	    if(ifchecke ==1  ) call getUT(0,eles)
     endif

	 if(ifield  >= 21 .and. ifield  <= 23) then  ! delta_j
	    ifcheck(3) = ifcheck(3) +1
		if(ifcheck(3) ==1) then
		   if(simtype(1) >=6 .and. simtype(1) <=8 ) then   ! hybrid simulation j= curlB
  	          do m=1,mblocks
		         dumm(m)%vector    = bv1(m)%vector  + bv0(m)%vector 
				 if(Ldevice) dumm(m)%vector = bv1(m)%vector
		      enddo
		      call get_j(dumm,delJ)
		      do m=1,mblocks
	             do n =1,3
                    dumm(m)%vector(n,:,:,:)  = delJ(m)%vector(n,:,:,:) * block(m)%node%h(3,n,:,:,:)
		         enddo
			  enddo
	          call becellecenter(dumm,delJ)
              if(Ldevice) call bconvert_rtheta(delJ)

		   else     ! otherwise, J=Ji+Je
              if(mod(box9%idiagbox,10) ==1) then
                 do m=1,mblocks
                    dumm(m)%vector = ji(m)%vector + je(m)%vector 
                 enddo
	             call bconcarE(dumm,delJ)
              else
                 do m=1,mblocks
	             do n=1,3
                    dumm(m)%vector(n,:,:,:) = (ji(m)%vector(n,:,:,:)+je(m)%vector(n,:,:,:)) *block(m)%node%h(3,n,:,:,:)
                 enddo
                 enddo
	             call becellecenter(dumm,delJ)
              endif
              if(Ldevice) call bconvert_rtheta(delJ)

		   endif
        endif
	 endif

	 if(ifield  == 7 .or. ifield  == 8 .or. ifield  == 9 ) then  
	    ifcheck(4) = ifcheck(4) +1
		if(ifcheck(4) ==1) then
  	       do m=1,mblocks
		      dumm(m)%vector    = bv1(m)%vector  + bv0(m)%vector 
		   enddo
		   call get_A(dumm,Av1)
        endif
	 endif


	 if(ifield  == 24 .or. ifield  == 25 .or. ifield  == 26 ) then  ! delta_B
	    ifcheck(5) = ifcheck(5) +1
		if(ifcheck(5) ==1) then
           if(abs(coordinate)  ==1 .or.abs(coordinate)  ==2 .or. abs(coordinate)  ==3 .or. &
		      abs(coordinate)  ==8 .or.abs(coordinate)  ==5 .or. Ldevice) then
		      do m=1,mblocks
	             do n =1,3
                    dumm(m)%vector(n,:,:,:) = bv1(m)%vector(n,:,:,:) * block(m)%node%h(4,n,:,:,:)
		         enddo
			  enddo
	          call bbcellecenter(dumm,delB)

		   else
  	          call field3(delB)
           endif
           if(Ldevice) call bconvert_rtheta(delB)

        endif
	 endif

	 if(ifield  == 27 .or. ifield  == 28 .or. ifield  == 29 ) then  ! delta_E
	    ifcheck(6) = ifcheck(6) +1
		if(ifcheck(6) ==1) then
           if(abs(coordinate)  ==1 .or.abs(coordinate)  ==2 .or. abs(coordinate)  ==3 .or. &
		      abs(coordinate)  ==8 .or.abs(coordinate)  ==5 .or. Ldevice) then
		      do m=1,mblocks
	             do n =1,3
                    dumm(m)%vector(n,:,:,:)  = ev1(m)%vector(n,:,:,:) / block(m)%node%h(3,n,:,:,:)
		         enddo
			  enddo
	          call becellecenter(dumm,delE)
		   else

           endif
           if(Ldevice) call bconvert_rtheta(delE)

        endif
	 endif

	 if(ifield  == 71) then  ! beta_i
        do m=1,mblocks
		   beta_i(m)%scalar = 0.
		   do n=1,kinds
		      beta_i(m)%scalar = beta_i(m)%scalar + block(m)%ini(n,:,:,:)%ni * &
			                     (block(m)%ini(n,:,:,:)%Ti(1)+block(m)%ini(n,:,:,:)%Ti(2)*2)/3.
		   enddo
		   temp(m)%scalar = bfld(m)%vector(1,:,:,:)**2+bfld(m)%vector(2,:,:,:)**2+bfld(m)%vector(3,:,:,:)**2
        enddo
        do m=1,mblocks
		   beta_i(m)%scalar = beta_i(m)%scalar/max(temp(m)%scalar,1.e-4_p2)*alpha
        enddo
	 endif



	 if((ifield  >= 97 .and. ifield  <= 100) .or. (ifield  >= 81 .and. ifield  <= 84) ) then  ! J_BE will be saved
	    ifcheck(8) = ifcheck(8) +1
		if(ifcheck(8) ==1 .and. dipole ) then
           call allocate_bvector(J_BE,8) 
!		   call get_JBE(J_BE)
		   call get_vdEvdB(j_be,eles)
        endif
	 endif

	 if(ifield  >= 30 .and. ifield  <= 32) then  ! beta_i
	    ifcheck(9) = ifcheck(9) +1
		if(ifcheck(9) ==1 ) then
           call allocate_bvector(v_s,3) 
           call allocate_bvector(cv_s,3) 
 	       call bcovconE(ev1,dumm)
  	       call crosscon(dumm,bv1,v_S) 
		   call bcovcarE(v_s,cv_s)
        endif
     endif

	 if((ifield  >= 101 .and. ifield  <= 106)  ) then  ! J_BE will be saved
	    ifcheck(10) = ifcheck(10) +1
		if(ifcheck(10) ==1 .and. dipole ) then
           call allocate_bvector(PPPe,6) 
		   call get_PPe(PPPe)
        endif
	 endif


  enddo


  do iplot = 1, box9%nfield
     ifield = box9%DiagQuantity(iplot)
     select case(ifield)
        case(1,2,3)                           !bx,by,bz
		   do m=1,mblocks
		      temp(m)%scalar = bfld(m)%vector(ifield,:,:,:)
		   enddo
        case(4)                               !ni total ni
		   do m=1,mblocks
		      temp(m)%scalar = ni(m)%scalar
			  if(mhd) temp(m)%scalar = ns(m)%scalar
		   enddo
           do mn=1,kinds_fi
	          do m=1,mblocks
                 temp(m)%scalar =   temp(m)%scalar  + ionf(mn)%block(m)%n
	          enddo
           enddo
        case(5)                               !ne
		   do m=1,mblocks
		      temp(m)%scalar = ne(m)%scalar 
		      if(hybrid) temp(m)%scalar = nif(m)%scalar
		   enddo
           do mn=1,kinds_fe
	          do m=1,mblocks
                 temp(m)%scalar =   temp(m)%scalar  + elef(mn)%block(m)%n
	          enddo
           enddo
        case(6)                               !phi
		   do m=1,mblocks
		      temp(m)%scalar = phi(m)%scalar
		   enddo
        case(7,8,9)                           !a1
		   do m=1,mblocks
		      temp(m)%scalar = av1(m)%vector(ifield-6,:,:,:)
		   enddo
        case(10)                              !Ue||
!          temp = uepar     !Jpe(3,:,:)
        case(11)                              !az
!          call cconvert(avt,vec1)
!          temp = vec1(3,:,:,:)
        case(12,13,14)                        !E
		   do m=1,mblocks
		      temp(m)%scalar = efld(m)%vector(ifield-11,:,:,:)
		   enddo
        case(15,16,17)                        !ui
		   do m=1,mblocks
		      temp(m)%scalar = ui(m)%vector(ifield-14,:,:,:)
			  if(mhd) temp(m)%scalar = v(m)%vector(ifield-14,:,:,:)
		   enddo
        case(18,19,20)                        !ue
		   do m=1,mblocks
		      temp(m)%scalar = ue(m)%vector(ifield-17,:,:,:)
		   enddo
        case(21,22,23)                        !J
		   do m=1,mblocks
		      temp(m)%scalar = delJ(m)%vector(ifield-20,:,:,:)
		   enddo
        case(24,25,26)                        !delB
		   do m=1,mblocks
		      temp(m)%scalar = delB(m)%vector(ifield-23,:,:,:)
		   enddo
        case(27,28,29)                        !delE
		   do m=1,mblocks
		      temp(m)%scalar = delE(m)%vector(ifield-26,:,:,:)
		   enddo
        case(30,31,32)                        !S     ;j0 X dB
!           allocate(delJ(3,nz0:nz1,mgrid0:mgrid1))
!		   call get_j(bv_bar,delJ)
!		   call cross(delJ,Bv1,vec1)
		   do m=1,mblocks
		      temp(m)%scalar = cv_S(m)%vector(ifield-29,:,:,:)
		   enddo
        case(33,34,35)                        !dj X dB
!           allocate(delJ(3,nz0:nz1,mgrid0:mgrid1))
!		   call get_j(bv1,delJ)
!		   call cross(delJ,Bv1,vec1)
!           temp = vec1(ifield-32,:,:)
        case(36,37,38)                        !divPpe/ne

!           temp = divPpe(ifield-35,:,:)/max(dene,1.d-4)
!	       call dconvert(temp,plot)
        case(39)                        !divPpe/ne
!           temp = dUedt(3,:,:)

! ----------------------------------------------plasma informations ------------------
        case(41,42,43)        !ui(1-3,:)
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%vi(ifield-40)
			  if(MHD) temp(m)%scalar = block(m)%ionf(1,:,:,:)%v(ifield-40)
		   enddo
        case(44)              !Ni
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ni
			  if(simtype(1)<0) temp(m)%scalar = block(m)%ionf(1,:,:,:)%n
		   enddo
        case(45,46)           !ui(1-3,:),Ti
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ti(ifield-44)
			  if(simtype(1)<0) temp(m)%scalar = block(m)%ionf(1,:,:,:)%t
		   enddo
        case(47,48,49)        !ui(1-3,:) kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%vi(ifield-46)
		   enddo
        case(50)              !ui(1-3,:),Ni, kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ni
		   enddo
        case(51,52)           !ui(1-3,:),Ti, kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ti(ifield-50)
		   enddo

        case(53,54,55)        !ui(1-3,:) kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%vi(ifield-52)
		   enddo
        case(56)              !ui(1-3,:),Ni, kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%ni
		   enddo
        case(57,58)           !ui(1-3,:),Ti, kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%ti(ifield-56)
		   enddo
! -------------------------------------------------------------------------
        case(59,60,61)        !ue(1-3,:), kind=1
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ve(ifield-58)
		   enddo
        case(62)              !Ne
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ne
		   enddo
        case(63,64)           !Te
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%te(ifield-62)
		   enddo
        case(65,66,67)        !ue(1-3,:),kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ve(ifield-64)
		   enddo
        case(68)              !Ne --- 2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ne
		   enddo
        case(69,70)           !Te ----2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%te(ifield-68)
		   enddo
        case(71)              !beta_i
		   do m=1,mblocks
		      temp(m)%scalar = beta_i(m)%scalar
		   enddo
        case(72)              !beta_e
		   do m=1,mblocks
		      temp(m)%scalar = beta_e(m)%scalar
		   enddo
! special diagnosis for chorus wave particle -wave interaction Landau dumping or cyctrol 
!       case(81,82)              ! for electron Three kinds 81-86
!		   do m=1,mblocks
!		      temp(m)%scalar = block(m)%ini(1,:,:,:)%edotv(ifield-80)
!		   enddo
!        case(83,84)              
!		   do m=1,mblocks
!		      temp(m)%scalar = block(m)%ini(2,:,:,:)%edotv(ifield-82)
!		   enddo
        case(81)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(1,:,:,:)
		   enddo
        case(82)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(3,:,:,:)
		   enddo
        case(83)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(2,:,:,:)
		   enddo
        case(84)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(4,:,:,:)
		   enddo

        case(85,86)              ! for electron Three kinds 81-86
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%edotv(ifield-84)
		   enddo

        case(87,88)              ! for ions first two kinds
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%edotv(ifield-86)
		   enddo

        case(89,90)              
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%edotv(ifield-88)
		   enddo

        case(91,92)              
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%edotv(ifield-90)
		   enddo

! ---------------------------------
        case(97)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(5,:,:,:)
		   enddo
        case(98)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(7,:,:,:)
		   enddo
        case(99)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(6,:,:,:)
		   enddo
        case(100)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(8,:,:,:)
		   enddo

! -------------------------------------------------------------------

        case(101)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(1,:,:,:)
		   enddo
        case(102)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(2,:,:,:)
		   enddo
        case(103)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(3,:,:,:)
		   enddo
        case(104)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(4,:,:,:)
		   enddo
        case(105)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(5,:,:,:)
		   enddo
        case(106)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(6,:,:,:)
		   enddo

! check divB=0
        case(611)
           call div_vectorB(bv1,temp)
        case(612)
           call bcovconE(ev1,efld)
           call div_vectorE(efld,temp)
		   do m=1,mblocks
		      temp(m)%scalar = temp(m)%scalar - eta * (qion * (ni(m)%scalar+nif(m)%scalar) + &
		                        qelectron * (ne(m)%scalar +nef(m)%scalar ))
           enddo
	 end select

     if(box9%opt ==  0) call dconvert(temp,plot,box9)
     if(box9%opt ==  1) call Adconvert(temp,plot,box9)
!    if(mod(idiagbox,10) ==0) call dconvert(temp,plot,cuv)
     if(mype==0) call check( nf90_put_var(ncid, varid(5+iplot), plot) )
  enddo

  if(mype==0) call check( nf90_close(ncid) )

  deallocate(plot)
  call allocate_bscalar(temp,-1) 
  call allocate_bscalar(beta_i,-1) 
  call allocate_bscalar(beta_e,-1) 
  call allocate_bvector(bfld,-1) 
  call allocate_bvector(efld,-1) 
  call allocate_bvector(ui,  -1) 
  call allocate_bvector(ue,  -1) 
  call allocate_bvector(delJ,-1) 
  call allocate_bvector(dumm,-1) 
  call allocate_bvector(av1, -1) 
  call allocate_bvector(delB,-1) 
  call allocate_bvector(delE,-1) 

  if(ifcheck(8) > 0) call allocate_bvector(J_be,-1) 
  if(ifcheck(9) > 0) then
     call allocate_bvector(v_s,-3) 
     call allocate_bvector(cv_s,-3) 
  endif
  if(ifcheck(10) > 0) call allocate_bvector(PPPe,-1) 


  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine diagnos_field_nc
! -------------------------------------------------------------------
! for global simulhation 
subroutine diagnos_field_nc1file(box9,iflag)
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  implicit none
  include 'mpif.h'
  type(diag_type)   :: box9
  real,    pointer,dimension(:,:,:) :: plot
  real(p2),pointer,dimension(:,:,:) :: test1
  real(p2),pointer,dimension(:,:,:,:) :: vec1,vec2,test3
  type(bvector_type), dimension(:), pointer    :: J_BE,PPPe
  real(p2) :: param(80)
  integer :: ic,i,iplot,ifield,m,ifchecki,ifchecke,n,ifcheck(30),iflag
  type(bscalar_type), dimension(:), pointer    :: temp,beta_i,beta_e
  type(bvector_type), dimension(:), pointer    :: bfld,efld,ui,ue,delJ,dumm,av1,delB,delE

  character*6 fldname(110)
  character*2 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(5),varid(50),dimids(4),icount(4),istart(4),icount1(1),istart1(1)
  logical    file_exist

  save ic
  data ic/0/

  interface 
      subroutine allocate_bscalar(c,m)
        use constants
		use grid_data
	    type(bscalar_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bscalar

      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector

      subroutine dconvert(a,b,c)
		use grid_data
	    type(bscalar_type), dimension(:), pointer    :: a
        real,    pointer,dimension(:,:,:)            :: b
	    type(diag_type)                              :: c
      end subroutine dconvert

    subroutine dconvert_pqw(a,b)
       use constants
	   use grid_data
	   type(bscalar_type), dimension(:), pointer    :: a
       real,pointer,dimension(:,:,:) :: b
    end subroutine dconvert_pqw

    subroutine cconvert(a,b)
       use constants
      real(p2),pointer,dimension(:,:,:) :: a,b
    end subroutine cconvert
    subroutine cross(a,b,c)
       use constants
       real(p2),pointer,dimension(:,:,:) :: a,b,c
    end subroutine cross
    subroutine getue(a)
       use constants
      real(p2),pointer,dimension(:,:,:) :: a
    end subroutine getue

    subroutine filter_vector(a,m)
       use constants
      real(p2),pointer,dimension(:,:,:) :: a
	  integer :: m
    end subroutine filter_vector

    subroutine becellecenter(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine becellecenter

    subroutine bbcellecenter(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bbcellecenter

    subroutine bconcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bconcarE

    subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bconcarB

    subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bcovcarE

    subroutine gather_scalar(a,b)
       use constants
	   use grid_data
	   type(bscalar_type), dimension(:), pointer    :: a
       real(p2),pointer,dimension(:,:,:) :: b
    end subroutine gather_scalar
    subroutine get_j(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine get_j

    subroutine get_A(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine get_A
    subroutine field3(a)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
    end subroutine field3

    subroutine get_JBE(a)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
    end subroutine get_JBE

    subroutine get_PPe(a)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
    end subroutine get_PPe

    subroutine get_vdevdb(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
	   type(particle_type), dimension(:), pointer    :: b
    end subroutine get_vdevdb

    subroutine getUT(m,a)
       use constants
	   use grid_data
	   type(particle_type), dimension(:), pointer    :: a
    end subroutine getUT

      subroutine Adconvert(a,b,c)
		use grid_data
	    type(bscalar_type), dimension(:), pointer    :: a
        real,    pointer,dimension(:,:,:)            :: b
	    type(diag_type)                              :: c
      end subroutine Adconvert


  end interface


  fldname = (/"Bx    ","By    ","Bz    ","Ni    ","Ne    ","phi   ","A1x   ","A1y   ","A1z   ","Uepar ", & !10
              "Avt   ","Ex    ","Ey    ","Ez    ","Upx   ","Upy   ","Upz   ","Uex   ","Uey   ","Uez   ", & !20
              "Jx    ","Jy    ","Jz    ","B1x   ","B1y   ","B1z   ","E1x   ","E1y   ","E1z   ","EcB1  ", & !30
			  "EcB2  ","EcB3  ","EcB4  ","JcB8  ","JcB9  ","divP1 ","divP2 ","divP3 ","duet1 ","duet2 ", & !40
			  "ui1x  ","ui1y  ","ui1z  ","Ni(1) ","Ti1par","Ti1per","ui2x  ","ui2y  ","ui2z  ","Ni(2) ", & !50
			  "Ti2par","Ti2per","ui3x  ","ui3y  ","ui3z  ","Ni(3) ","T31par","Ti3per","ue1x  ","ue1y  ", & !60
			  "ue1z  ","Ne(1) ","Te1par","Te1per","ue2x  ","ue2y  ","ue2z  ","Ne(2) ","Te2par","Te2per", & !70
			  "betai ","betae ","divB  ","divE  ","ue2x  ","ue2y  ","ue2z  ","Ne(2) ","Te2par","Te2per", & !80
			  "edv11 ","edv12 ","edv21 ","edv22 ","edv31 ","edv32 ","edv11i","edv12i","edv21i","edv22i", & !90
			  "edv31i","edv32i","edv21 ","edv22 ","edv31 ","edv32 ","JB1   ","JE1   ","JB2   ","JE2   ", & !100
			  "PPe11 ","PPe22 ","PPe33 ","PPe12 ","PPe13 ","PPe23 ","JB11  ","JE11  ","JB21  ","JE21  "/)  !110

!  fldname = (/"Bx","By","Bz","Ni","Ne","phi","A1x","A1y","A1z","Uepr",&
!              "Avt","Ex","Ey","Ez","Upx","Upy","Upz","Uex","Uey","Uez", &
!              "Jx","Jy","Jz","B1x","B1y","B1z","JcB1","JcB2","JcB3","JcB4",&
!			  "JcB5","JcB6","JcB7","JcB8","JcB9","divP1","divP2","divP3","duet1","duet2", &
!			  "ui1x","ui1y","ui1z","Ti1par","Ti1per","ui2x","ui2y","ui2z","Ti2par","Ti2per",&
!			  "ue1x","ue1y","ue1z","Te1par","Te1per","ue2x","ue2y","ue2z","Te2par","Te2per",&
!			  "divB","divE"/)

  if(iflag ==0 ) then
     idiagf  = idiagf + 1
  endif

  if(iflag ==1 ) then
     idiagf1 = idiagf1 + 1
  endif	  

  if(iflag ==2 ) then
     idiagf2 = idiagf2 + 1
  endif	  
! --------------------------------------------------------------
  if(iflag ==0 )  ic  = idiagf  !istep/ndiagf  + 1   
  if(iflag ==1 )  ic  = idiagf1 !istep/ndiagf1 + 1   
  if(iflag ==2 )  ic  = idiagf2 !istep/ndiagf2 + 1   
 
  if(ic > 20000) then
     if(iflag ==0 ) ic0_save(5)  = ic0_save(5) + 1
     if(iflag ==1 ) ic0_save(6)  = ic0_save(6) + 1
     if(iflag ==2 ) ic0_save(7)  = ic0_save(7) + 1
  endif
  
  if(iflag == 0) write(ct,'(I2.2)')ic0_save(5) + 1 
  if(iflag == 1) write(ct,'(I2.2)')ic0_save(6) + 1 
  if(iflag == 2) write(ct,'(I2.2)')ic0_save(7) + 1 
! --------------------------------------------------------------
  if(iflag == 0) filename = 'field1file'//ct//'.nc'
  if(iflag == 1) filename = 'fieldfine1file'//ct//'.nc'
  if(iflag == 2) filename = 'fieldbox1file'//ct//'.nc'

  inquire(file=filename,exist=file_exist)
  if(.not. file_exist .or. (ic == 1 .and. irun ==0)) then
      if(iflag ==0 )  idiagf  = 1  
      if(iflag ==1 )  idiagf1 = 1  
      if(iflag ==2 )  idiagf2 = 1  
  endif

  if(iflag ==0 )  ic  = idiagf  !istep/ndiagf  + 1   
  if(iflag ==1 )  ic  = idiagf1 !istep/ndiagf1 + 1   
  if(iflag ==2 )  ic  = idiagf2 !istep/ndiagf2 + 1   

  icount1 = (/1/)
  istart1 = (/ic/)

  param     = 0.
  param(1)  = stime  ! has been scaled to omega_i
  param(2)  = aelectron
  param(3)  = aion
  param(4)  = br0
  param(5)  = bt0
  param(6)  = bz0
  param(7)  = betae
  param(8)  = deltax !box9%dxd
  param(9)  = deltay !box9%dyd
  param(10) = deltaz !box9%dzd
  param(11) = kinds
! real unit ----
  param(12) = ubfield   !ulength/rhos  !because xdiag is in unit of rhos
  param(13) = uefield   !utime  /alfventime
  param(14) = udensity  !uspeed
  param(15) = uspeed    !uspeed
  param(16) = utime     !uspeed
  param(17) = utemperature
  param(18) = uj0
  param(19) = uphi
  param(20) = ulength*rhos
  param(21) = cva
  param(22) = lamda_i/rhos
  param(23) = lamda_e/rhos
  param(24) = uRe/rhos
  param(25) = lamda_de/rhos
  param(29) = va
  param(30) = alpha

! some key parameters
  param(41) = omega_lh/omega_i
  param(75-box9%nfield+1:75)= DiagQuantity(1:box9%nfield) 
  param(77) = simtype(2)
  param(78) = simtype(1)
  param(79) = coordinate
  param(80) = box9%nfield   ! number of field quantity are recorded



  istart = (/ 1, 1, 1, ic /)
  if(mype==0 ) then
!     if(mod(idiagbox9,10) ==1) then
        icount = (/box9%nxd, box9%nyd, box9%nzd, 1 /)
!       filename = 'field.nc'
!        if(iflag == 0) filename = 'field1file'//ct//'.nc'
!        if(iflag == 1) filename = 'fieldfine1file'//ct//'.nc'
!        if(iflag == 2) filename = 'fieldbox1file'//ct//'.nc'
!       Create the file. 
        if(ic == 1) then
        call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
        call check( nf90_def_dim(ncid, 'param', 80, dimid(1)) )
        call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(1), varid(1)) )
        call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

        call check( nf90_def_dim(ncid, 'X', box9%nxd, dimid(2)) )
        call check( nf90_def_var(ncid, 'X', NF90_REAL, dimid(2), varid(2)) )
        call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

        call check( nf90_def_dim(ncid, 'Y', box9%nyd, dimid(3)) )
        call check( nf90_def_var(ncid, 'Y', NF90_REAL, dimid(3), varid(3)) )
        call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

        call check( nf90_def_dim(ncid, 'Z', box9%nzd, dimid(4)) )
        call check( nf90_def_var(ncid, 'Z', NF90_REAL, dimid(4), varid(4)) )
        call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

        call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
        call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )


        dimids=(/dimid(2),dimid(3),dimid(4),dimid(5)/)
        do iplot = 1, box9%nfield
           ifield = box9%DiagQuantity(iplot)
           call check( nf90_def_var(ncid, fldname(ifield), NF90_REAL, dimids, varid(5+iplot)) )
           call check( nf90_put_att(ncid, varid(5+iplot), UNITS, 'non') )
        enddo

        call check( nf90_enddef(ncid) )

        call check( nf90_put_var(ncid, varid(1), param) )
        call check( nf90_put_var(ncid, varid(2), box9%xdiag) )
        call check( nf90_put_var(ncid, varid(3), box9%ydiag) )
        call check( nf90_put_var(ncid, varid(4), box9%zdiag) )
        call check( nf90_put_var(ncid, varid(5), stime, start = istart1) )
        
		else
!       ic > 1 open exist netCDF file
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime,start = istart1) )

          do iplot = 1, box9%nfield
             ifield = box9%DiagQuantity(iplot)
             call check( nf90_inq_varid(ncid, fldname(ifield), varid(5+iplot)) )
          enddo
        endif
!	 endif
  endif


  allocate(plot(box9%nxd,box9%nyd,box9%nzd))

  call allocate_bscalar(temp,1) 
  call allocate_bscalar(beta_i,1) 
  call allocate_bscalar(beta_e,1) 
  call allocate_bvector(bfld,3) 
  call allocate_bvector(efld,3) 
  call allocate_bvector(ui,  3) 
  call allocate_bvector(ue,  3) 
  call allocate_bvector(delJ,3) 
  call allocate_bvector(dumm,3) 
  call allocate_bvector(av1, 3) 
  call allocate_bvector(delB,3) 
  call allocate_bvector(delE,3) 

! (1) b-field -------
  if(mod(box9%idiagbox,10) ==1  .or. mod(box9%idiagbox,10) ==3) then
     do m=1,mblocks
        dumm(m)%vector = bv1(m)%vector  + bv0(m)%vector 
     enddo
	 call bconcarB(dumm,bfld)
  else
     do m=1,mblocks
	    do n=1,3
           dumm(m)%vector(n,:,:,:) = (bv1(m)%vector(n,:,:,:)  + bv0(m)%vector(n,:,:,:))* &
		                             block(m)%node%h(4,n,:,:,:) 
        enddo
     enddo
	 call bbcellecenter(dumm,bfld)
  endif
! (2) e-field -------
  if(mod(box9%idiagbox,10) ==1) then
     do m=1,mblocks
        dumm(m)%vector = ev1(m)%vector 
     enddo
	 call bcovcarE(dumm,efld)
  else
     do m=1,mblocks
	    do n=1,3
           dumm(m)%vector(n,:,:,:) =  ev1(m)%vector(n,:,:,:)/block(m)%node%h(3,n,:,:,:)   
        enddo
     enddo
     call becellecenter(dumm,efld)
  endif
! (3) ui-field -------
  if(mod(box9%idiagbox,10) ==1) then
     do m=1,mblocks
        dumm(m)%vector = ji(m)%vector 
     enddo
	 call bconcarE(dumm,ui)
  else
     do m=1,mblocks
	    do n=1,3
           dumm(m)%vector(n,:,:,:) = ji(m)%vector(n,:,:,:)*block(m)%node%h(3,n,:,:,:)
        enddo
     enddo
	 call becellecenter(dumm,ui)
  endif

  do m=1,mblocks
	 do n=1,3
        if(.not. df .and. simtype(1) /= 6 .and. simtype(1) /= 7 ) &
		    ui(m)%vector(n,:,:,:) = ui(m)%vector(n,:,:,:) /max(ni(m)%scalar,1.e-4_p2)
        if(df)ui(m)%vector(n,:,:,:) = ui(m)%vector(n,:,:,:) /max(ni0(m)%scalar,1.e-4_p2)
	 enddo
  enddo
! (4) ue-field --------------
  if(simtype(1) == 6 .or. simtype(1) == 7 .or. simtype(1) == 66 .or. simtype(1) == 67) then
!     do m=1,mblocks
!	    ue(m)%vector(1,:,:,:) = block(m)%ionf(1,:,:,:)%v(1)
!	    ue(m)%vector(2,:,:,:) = block(m)%ionf(1,:,:,:)%v(2)
!	    ue(m)%vector(3,:,:,:) = block(m)%ionf(1,:,:,:)%v(3)
!	 enddo
  else
     if(mod(box9%idiagbox,10) ==1) then
        do m=1,mblocks
           dumm(m)%vector = je(m)%vector/qelectron 
        enddo
	    call bcovcarE(dumm,ue)
     else
        do m=1,mblocks
	       do n=1,3
              dumm(m)%vector(n,:,:,:) = je(m)%vector(n,:,:,:)/qelectron !/block(m)%node%h(1,n,:,:,:)
           enddo
        enddo
	    call becellecenter(dumm,ue)
     endif

     do m=1,mblocks
	    do n=1,3
           if(.not. df)ue(m)%vector(n,:,:,:) = ue(m)%vector(n,:,:,:) /max(ne(m)%scalar,1.e-4_p2)
           if(df)ue(m)%vector(n,:,:,:) = ue(m)%vector(n,:,:,:) /max(ne0(m)%scalar,1.e-4_p2)
	    enddo
     enddo
  endif

  ifchecki = 0
  ifchecke = 0
  ifcheck  = 0
  do iplot = 1,box9%nfield
     ifield = box9%DiagQuantity(iplot)
	 if(ifield >= 40 .and. ifield <= 58 ) then
	    ifchecki =   ifchecki + 1
	    if(ifchecki ==1) call getUT(0,ions)
     endif
	 if((ifield >= 59 .and. ifield<=70) .or.(ifield >= 81 .and. ifield<=86)  ) then
	    ifchecke =   ifchecke + 1
	    if(ifchecke ==1) call getUT(0,eles)
     endif

	 if(ifield  == 21 .or. ifield  == 22 .or. ifield  == 23 ) then  ! delta_j
	    ifcheck(3) = ifcheck(3) +1
		if(ifcheck(3) ==1) then
		   if(simtype(1) >=6) then   ! hybrid simulation j= curlB
  	          do m=1,mblocks
		         dumm(m)%vector    = bv1(m)%vector  + bv0(m)%vector 
		      enddo
		      call get_j(dumm,delJ)
		   else     ! otherwise, J=Ji+Je
              if(mod(box9%idiagbox,10) ==1) then
                 do m=1,mblocks
                    dumm(m)%vector = ji(m)%vector + je(m)%vector 
                 enddo
	             call bconcarE(dumm,delJ)
              else
                 do m=1,mblocks
	             do n=1,3
                    dumm(m)%vector(n,:,:,:) = (ji(m)%vector(n,:,:,:)+je(m)%vector(n,:,:,:)) *block(m)%node%h(3,n,:,:,:)
                 enddo
                 enddo
	             call becellecenter(dumm,delJ)
              endif
		   endif
        endif
	 endif

	 if(ifield  == 7 .or. ifield  == 8 .or. ifield  == 9 ) then  ! delta_j
	    ifcheck(4) = ifcheck(4) +1
		if(ifcheck(4) ==1) then
  	       do m=1,mblocks
		      dumm(m)%vector    = bv1(m)%vector  + bv0(m)%vector 
		   enddo
		   call get_A(dumm,Av1)
        endif
	 endif
	 if(ifield  == 24 .or. ifield  == 25 .or. ifield  == 26 ) then  ! delta_B
	    ifcheck(5) = ifcheck(5) +1
		if(ifcheck(5) ==1) then
           if(abs(coordinate)  ==1 .or.abs(coordinate)  ==2 .or. abs(coordinate)  ==3 .or. abs(coordinate)  ==8 ) then
		      do m=1,mblocks
	             do n =1,3
                    dumm(m)%vector(n,:,:,:) = bv1(m)%vector(n,:,:,:) * block(m)%node%h(4,n,:,:,:)
		         enddo
			  enddo
	          call bbcellecenter(dumm,delB)

		   else
  	          call field3(delB)
           endif
        endif
	 endif

	 if(ifield  == 27 .or. ifield  == 28 .or. ifield  == 29 ) then  ! delta_E
	    ifcheck(6) = ifcheck(6) +1
		if(ifcheck(6) ==1) then
           if(abs(coordinate)  ==1 .or.abs(coordinate)  ==2 .or. abs(coordinate)  ==3 .or. abs(coordinate)  ==8 ) then
		      do m=1,mblocks
	             do n =1,3
                    dumm(m)%vector(n,:,:,:)  = ev1(m)%vector(n,:,:,:) / block(m)%node%h(3,n,:,:,:)
		         enddo
			  enddo
	          call becellecenter(dumm,delE)
		   else

           endif
        endif
	 endif

	 if(ifield  == 71) then  ! beta_i
        do m=1,mblocks
		   beta_i(m)%scalar = 0.
		   do n=1,kinds
		      beta_i(m)%scalar = beta_i(m)%scalar + block(m)%ini(n,:,:,:)%ni * &
			                     (block(m)%ini(n,:,:,:)%Ti(1)+block(m)%ini(n,:,:,:)%Ti(2)*2)/3.
		   enddo
		   temp(m)%scalar = bfld(m)%vector(1,:,:,:)**2+bfld(m)%vector(2,:,:,:)**2+bfld(m)%vector(3,:,:,:)**2
        enddo
        do m=1,mblocks
		   beta_i(m)%scalar = beta_i(m)%scalar/max(temp(m)%scalar,1.e-4_p2)*alpha
        enddo
	 endif

	 if((ifield  >= 97 .and. ifield  <= 100) .or. (ifield  >= 81 .and. ifield  <= 84) ) then  ! J_BE will be saved
	    ifcheck(8) = ifcheck(8) +1
		if(ifcheck(8) ==1 .and. dipole ) then
           call allocate_bvector(J_BE,8) 
		   call get_vdEvdB(j_be,eles)
        endif
	 endif

	 if((ifield  >= 101 .and. ifield  <= 106)  ) then  ! J_BE will be saved
	    ifcheck(10) = ifcheck(10) +1
		if(ifcheck(10) ==1 .and. dipole ) then
           call allocate_bvector(PPPe,6) 
		   call get_PPe(PPPe)
        endif
	 endif

  enddo

  do iplot = 1, box9%nfield
     ifield = box9%DiagQuantity(iplot)
     select case(ifield)
        case(1,2,3)                           !bx,by,bz
		   do m=1,mblocks
		      temp(m)%scalar = bfld(m)%vector(ifield,:,:,:)
		   enddo
        case(4)                               !ni
		   do m=1,mblocks
		      temp(m)%scalar = ni(m)%scalar
		   enddo
        case(5)                               !ne
		   do m=1,mblocks
		      temp(m)%scalar = ne(m)%scalar
		      if(simtype(1) == 6 .or. simtype(1) == 7 .or. simtype(1) == 66 .or. simtype(1) == 67) &
              temp(m)%scalar = nif(m)%scalar
		   enddo
        case(6)                               !phi
		   do m=1,mblocks
		      temp(m)%scalar = phi(m)%scalar
		   enddo
        case(7,8,9)                           !a1
!          temp = av1(ifield-6,:,:,:)
        case(10)                              !Ue||
!          temp = uepar     !Jpe(3,:,:)
        case(11)                              !az
!          call cconvert(avt,vec1)
!          temp = vec1(3,:,:,:)
        case(12,13,14)                        !E
		   do m=1,mblocks
		      temp(m)%scalar = efld(m)%vector(ifield-11,:,:,:)
		   enddo
        case(15,16,17)                        !ui
		   do m=1,mblocks
		      temp(m)%scalar = ui(m)%vector(ifield-14,:,:,:)
		   enddo
        case(18,19,20)                        !ue
		   do m=1,mblocks
		      temp(m)%scalar = ue(m)%vector(ifield-17,:,:,:)
		   enddo
        case(21,22,23)                        !J
		   do m=1,mblocks
		      temp(m)%scalar = delJ(m)%vector(ifield-20,:,:,:)
		   enddo
        case(24,25,26)                        !delB
		   do m=1,mblocks
		      temp(m)%scalar = delB(m)%vector(ifield-23,:,:,:)
		   enddo
        case(27,28,29)                        !delE
		   do m=1,mblocks
		      temp(m)%scalar = delE(m)%vector(ifield-26,:,:,:)
		   enddo
        case(30,31,32)                        !j0 X dB
!           allocate(delJ(3,nz0:nz1,mgrid0:mgrid1))
!		   call get_j(bv_bar,delJ)
!		   call cross(delJ,Bv1,vec1)
!           temp = vec1(ifield-29,:,:)
        case(33,34,35)                        !dj X dB
!           allocate(delJ(3,nz0:nz1,mgrid0:mgrid1))
!		   call get_j(bv1,delJ)
!		   call cross(delJ,Bv1,vec1)
!           temp = vec1(ifield-32,:,:)
        case(36,37,38)                        !divPpe/ne

!           temp = divPpe(ifield-35,:,:)/max(dene,1.d-4)
!	       call dconvert(temp,plot)
        case(39)                        !divPpe/ne
!           temp = dUedt(3,:,:)

! ----------------------------------------------plasma informations ------------------
        case(41,42,43)        !ui(1-3,:)
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%vi(ifield-40)
		   enddo
        case(44)              !Ni
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ni
		   enddo
        case(45,46)           !ui(1-3,:),Ti
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ti(ifield-44)
		   enddo
        case(47,48,49)        !ui(1-3,:) kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%vi(ifield-46)
		   enddo
        case(50)              !ui(1-3,:),Ni, kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ni
		   enddo
        case(51,52)           !ui(1-3,:),Ti, kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ti(ifield-50)
		   enddo

        case(53,54,55)        !ui(1-3,:) kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%vi(ifield-52)
		   enddo
        case(56)              !ui(1-3,:),Ni, kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%ni
		   enddo
        case(57,58)           !ui(1-3,:),Ti, kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%ti(ifield-56)
		   enddo
! -------------------------------------------------------------------------
        case(59,60,61)        !ue(1-3,:), kind=1
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ve(ifield-58)
		   enddo
        case(62)              !Ne
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ne
		   enddo
        case(63,64)           !Te
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%te(ifield-62)
		   enddo
        case(65,66,67)        !ue(1-3,:),kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ve(ifield-64)
		   enddo
        case(68)              !Ne --- 2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ne
		   enddo
        case(69,70)           !Te ----2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%te(ifield-68)
		   enddo
        case(71)              !beta_i
		   do m=1,mblocks
		      temp(m)%scalar = beta_i(m)%scalar
		   enddo
        case(72)              !beta_e
		   do m=1,mblocks
		      temp(m)%scalar = beta_e(m)%scalar
		   enddo
! special diagnosis for chorus wave particle -wave interaction Landau dumping or cyctrol 
        case(81)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(1,:,:,:)
		   enddo
        case(82)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(3,:,:,:)
		   enddo
        case(83)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(2,:,:,:)
		   enddo
        case(84)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(4,:,:,:)
		   enddo


        case(85,86)              ! for electron Three kinds 81-86
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%edotv(ifield-84)
		   enddo

! ---------------------------------
        case(97)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(5,:,:,:)
		   enddo
        case(98)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(7,:,:,:)
		   enddo
        case(99)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(6,:,:,:)
		   enddo
        case(100)              
		   do m=1,mblocks
		      temp(m)%scalar = J_be(m)%vector(8,:,:,:)
		   enddo

! -------------------------------------------------------------------

        case(101)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(1,:,:,:)
		   enddo
        case(102)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(2,:,:,:)
		   enddo
        case(103)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(3,:,:,:)
		   enddo
        case(104)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(4,:,:,:)
		   enddo
        case(105)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(5,:,:,:)
		   enddo
        case(106)              
		   do m=1,mblocks
		      temp(m)%scalar = PPPe(m)%vector(6,:,:,:)
		   enddo

! check divB=0
        case(611)
           call div_vectorB(bv1,temp)
        case(612)
           call bcovconE(ev1,efld)
           call div_vectorE(efld,temp)
		   do m=1,mblocks
		      temp(m)%scalar = temp(m)%scalar - eta * (qion * (ni(m)%scalar+nif(m)%scalar) + &
		                        qelectron * (ne(m)%scalar +nef(m)%scalar ))
           enddo
	 end select

     if(box9%opt ==  0) call dconvert(temp,plot,box9)
     if(box9%opt ==  1) call Adconvert(temp,plot,box9)
     if(mype==0) call check( nf90_put_var(ncid, varid(5+iplot), plot, start = istart, &
                              count = icount) )

  enddo

  if(mype==0) call check( nf90_close(ncid) )

  deallocate(plot)
  call allocate_bscalar(temp,-1) 
  call allocate_bscalar(beta_i,-1) 
  call allocate_bscalar(beta_e,-1) 
  call allocate_bvector(bfld,-1) 
  call allocate_bvector(efld,-1) 
  call allocate_bvector(ui,  -1) 
  call allocate_bvector(ue,  -1) 
  call allocate_bvector(delJ,-1) 
  call allocate_bvector(dumm,-1) 
  call allocate_bvector(av1, -1) 
  call allocate_bvector(delB,-1) 
  call allocate_bvector(delE,-1) 
  if(ifcheck(8) > 0) call allocate_bvector(J_be,-1) 
  if(ifcheck(10) > 0) call allocate_bvector(PPPe,-1) 
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine diagnos_field_nc1file
! -------------------------------------------------------------------

subroutine output4mc_nc
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use global_case
  use output_mc
  implicit none
  include 'mpif.h'
! real,pointer,dimension(:,:,:) :: en,ux,uy,uz,tpar,tper
  type(bscalar_type), dimension(:), pointer     :: temp,temp1
  type(bvector_type), dimension(:), pointer     :: Tpar,Tper
  real,pointer,dimension(:,:,:) :: plot3
  real,pointer,dimension(:,:,:,:) :: ff,fft,ffp
  real,pointer,dimension(:) :: coorx,coory,coorz
  character*5 ct
  character*10 fldname(10)
  integer :: iplot,ifield
  integer :: m,n,i,j,k,index,L,imc1,imc2,imc3,imc4,jmc1,jmc2,jmc3,jmc4,IERROR
  real ::    param(50),factor,x,y,z,vx,vy,vz,r,phi0_,ee,angle,bx,by,bz,vtot,vpar,dmc1,dmc2,dmc3,dmc4, &
             amass
  real*8  :: xyz(3),p,q,w
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(20),varid(30),dimids4(4),dimids2(2),dimids3(3)

  interface 
      subroutine allocate_bscalar(c,m)
        use constants
		use grid_data
	    type(bscalar_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bscalar
      subroutine allocate_bvector(c,m)
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine dconvert(a,b,c)
		use grid_data
	    type(bscalar_type), dimension(:), pointer    :: a
        real,    pointer,dimension(:,:,:)            :: b
	    type(diag_type)                              :: c
      end subroutine dconvert
      subroutine bcell2ecenter(a,b)
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcell2ecenter

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw

  end interface

! if(ionsolver < 0) return
! new output for MC
! save ion information of f(E,alpha,equator) where E-energy and cos(alpha)-pich angle;   and density 
  idiagfRmc = idiagfRmc + 1
  write(ct,'(I5.5)')IDIAGfRmc 
  write(ct,'(I5.5)')int(stime + 2.*dt)
! for netCDF file generation

! fldname = (/"PSD","num_par","density","Tpar","Tper","Bx","By","Bz","potential","Jpar"/)
  fldname = (/"PSD      ",&
              "num_par  ",&
			  "density  ",&
			  "Tpar     ",&
			  "Tper     ",&
			  "Bx       ",&
			  "By       ",&
			  "Bz       ",&
			  "potential",&
			  "Jpar     "  /)

  if(mype==0) then
     filename = "eqdata"//ct//".nc"
!    Create the file. 
     call check( nf90_create(FILENAME, nf90_clobber, ncid) )

! define dimensions
     call check( nf90_def_dim(ncid, 'time', 1, dimid(1)) )
     call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(1), varid(1)) )
     call check( nf90_put_att(ncid,  varid(1), UNITS, "s") )

     call check( nf90_def_dim(ncid, 'logE', n_E, dimid(2)) )
     call check( nf90_def_var(ncid, 'logE', NF90_REAL, dimid(2), varid(2)) )
     call check( nf90_put_att(ncid,  varid(2), UNITS, "log(eV)") )

     call check( nf90_def_dim(ncid, 'alpha', n_alpha, dimid(3)) )
     call check( nf90_def_var(ncid, 'alpha', NF90_REAL, dimid(3), varid(3)) )
     call check( nf90_put_att(ncid,  varid(3), UNITS, "rad") )

     call check( nf90_def_dim(ncid, 'r', n_r, dimid(4)) )
     call check( nf90_def_var(ncid, 'r', NF90_REAL, dimid(4), varid(4)) )
     call check( nf90_put_att(ncid,  varid(4), UNITS, "Re") )

     call check( nf90_def_dim(ncid, 'mlt', n_phi, dimid(5)) )
     call check( nf90_def_var(ncid, 'mlt', NF90_REAL, dimid(5), varid(5)) )
     call check( nf90_put_att(ncid,  varid(5), UNITS, "rad") )

     call check( nf90_def_dim(ncid, 'x', mcb%nxd, dimid(6)) )
     call check( nf90_def_var(ncid, 'x', NF90_REAL, dimid(6), varid(6)) )
     call check( nf90_put_att(ncid,  varid(6), UNITS, "Re") )

     call check( nf90_def_dim(ncid, 'y', mcb%nyd, dimid(7)) )
     call check( nf90_def_var(ncid, 'y', NF90_REAL, dimid(7), varid(7)) )
     call check( nf90_put_att(ncid,  varid(7), UNITS, "Re") )

     call check( nf90_def_dim(ncid, 'z', mcb%nzd, dimid(8)) )
     call check( nf90_def_var(ncid, 'z', NF90_REAL, dimid(8), varid(8)) )
     call check( nf90_put_att(ncid,  varid(8), UNITS, "Re") )

     call check( nf90_def_dim(ncid, 'latitude', mion, dimid(9)) )
     call check( nf90_def_var(ncid, 'latitude', NF90_REAL, dimid(9), varid(9)) )
     call check( nf90_put_att(ncid,  varid(9), UNITS, "rad") )

     call check( nf90_def_dim(ncid, 'longitude', nion, dimid(10)) )
     call check( nf90_def_var(ncid, 'longitude', NF90_REAL, dimid(10), varid(10)) )
     call check( nf90_put_att(ncid,  varid(10), UNITS, "rad") )

     allocate (COORX(mion),COORY(nion))
     do m = 1, mion
        coorx(m) = (m-1.)*dtheta
     enddo
     do m = 1, nion
        coory(m) = (m-1.)*dphi
     enddo


     dimids4=(/dimid(2),dimid(3),dimid(4),dimid(5)/)

     call check( nf90_def_var(ncid, fldname(1), NF90_REAL, dimids4, varid(11)) )  !psd
     call check( nf90_put_att(ncid, varid(11), UNITS, "1/cm^3log(E)rd") )
     call check( nf90_def_var(ncid, fldname(2), NF90_REAL, dimids4, varid(12)) )  !particles number
     call check( nf90_put_att(ncid, varid(12), UNITS, "#") )

     dimids2=(/dimid(4),dimid(5)/)
     call check( nf90_def_var(ncid, fldname(3), NF90_REAL, dimids2, varid(13)) )  !psd
     call check( nf90_put_att(ncid, varid(13), UNITS, "1/cm^3") )
     call check( nf90_def_var(ncid, fldname(4), NF90_REAL, dimids2, varid(14)) )  !psd
     call check( nf90_put_att(ncid, varid(14), UNITS, "eV") )
     call check( nf90_def_var(ncid, fldname(5), NF90_REAL, dimids2, varid(15)) )  !psd
     call check( nf90_put_att(ncid, varid(15), UNITS, "eV") )

     dimids3=(/dimid(6),dimid(7),dimid(8)/)
	 do m=1,3
        call check( nf90_def_var(ncid, fldname(5+m), NF90_REAL, dimids3, varid(15+m)) )  !psd
        call check( nf90_put_att(ncid, varid(15+m), UNITS, "nT") )
     enddo

     dimids2=(/dimid(9),dimid(10)/)
     call check( nf90_def_var(ncid, fldname(9), NF90_REAL, dimids2, varid(19)) )  !psd
     call check( nf90_put_att(ncid, varid(19), UNITS, "V") )

     call check( nf90_def_var(ncid, fldname(10), NF90_REAL, dimids2, varid(20)) )  !psd
     call check( nf90_put_att(ncid, varid(20), UNITS, "muA/m") )


     call check( nf90_enddef(ncid) )

     call check( nf90_put_var(ncid, varid(1), stime*utime) )
     call check( nf90_put_var(ncid, varid(2), array_e) )
     call check( nf90_put_var(ncid, varid(3), array_alpha) )
     call check( nf90_put_var(ncid, varid(4), array_r) )
     call check( nf90_put_var(ncid, varid(5), array_phi) )
     call check( nf90_put_var(ncid, varid(6), mcb%xdiag) )
     call check( nf90_put_var(ncid, varid(7), mcb%ydiag) )
     call check( nf90_put_var(ncid, varid(8), mcb%zdiag) )
     call check( nf90_put_var(ncid, varid(9), coorx) )
     call check( nf90_put_var(ncid, varid(10),coory) )
	 deallocate(coorx,coory)
!	 call check( nf90_close(ncid) )
  endif

! call getplasmadistribution(2)
  allocate(ff(n_e,n_alpha,n_r,n_phi),fft(n_e,n_alpha,n_r,n_phi),ffp(n_e,n_alpha,n_r,n_phi))
  ff = 0.
  do m=1,mblocks
     dO  L = 1, block(m)%mi
         index  =   block(m)%ion(L)%kind 
         if(index <= kinds) then
            Factor= block(m)%ion(L)%w * fraci(index) * qions(index)
            amass = aion*mions(index)

	        p     = block(m)%ion(L)%p(1)
	        q     = block(m)%ion(L)%p(2) 
	        w     = block(m)%ion(L)%p(3)
			i     = int(p) + 1
			j     = int(q) + 1
			k     = int(w) + 1
		    xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
            x     = xyz(1)
			y     = xyz(2)
			z     = xyz(3)
		    r   = sqrt(x**2+y**2+z**2)
            if(r >= r_min .and. r<r_max .and. abs(z)<= 1.) then
		       phi0_ = acos(max(min(x/max(sqrt(x**2+y**2),1.e-4),1.),-1.))
		       if(y < 0.) phi0_ = pi*2.-phi0_
	 	       bx    = bcar(m)%vector(1,i,j,k)
		       by    = bcar(m)%vector(2,i,j,k)
		       bz    = bcar(m)%vector(3,i,j,k)
	           vx    = block(m)%ion(L)%v(1)/ amass
	           vy    = block(m)%ion(L)%v(2)/ amass
	           vz    = block(m)%ion(L)%v(3)/ amass
		       vpar  = (vx*bx+vy*by+vz*bz)/max(sqrt(bx**2+by**2+bz**2),1.e-4)
		       vtot  = sqrt(vx**2+vy**2+vz**2)
		       
			   angle = acos( max(min(abs(vpar)/max(vtot,1.e-6),1.),-1.))   ![0., pi/2]
		       ee    = alog10(abs(0.5*vtot**2 * real(ut0)))
               if(ee<e_min .or. ee>e_max) goto 881

               imc1  = min(max(int(max(ee-e_min,0.)/delta_e+1.),1),n_e-1)
		       dmc1  = imc1-(ee - e_min)/delta_e  

               imc2  = min(max(int(max(angle,0.)/delta_alpha+1.),1),n_alpha-1)
		       dmc2  = imc2-(angle - 0.)/delta_alpha 

               imc3  = min(max(int(max(r-r_min,0.)/delta_r+1.),1),n_r-1)
		       dmc3  = imc3-(r - r_min)/delta_r  

               imc4  = min(max(int(max(phi0_,0.)/delta_phi+1.),1),n_phi-1)
		       dmc4  = imc4-(phi0_ - 0.)/delta_phi  

               ff(imc1,imc2,imc3,imc4) = ff(imc1,imc2,imc3,imc4) + dmc1*dmc2*dmc3*dmc4*Factor	          
               ff(imc1,imc2,imc3,imc4+1) = ff(imc1,imc2,imc3,imc4+1) + dmc1*dmc2*dmc3*(1.-dmc4)*Factor	          
               ff(imc1,imc2,imc3+1,imc4) = ff(imc1,imc2,imc3+1,imc4) + dmc1*dmc2*(1.-dmc3)*dmc4*Factor	          
               ff(imc1,imc2,imc3+1,imc4+1) = ff(imc1,imc2,imc3+1,imc4+1) + dmc1*dmc2*(1.-dmc3)*(1.-dmc4)*Factor	          

               ff(imc1,imc2+1,imc3,imc4) = ff(imc1,imc2+1,imc3,imc4) + dmc1*(1.-dmc2)*dmc3*dmc4*Factor	          
               ff(imc1,imc2+1,imc3,imc4+1) = ff(imc1,imc2+1,imc3,imc4+1) + dmc1*(1.-dmc2)*dmc3*(1.-dmc4)*Factor	          
               ff(imc1,imc2+1,imc3+1,imc4) = ff(imc1,imc2+1,imc3+1,imc4) + dmc1*(1.-dmc2)*(1.-dmc3)*dmc4*Factor	          
               ff(imc1,imc2+1,imc3+1,imc4+1) = ff(imc1,imc2+1,imc3+1,imc4+1) + dmc1*(1.-dmc2)*(1.-dmc3)*(1.-dmc4)*Factor	          

               ff(imc1+1,imc2,imc3,imc4) = ff(imc1+1,imc2,imc3,imc4) + (1.-dmc1)*dmc2*dmc3*dmc4*Factor	          
               ff(imc1+1,imc2,imc3,imc4+1) = ff(imc1+1,imc2,imc3,imc4+1) + (1.-dmc1)*dmc2*dmc3*(1.-dmc4)*Factor	          
               ff(imc1+1,imc2,imc3+1,imc4) = ff(imc1+1,imc2,imc3+1,imc4) + (1.-dmc1)*dmc2*(1.-dmc3)*dmc4*Factor	          
               ff(imc1+1,imc2,imc3+1,imc4+1) = ff(imc1+1,imc2,imc3+1,imc4+1) + dmc1*dmc2*(1.-dmc3)*(1.-dmc4)*Factor	          

               ff(imc1+1,imc2+1,imc3,imc4) = ff(imc1+1,imc2+1,imc3,imc4) + (1.-dmc1)*(1.-dmc2)*dmc3*dmc4*Factor	          
               ff(imc1+1,imc2+1,imc3,imc4+1) = ff(imc1+1,imc2+1,imc3,imc4+1) + (1.-dmc1)*(1.-dmc2)*dmc3*(1.-dmc4)*Factor	          
               ff(imc1+1,imc2+1,imc3+1,imc4) = ff(imc1+1,imc2+1,imc3+1,imc4) + (1.-dmc1)*(1.-dmc2)*(1.-dmc3)*dmc4*Factor	          
               ff(imc1+1,imc2+1,imc3+1,imc4+1) = ff(imc1+1,imc2+1,imc3+1,imc4+1) +  &
		                                     (1.-dmc1)*(1.-dmc2)*(1.-dmc3)*(1.-dmc4)*Factor	  
881            continue											         
            endif

         endif
     enddo
  enddo

  CALL MPI_ALLREDUCE(ff,fft,size(ff),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  fft(:,:,:,1)         = fft(:,:,:,1) + fft(:,:,:,n_phi)
  fft(:,:,:,n_phi)     = fft(:,:,:,1)
  fft(1,:,:,:)         = 2.*fft(1,:,:,:)
  fft(n_e,:,:,:)       = 2.*fft(n_e,:,:,:)
  fft(:,1,:,:)         = 2.*fft(:,1,:,:)
  fft(:,n_alpha,:,:)   = 2.*fft(:,n_alpha,:,:)
  ffp                  = fft
  do i=1,n_r
     do j=1,n_phi
	    fft(:,:,i,j) = fft(:,:,i,j)/vol_mc(i,j)/udensity
     enddo
  enddo


  if( mype == 0) then
      call check( nf90_put_var(ncid, varid(11),fft) )
      call check( nf90_put_var(ncid, varid(12),ffp) )
  endif
  deallocate (ff,fft,ffp)

  call allocate_bscalar(temp,1)
  call allocate_bscalar(temp1,1)
! write ni, Tpar Tper at equator
  allocate(plot3(n_r,n_phi,1))
! density
  do m=1,mblocks
	 temp(m)%scalar   = 0.
	 do n=1,kinds
	    temp(m)%scalar  = temp(m)%scalar + block(m)%ini(n,:,:,:)%ni * udensity !/ens
     enddo
  enddo
  call dconvert(temp,plot3,MCeq)
  if(mype==0) call check( nf90_put_var(ncid, varid(13),plot3(:,:,1)) )	  
! temperature
  do m=1,mblocks
	 temp(m)%scalar   = 0.
	 temp1(m)%scalar  = 0.
	 do n=1,kinds
	    temp(m)%scalar   = temp(m)%scalar  + block(m)%ini(n,:,:,:)%ni * block(m)%ini(n,:,:,:)%Ti(1)*uT0
	    temp1(m)%scalar  = temp1(m)%scalar + block(m)%ini(n,:,:,:)%ni
     enddo
     temp(m)%scalar      = temp(m)%scalar/max(temp1(m)%scalar,1.e-4_p2)
  enddo
  call dconvert(temp,plot3,MCeq)
  if(mype==0) call check( nf90_put_var(ncid, varid(14),plot3(:,:,1)) )	  
  do m=1,mblocks
	 temp(m)%scalar   = 0.
	 temp1(m)%scalar  = 0.
	 do n=1,kinds
	    temp(m)%scalar   = temp(m)%scalar  + block(m)%ini(n,:,:,:)%ni * block(m)%ini(n,:,:,:)%Ti(2)*uT0
	    temp1(m)%scalar  = temp1(m)%scalar + block(m)%ini(n,:,:,:)%ni
     enddo
     temp(m)%scalar      = temp(m)%scalar/max(temp1(m)%scalar,1.e-4_p2)
  enddo
  call dconvert(temp,plot3,MCeq)
  if(mype==0) call check( nf90_put_var(ncid, varid(15),plot3(:,:,1)) )	  

  deallocate(plot3)
  allocate(plot3(nrpointB,ntpointB,nzpointB))

! record B ----------------------------------------------------------
! magnetic field
  do n = 1,3
     do m=1,mblocks
	    temp(m)%scalar = bcar(m)%vector(n,:,:,:)*ubfield !block(m)%ini%localB0(n,:,:,:)*ubfield
     enddo
	 call dconvert(temp,plot3,MCB)
     if(mype==0) call check( nf90_put_var(ncid, varid(15+n),plot3) )		  
  enddo

  call allocate_bscalar(temp,-1)
  call allocate_bscalar(temp1,-1)
  deallocate(plot3)
! record phi on ionosphere

  if(mype==0) then
     call check( nf90_put_var(ncid, varid(19),phi_ion(1:mion,1:nion)*uphi ))	
     call check( nf90_put_var(ncid, varid(20),Ajpar(1:mion,1:nion)*uj0/scaling) )	

     call check( nf90_close(ncid) )
  endif

  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output4mc_nc


!=================================================================================================
subroutine diagnos_particle66_nc
  use netcdf
  use global_parameters
  use grid_data
  use vector_functions
  use diagnos
  implicit none
  include 'mpif.h'
  integer :: nions
  real,pointer,dimension(:,:) :: send,recv
  integer,pointer,dimension(:) :: num_recv,request
  type(bvector_type), dimension(:), pointer    :: Beff,bv1_,ev1_
  real(p2),pointer,dimension(:,:) :: qvp_writi,qvp_write,dummy
  integer :: num_diag,num_particle_diag,num_particle_diag_pe,num,num_skip,m,L,icount,ierror
  integer :: nregion,n,kind,loop,Li_tot,Le_tot,Lt,i,j,k
  real(p2) :: r,theta,zeta,vx,vy,vz,x,y,z,vper,amass,p,q,w
  real(p2) :: vb(3),vb1(3),ve1(3),vv(3),uu(3),vpar,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,angle,csi, &
              bw(3),vw(3),x2b(3,3),vbper1(3),vbper2(3),gamma,vpar_,vper_,vdb,vde
  integer  :: ii,jj,kk,ip,jp,kp,ier,nbuff,nbuffi,nbuffe,nvec

  integer,pointer,dimension(:)   :: L_start,L_ok
  real,pointer,dimension(:,:) :: plot
  real,pointer,dimension(:,:) :: bb 
  real     :: bs(3),bsa(3)
  real  :: param(180)
  character*6 fldname(74)
  character*4 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(4),varid(30),dimids(2)
  logical :: yes
  integer :: status(mpi_status_size)
  real(8) :: xyz(3)


  interface
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

  end interface
!------------------------
  nregion = int(pdomain(1)+0.001)
  if(n_pregion <1) return
  allocate(bb(n_pregion,3))

  call allocate_bvector(beff,3)
  call allocate_bvector(bv1_,3)
  call allocate_bvector(ev1_,3)
  do m=1,mblocks
    Beff(m)%vector = bv0(m)%vector + bv1(m)%vector
  enddo
  call bconcarB(beff,Bcar)
  call bconcarB(bv1,Bv1_)
  call bcovcarE(ev1,ev1_)
  call allocate_bvector(beff,-3)

    do m=1,n_pregion
       i = min(max(int(domain_pregion(m,1)+domain_pregion(m,2))/2,1),nx)
       j = min(max(int(domain_pregion(m,3)+domain_pregion(m,4))/2,1),nx)
       k = min(max(int(domain_pregion(m,5)+domain_pregion(m,6))/2,1),nx)
       bs = 0.
	   do n=1,mblocks
          if(i>=block(n)%i0 .and. i<=block(n)%i1 .and. &
		     j>=block(n)%j0 .and. j<=block(n)%j1 .and. &
			 k>=block(n)%k0 .and. k<=block(n)%k1) bs = bcar(n)%vector(:,i,j,k)
	   enddo

       CALL MPI_ALLREDUCE(bs,bsa,3,mpi_real,MPI_SUM,myComm1,IERROR)
       bb(m,:) = bsa
    enddo

  nvec    = 12

  if(mype==0) then 
     nbuffi = 10000
     nbuffe = 10000
  	 allocate(QVP_writi(nvec,nbuffi)) !nbuff
  	 allocate(QVP_write(nvec,nbuffe)) !nbuff
  endif


! do ions first
  yes     = .false.
  do i = 1,int(pdomain(3))
     if(pdomain(3+i) < 10) yes = .true.
  enddo


!  ---------- do ions -----------------------------------
  if(yes) then
     nbuff = max(maxval(block%mi),2000)*2
     allocate(send(nvec, nbuff),num_recv(numberpe),request(numberpe))

     num =0
     do m=1,mblocks
        do L = 1,block(m)%mi
           kind  = block(m)%ion(L)%kind
	       amass = mions(kind)*aion
           p     = block(m)%ion(L)%p(1)
           q     = block(m)%ion(L)%p(2)
           w     = block(m)%ion(L)%p(3)
   	       vx    = block(m)%ion(L)%v(1)/amass
	       vy    = block(m)%ion(L)%v(2)/amass
	       vz    = block(m)%ion(L)%v(3)/amass
           do n = 1,nregion
			  if(p>= domain_pregion(n,1) .and. p<= domain_pregion(n,2) .and. & 
			     q>= domain_pregion(n,3) .and. q<= domain_pregion(n,4) .and. & 
			     w>= domain_pregion(n,5) .and. w<= domain_pregion(n,6) )then
                 num = num +1

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1
                 vb    = bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
				 vv      = block(m)%ion(L)%v/amass
				 vpar    = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)
				 if(num >=  nbuff) then
                    allocate (dummy(nvec,num-1))
					dummy   = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(nvec,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
                 endif
				 send(1,num) = block(m)%ion(L)%kind
				 send(2,num) = p
				 send(3,num) = q
				 send(4,num) = w
				 send(5,num) = block(m)%ion(L)%v(1)/amass
				 send(6,num) = block(m)%ion(L)%v(2)/amass
				 send(7,num) = block(m)%ion(L)%v(3)/amass
				 send(8,num) = vpar
				 send(9,num) = block(m)%ion(L)%w *block(m)%ion(L)%ww
				 goto 110
		       endif  
			enddo  !end of n --region
110			continue
	    enddo ! end of L
     enddo ! end of mblock

     call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)
     if(mype > 0 .and. num > 0) then
	       call mpi_isend(send(:,1:num),nvec*num, mpi_real,0, &
		                           mype, MPI_COMM_WORLD, request(mype),ier)

     endif
     if(mype==0) then
	    Li_tot = 0
        do i=1,num
		   Li_tot = Li_tot +1
		   if(Li_tot >= nbuffi) then
              allocate(dummy(nvec,Li_tot-1))
			  dummy = QVP_writi(:,1:Li_tot-1)
			  deallocate(QVP_writi)
			  nbuffi = nbuffi*1.5
			  allocate(QVP_writi(nvec,nbuffi))
			  QVP_writi(:,1:Li_tot-1) = dummy
			  deallocate(dummy)
		   endif
           QVP_writi(:,Li_tot) = send(:,i)
        enddo
     endif

     do L=1,numberpe-1
	    if(num_recv(L+1) > 0) then
	       if(mype==0) then
		      allocate(recv(nvec,num_recv(L+1)))
              call mpi_recv(recv,nvec*num_recv(L+1), mpi_real,L, &
		                           L, MPI_COMM_WORLD, status, ier)
		  
		      do i=1,num_recv(L+1)
		         Li_tot = Li_tot +1
		         if(Li_tot >= nbuffi) then
                    allocate(dummy(nvec,Li_tot-1))
			        dummy = QVP_writi(:,1:Li_tot-1)
			        deallocate(QVP_writi)
			        nbuffi = nbuffi*1.5
			        allocate(QVP_writi(nvec,nbuffi))
			        QVP_writi(:,1:Li_tot-1) = dummy
			        deallocate(dummy)
		         endif
                 QVP_writi(:,Li_tot) = recv(:,i)
              enddo
			  deallocate(recv)
           elseif(mype==L)then
		      call MPI_Wait(request(L), status, ier)
		   endif
		endif
	 enddo 
	 deallocate(send,request,num_recv)
  endif

! now do electrons
! num_skip = max(1,me/NUM_PARTICLE_DIAG_PE)
! if(qsload == 1) num_skip = 1


  yes     = .false.
  do i = 1,int(pdomain(3))
     if(pdomain(3+i) > 10) yes = .true.
  enddo
! ------------------ do electrons -----------------------------------
  if(yes) then
     nbuff = max(maxval(block%me),2000)*2
     allocate(send(nvec, nbuff),num_recv(numberpe),request(numberpe))
     num =0
     do m=1,mblocks
        do L = 1,block(m)%me
           p     = block(m)%ele(L)%p(1)
           q     = block(m)%ele(L)%p(2)
           w     = block(m)%ele(L)%p(3)
   	       vx    = block(m)%ele(L)%v(1)/aelectron
	       vy    = block(m)%ele(L)%v(2)/aelectron
	       vz    = block(m)%ele(L)%v(3)/aelectron
           do n = 1,nregion
			  if(p>= domain_pregion(n,1) .and. p<= domain_pregion(n,2) .and. & 
			     q>= domain_pregion(n,3) .and. q<= domain_pregion(n,4) .and. & 
			     w>= domain_pregion(n,5) .and. w<= domain_pregion(n,6) )then
                 num = num +1

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1
                 vb    = bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  


                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  


                 ve1   = ev1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             ev1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 ev1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 ev1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 ev1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             ev1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 ev1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 ev1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  


                 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
				 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)

                 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2
				 uu       = block(m)%ele(L)%v(:)/aelectron
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma 
				
				 vpar     = vb(1)*uu(1) + vb(2)*uu(2) + vb(3)*uu(3)
				 vde      = a_dot_b(vv,ve1)
				 vdb      = a_dot_b(vv,vb1)

				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)

                 angle    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 angle    = mod(angle+pi2,pi2)

                 
				 if(num >=  nbuff) then
                    allocate (dummy(nvec,num-1))
					dummy   = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(nvec,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
                 endif
				 send(1,num) = block(m)%ele(L)%kind
				 send(2,num) = p
				 send(3,num) = q
				 send(4,num) = w
				 send(5,num) = block(m)%ele(L)%v(1)/aelectron
				 send(6,num) = block(m)%ele(L)%v(2)/aelectron
				 send(7,num) = block(m)%ele(L)%v(3)/aelectron
				 send(8,num) = vpar
				 send(9,num) = block(m)%ele(L)%w *block(m)%ele(L)%ww
				 send(10,num)= csi
				 send(11,num)= vde
				 send(12,num)= vdb
				 goto 111
		       endif  
			enddo  !end of n --region
111			continue
	    enddo ! end of L
     enddo ! end of mblock


     call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)
     if(mype > 0 .and. num > 0) then
	       call mpi_isend(send(:,1:num),nvec*num, mpi_real,0, &
		                           mype, MPI_COMM_WORLD, request(mype),ier)

     endif

     if(mype==0) then
	    Le_tot = 0
        do i=1,num
		   Le_tot = Le_tot +1
		   if(Le_tot >= nbuffe) then
              allocate(dummy(nvec,Le_tot-1))
			  dummy = QVP_write(:,1:Le_tot-1)
			  deallocate(QVP_write)
			  nbuffe = nbuffe*1.5
			  allocate(QVP_write(nvec,nbuffe))
			  QVP_write(:,1:Le_tot-1) = dummy
			  deallocate(dummy)
		   endif
           QVP_write(:,Le_tot) = send(:,i)
        enddo
     endif

     do L=1,numberpe-1
	    if(num_recv(L+1) > 0) then
	       if(mype==0) then
		      allocate(recv(nvec,num_recv(L+1)))
              call mpi_recv(recv,nvec*num_recv(L+1), mpi_real,L, &
		                           L, MPI_COMM_WORLD, status, ier)
		  
		      do i=1,num_recv(L+1)
		         Le_tot = Le_tot +1
		         if(Le_tot >= nbuffe) then
                    allocate(dummy(nvec,Le_tot-1))
			        dummy = QVP_write(:,1:Le_tot-1)
			        deallocate(QVP_write)
			        nbuffe = nbuffe*1.5
			        allocate(QVP_write(nvec,nbuffe))
			        QVP_write(:,1:Le_tot-1) = dummy
			        deallocate(dummy)
		         endif
                 QVP_write(:,Le_tot) = recv(:,i)
              enddo
			  deallocate(recv)
           elseif(mype==L)then
		      call MPI_Wait(request(L), status, ier)
		   endif
		endif
	 enddo 
	 deallocate(send,request,num_recv)
  endif

  call allocate_bvector(bv1_,-3)
  call allocate_bvector(ev1_,-3)

! ---------------------------------------------------------------
  param      = 0.
  param(1)   = stime   ! has been scaled to omega_i
  param(2)   = Li_tot  ! has been scaled to omega_i
  param(3)   = Le_tot  ! has been scaled to omega_i

  param(4)   = nvec
  param(5)   = kinds
  param(6)   = betae
  param(7)   = nx
  param(8)   = ny
  param(9)   = nz

  param(10:kinds+9)   = fraci(1:kinds)
  param(20:kinds+19)  = frace(1:kinds)
  param(30)   = nregion
  param(31:30+nregion)   = domain_pregion(1:nregion,1)
  param(41:40+nregion)   = domain_pregion(1:nregion,2)
  param(51:50+nregion)   = domain_pregion(1:nregion,3)
  param(61:60+nregion)   = domain_pregion(1:nregion,4)
  param(71:70+nregion)   = domain_pregion(1:nregion,5)
  param(81:80+nregion)   = domain_pregion(1:nregion,6)
  param(91:90+nregion)   = bb(1:nregion,1)
  param(101:100+nregion) = bb(1:nregion,2)
  param(111:110+nregion) = bb(1:nregion,3)


  param(121)  = ulength*deltax
  param(122)  = ulength*deltay
  param(123)  = ulength*deltaz
  param(124)  = uspeed
  param(125)  = Va
  param(126)  = cspeed
  param(127)  = uT0

! ------------------------------
  idiagP = idiagP + 1
  write(ct,'(I4.4)')IDIAGP 
  if(Li_tot+Le_tot == 0 ) goto 831


  if(mype==0) then
        filename = 'particle'//ct//'.nc'
        call check( nf90_create(FILENAME, nf90_clobber, ncid) )

        call check( nf90_def_dim(ncid, 'param', 180, dimid(1)) )
        call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(1), varid(1)) )
        call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

        call check( nf90_def_dim(ncid, 'var', nvec, dimid(2)) )
        call check( nf90_def_var(ncid, 'var', NF90_REAL, dimid(2), varid(2)) )
        call check( nf90_put_att(ncid,  varid(2), UNITS, 'non') )

        call check( nf90_def_dim(ncid, 'nion', Li_tot, dimid(3)) )
        call check( nf90_def_var(ncid, 'nion', NF90_REAL, dimid(3), varid(3)) )
        call check( nf90_put_att(ncid,  varid(3), UNITS, 'non') )

        call check( nf90_def_dim(ncid, 'nele', Le_tot, dimid(4)) )
        call check( nf90_def_var(ncid, 'nele', NF90_REAL, dimid(4), varid(4)) )
        call check( nf90_put_att(ncid,  varid(4), UNITS, 'non') )


        dimids=(/dimid(2),dimid(3)/)
        call check( nf90_def_var(ncid, 'ion', NF90_REAL, dimids, varid(5)) )
        call check( nf90_put_att(ncid, varid(5), UNITS, 'non') )
        dimids=(/dimid(2),dimid(4)/)
        call check( nf90_def_var(ncid, 'ele', NF90_REAL, dimids, varid(6)) )
        call check( nf90_put_att(ncid, varid(6), UNITS, 'non') )


        call check( nf90_enddef(ncid) )

        call check( nf90_put_var(ncid, varid(1), param) )
		allocate(plot(nvec,Li_tot))
		plot = QVP_writi(:,1:Li_tot)
        call check( nf90_put_var(ncid, varid(5), plot)  )

        deallocate(plot)
		allocate(plot(nvec,Le_tot))
		plot = QVP_write(:,1:Le_tot)

        call check( nf90_put_var(ncid, varid(6), plot) )
		deallocate(plot)

        call check( nf90_close(ncid) )
  endif
831  continue
   
  if(mype ==0 ) deallocate(qvp_writi,qvp_write)
  deallocate(bb)

  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine diagnos_particle66_nc





! -------------------------------------------------------------------
! write E& B (x,t) data for test particle simulations 
subroutine write_field_nc1file()
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  implicit none
  include 'mpif.h'
  real(p2),  pointer,dimension(:,:,:,:) :: plot3
  real,  pointer,dimension(:,:,:,:) :: plot1
  integer :: ic,i,iplot,m,n,iflag
  character*1 fldname(2)
  character*2 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(5),varid(30),dimids(5),icount(5),istart(5),icount1(1),istart1(1)

  save ic
  data ic/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine gather_vector(c,d)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
        real(p2),pointer,dimension(:,:,:,:) :: d
      end subroutine gather_vector
  end interface

  fldname = (/"B","E"/)

  if(abs(mod(simtype(1),10)) <=5 ) then
     write(ct,'(I2.2)')(istep-1)/20000
     ic      = mod((istep-1),20000) + 1   !max(mod(ic,20000),1) 
  else
     write(ct,'(I2.2)')(istep-1)/500
     ic      = mod((istep-1),500) + 1   !max(mod(ic,20000),1) 
  endif

  icount1 = (/1/)
  istart1 = (/ic/)


  istart = (/1, 1, 1, 1, ic /)
  if(mype==0 ) then
        icount = (/3, nxp, nyp, nzp, 1 /)
		filename = 'EandB'//ct//'.nc'

!       Create the file. 
        if(ic == 1) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'X', nxp, dimid(2)) )
          call check( nf90_def_var(ncid, 'X', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'Y', nyp, dimid(3)) )
          call check( nf90_def_var(ncid, 'Y', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'Z', nzp, dimid(4)) )
          call check( nf90_def_var(ncid, 'Z', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          dimids=(/dimid(1),dimid(2),dimid(3),dimid(4),dimid(5)/)
          do iplot = 1, 2
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(5+iplot)) )
           call check( nf90_put_att(ncid, varid(5+iplot), UNITS, 'non') )
          enddo

          call check( nf90_enddef(ncid) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 2
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(5+iplot)) )
          enddo
        endif
  endif


  allocate(plot3(3,nxp,nyp,nzp),plot1(3,nxp,nyp,nzp))


  do iplot = 1, 2
     if(iplot==1) call gather_vector(bv1,plot3)
     if(iplot==2) call gather_vector(ev1,plot3)

	 plot1  = plot3

     if(mype==0) call check( nf90_put_var(ncid, varid(5+iplot), plot1, start = istart, &
                              count = icount) )

  enddo

  if(mype==0) call check( nf90_close(ncid) )

  deallocate(plot3,plot1)

  ic  = ic +1
  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine write_field_nc1file
! -------------------------------------------------------------------------------------------

! -------------------------------------------------------------------
! write E& B (x,t) data for test particle simulations 
subroutine read_field_nc1file()
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  implicit none
  include 'mpif.h'
  real(p2),  dimension(:,:,:,:,:), allocatable :: Bsave
  real,  pointer,dimension(:,:,:,:,:) :: EBsave0,dummy
  real,  pointer,dimension(:,:,:,:,:) :: send,recv
  integer :: ic,i,iplot,m,n,iflag
  character*1 fldname(2)
  character*2 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: dimid(5),varid(30),dimids(5),icount(5),istart(5),icount1(1),istart1(1),icount2
  logical :: transferx(3)

  character(len=80) :: infile, dname, attname, attvalue, varname
  character(len=80), dimension(:), allocatable :: dimname
  integer :: ncid, ndims, nvars, nglobalatts, unlimdimid,ierror,ier,num,k,itime
  integer :: dimlen, attlen, att_type
  integer :: vardims, varnatts, vartype
  integer, dimension(nf90_max_var_dims) :: vardimids
  integer, dimension(:), allocatable :: dimsize
  integer :: status(mpi_status_size)

  save dimsize,transferx
  save ic
  data ic/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
  end interface

!  transfer = .false.
! transfer = .true.

  if(abs(mod(simtype(1),10)) <=5 ) then
     write(ct,'(I2.2)')(istep-1)/20000
     ic      = mod((istep-1),20000) + 1   !max(mod(ic,20000),1) 
  else
     write(ct,'(I2.2)')(istep-1)/500
     ic      = mod((istep-1),500) + 1   !max(mod(ic,20000),1) 
  endif

  if(ic > 1 .and. ic > dimsize(5)) call end_run

  if(ic == 1 ) then
     if(mype == 0 ) then
		filename = '../tps/EandB'//ct//'.nc'
        call check( nf90_open(filename, nf90_nowrite, ncid) )
        call check( nf90_inquire(ncid, nDims, nVars, nGlobalAtts, unlimdimid) )

 	    if(allocated (dimsize)) deallocate(dimsize)
        allocate(dimsize(ndims))

        do n = 1,ndims
           call check( nf90_inquire_dimension(ncid, n, dname, dimlen) )
           write(*,*) '>>> dimension name = ',trim(dname), ' length = ',dimlen
           dimsize(n) = dimlen
        enddo

        do n = 1,nvars
           call check( nf90_inquire_variable(ncid, n, varname, vartype,  &
                    vardims, vardimids, varnatts) )
           !write(*,*,advance='no') '>>> variable name = ',varname, '('
            write(*,*) '>>> variable name = ',trim(varname)
        enddo

	    allocate(Bsave(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
		allocate(EBsave0(6,dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
        call check( nf90_get_var(ncid, 6, Bsave) )
		EBsave0(1:3,:,:,:,:) = Bsave
        call check( nf90_get_var(ncid, 7, Bsave) )
		EBsave0(4:6,:,:,:,:) = Bsave

		deallocate(Bsave)

        call check( nf90_close(ncid) )
     endif   ! end of mype ==0

! --------- need 2-D test particle but data here is in 1D -----------------------
     if(mype == 0) then
	    transferx(1) = dimsize(2) /= nxp .and. dimsize(2)==1
	    transferx(2) = dimsize(3) /= nyp .and. dimsize(3)==1
	    transferx(3) = dimsize(4) /= nzp .and. dimsize(4)==1
	 endif
! --------------------------------------------------------------------------------
     call MPI_Barrier( mpi_comm_world, ierror)
     call mpi_bcast(ndims,1,MPI_integer,0,mpi_comm_world,ierror)
     call mpi_bcast(transferx,3,MPI_logical,0,mpi_comm_world,ierror)
	 if(mype /=0 )then
 	    if(allocated (dimsize)) deallocate(dimsize)
	    allocate(dimsize(ndims))
     endif
     call mpi_bcast(dimsize,ndims,MPI_integer,0,mpi_comm_world,ierror)
	 do m=1,mblocks
	    if(allocated(block(m)%EBsave)) deallocate(block(m)%EBsave)
		if(transferx(1)) then
		   allocate(block(m)%EBsave(6,1, &
		                           block(m)%ny0:block(m)%ny1, &
						           block(m)%nz0:block(m)%nz1,dimsize(5)))

		else
		   allocate(block(m)%EBsave(6,block(m)%nx0:block(m)%nx1, &
		                           block(m)%ny0:block(m)%ny1, &
						           block(m)%nz0:block(m)%nz1,dimsize(5)))
        endif
	 enddo


     call MPI_Barrier( mpi_comm_world, ierror)

     if(mype==0) then
        do m=1,mblocks
		   if(transferx(1)) then
              block(m)%EBsave(:,:, &
		                     max(block(m)%ny0,1):block(m)%ny1, &
						     max(block(m)%nz0,1):block(m)%nz1,:) = &
				       EBsave0(:,:, &
		                     max(block(m)%ny0,1):block(m)%ny1, &
						     max(block(m)%nz0,1):block(m)%nz1,:)

		   else
              block(m)%EBsave(:,max(block(m)%nx0,1):block(m)%nx1, &
		                     max(block(m)%ny0,1):block(m)%ny1, &
						     max(block(m)%nz0,1):block(m)%nz1,:) = &
				       EBsave0(:,max(block(m)%nx0,1):block(m)%nx1, &
		                     max(block(m)%ny0,1):block(m)%ny1, &
						     max(block(m)%nz0,1):block(m)%nz1,:)
           endif
		enddo
     endif

     call MPI_Barrier( mpi_comm_world, ierror)

     do itime = 1, dimsize(5),1000

     if(mype==0) then
! sending informations --------------
		do k=2,numberpe
           num = Eachblocks(k)
		   do m=1,num
		      if(transferx(1)) then
                 allocate(send(6,1, &
			                  Eachgrid(k,m,3):Eachgrid(k,m,4), &
			                  Eachgrid(k,m,5):Eachgrid(k,m,6), &
			                  1000))
			  
			     send = EBsave0(:,:, &
			                   Eachgrid(k,m,3):Eachgrid(k,m,4), &
			                   Eachgrid(k,m,5):Eachgrid(k,m,6), &
			                   itime:itime+999)

			  else
                 allocate(send(6,Eachgrid(k,m,1):Eachgrid(k,m,2), &
			                  Eachgrid(k,m,3):Eachgrid(k,m,4), &
			                  Eachgrid(k,m,5):Eachgrid(k,m,6), &
			                  1000))
			  
			     send = EBsave0(:,Eachgrid(k,m,1):Eachgrid(k,m,2), &
			                   Eachgrid(k,m,3):Eachgrid(k,m,4), &
			                   Eachgrid(k,m,5):Eachgrid(k,m,6), &
			                   itime:itime+999)
              endif
			  call mpi_send(send, size(send), mpi_real, k-1, m, mpi_comm_world, ier)
              deallocate(send)
		   enddo
		enddo
	 endif

     if(mype/=0) then
	   do m=1,mblocks
	      if(transferx(1)) then
             allocate (recv(6,1, &
		                   max(block(m)%ny0,1):block(m)%ny1, &
						   max(block(m)%nz0,1):block(m)%nz1,1000))
		  
             call mpi_recv(recv,size(recv),mpi_real,0,m,mpi_comm_world,status,ier)
             block(m)%EBsave(:,:, &
		                    max(block(m)%ny0,1):block(m)%ny1, &
						    max(block(m)%nz0,1):block(m)%nz1,itime:itime+9)  = recv

		  else
             allocate (recv(6,max(block(m)%nx0,1):block(m)%nx1, &
		                   max(block(m)%ny0,1):block(m)%ny1, &
						   max(block(m)%nz0,1):block(m)%nz1,1000))
		  
             call mpi_recv(recv,size(recv),mpi_real,0,m,mpi_comm_world,status,ier)
             block(m)%EBsave(:,max(block(m)%nx0,1):block(m)%nx1, &
		                    max(block(m)%ny0,1):block(m)%ny1, &
						    max(block(m)%nz0,1):block(m)%nz1,itime:itime+9)  = recv
          endif
          deallocate(recv)
	   enddo
	 endif

     call MPI_Barrier( mpi_comm_world, ierror)

     enddo  !end of itime -------
	 if(mype==0) deallocate(EBsave0)
     if(mype==0) write(*,*) '>>> readding done and send & recv --------'
  endif


  do m=1,mblocks
     if(transferx(1)) then
	    do i=1,nxp
           bdriver(m)%vector(:,i,:,:) = block(m)%EBsave(1:3,1,:,:,ic)
           edriver(m)%vector(:,i,:,:) = block(m)%EBsave(4:6,1,:,:,ic)
        enddo
	 else
        bdriver(m)%vector = block(m)%EBsave(1:3,:,:,:,ic)
        edriver(m)%vector = block(m)%EBsave(4:6,:,:,:,ic)
	 endif
  enddo


  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine read_field_nc1file
! -------------------------------------------------------------------------------------------
subroutine output_fetps_nc1file(kind)  ! test particles
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  implicit none
  include 'mpif.h'
  type(TestParticle_type) :: TP
  real,  pointer,dimension(:) :: array_e,array_a
  real,  pointer,dimension(:,:) :: flux,fluxt
  integer :: ic,ic0(4),i,iplot,m,n,iflag
  character*3 fldname(2)
  character*3 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  :: param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,angle,vtot
  real(p2)  :: amass,wx1,wy1,wz1,wx0,wy0,wz0,bc(3)
  real(p2)  :: uu(3),gamma,beq0,B !,vvmin,vvmax,dv
  integer   :: kind,ii,jj,kk,ip,jp,kp
  integer :: j,k,L,imc1,imc2,index,IERROR
  logical    file_exist
  real*8    :: xyz(3),p,q,w

  save ic,ic0
  data ic/0/,ic0/0,0,0,0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
 
   
!  vvmin      = -3000/uspeed
!  vvmax      =  3000/uspeed

!  dv         =  (vvmax-vvmin)/real(nv-1)
  allocate(array_e(nv),array_a(na))
  do m=1,nv
     array_e(m) = emin + de*(m-1.)
  enddo
  do m=1,na
     array_a(m) = amin + dalpha*(m-1.)
  enddo

  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed




  param      = 0.
  param(1)   = nV
  param(2)   = na
  param(3)   = Emin
  param(4)   = Emax
  param(5)   = amin
  param(6)   = amax
  param(7)   = vvmin(kind)
  param(8)   = vvmax(kind)
  param(9)   = dE
  param(10)  = dalpha
  param(11)  = dv(kind)
  param(21)  = utime
  param(22)  = uspeed
  param(23)  = ut0
  param(24)  = cspeed

  fldname = (/"PSD"," E "/)

  write(ct,'(I3.3)')((istep-istep0)/ndiagfe - 1)/100000
  ic      = mod(((istep-istep0)/ndiagfe-1),100000) + 1   !max(mod(ic,20000),1) 
  ic0(kind)     = ic0(kind)+1

  if(kind==1) filename = 'feP'//ct//'.nc'
  if(kind==2) filename = 'feHe'//ct//'.nc'
  if(kind==3) filename = 'feO'//ct//'.nc'
  if(kind==4) filename = 'feE'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if(.not. file_exist  .or. ic ==0  .or. (ic0(kind)==1 .and. irun ==0)) then
      ic0_fe(kind) = ic - 1   
  endif
  ic = ic - ic0_fe(kind)

  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, ic /)
  if(mype==0 ) then
        icount = (/nv, na, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )


          call check( nf90_def_dim(ncid, 'nE', nv, dimid(3)) )
          call check( nf90_def_var(ncid, 'nE', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'eV') )

          call check( nf90_def_dim(ncid, 'nalpha', na, dimid(4)) )
          call check( nf90_def_var(ncid, 'nalpha', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, ' ') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )

          dimids=(/dimid(3),dimid(4),dimid(5)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(3), array_e) )
          call check( nf90_put_var(ncid, varid(4), array_a) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif


  allocate(flux(nv,na),fluxt(nv,na))
  flux = 0.

  do m=1,mblocks
	 if(kind ==1) TP = block(m)%TPi
	 if(kind ==2) TP = block(m)%TPHe
     if(kind ==3) TP = block(m)%TPO
     if(kind ==4) TP = block(m)%TPe
     dO  L = 1, TP%mi
         Factor= TP%qv(9,L)
         amass  = TP%amass

         p        = TP%qv(1,L)
         q        = TP%qv(2,L)
         w        = TP%qv(3,L)

	     wx1      = p +1.0 
         ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	     ip       = min(block(m)%nx1,ii+1)
	     wx1      = wx1 - ii
	     wx0      = 1.0 - wx1

	     wy1      = q +1.0 
         jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	     jp       = min(block(m)%ny1,jj+1)
	     wy1      = wy1 - jj
	     wy0      = 1.0 - wy1

         wz1      = w +1.0 
         kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	     kp       = min(block(m)%nz1,kk+1)
	     wz1      = wz1 - kk
	     wz0      = 1.0 - wz1

	 	 bx       = bcar(m)%vector(1,ii,jj,kk)
		 by       = bcar(m)%vector(2,ii,jj,kk)
		 bz       = bcar(m)%vector(3,ii,jj,kk)
		 B        = sqrt(bx**2+by**2+bz**2)

         uu(1:3)  = TP%qv(4:6,L)/amass
	     gamma    = sqrt(1.+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	     vpar     = (uu(1)*bx+uu(2)*by+uu(3)*bz)/max(sqrt(bx**2+by**2+bz**2),1.e-4)
		 vtot     = sqrt(uu(1)**2+uu(2)**2+uu(3)**2)
		 vper     = sqrt(abs(vtot**2 - vpar**2))
		 angle    = acos( max(min(abs(vpar)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2]
! to get equatorial pitch angle for dipole case

         if(dipole) then
		    xyz           = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
            x             = xyz(1)
            y             = xyz(2)
            z             = xyz(3)
			bc            = dipole3DB(x,y,z)
		    B             = sqrt(bc(1)**2 + bc(2)**2 + bc(3)**2)

		    xyz           = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !get magnetic field at equator
            x             = xyz(1)
            y             = xyz(2)
            z             = xyz(3)
			bc            = dipole3DB(x,y,z)
		    Beq0          = sqrt(bc(1)**2 + bc(2)**2 + bc(3)**2)

			vper          = vper * sqrt(Beq0/max(B,1.e-6))
		    angle         = asin( max(min(abs(vper)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2]
		 endif

	     if(gamma <1.e-6) then
	        ee       = 0.5*amass * (uu(1)**2+uu(2)**2+uu(3)**2) *uT0
	     else
	        ee       = (gamma-1.)*cspeed**2*amass * uT0
	     endif
		 ee        = log10(max(ee,1.e-10_p2))

         if(ee<emin .or. ee>emax) goto 881
            imc1  = min(max(int(max(ee-emin,0.)/de+1.),1),nv-1)
		    dmc1  = imc1-(ee - emin)/de  

            imc2  = min(max(int(max(angle-amin,0.)/dalpha+1.),1),na-1)
		    dmc2  = imc2-(angle - amin)/dalpha
            flux(imc1,imc2)     = flux(imc1,imc2) + dmc1*dmc2*Factor	          
            flux(imc1,imc2+1)   = flux(imc1,imc2+1) + dmc1*(1.-dmc2)*Factor	          
            flux(imc1+1,imc2)   = flux(imc1+1,imc2) + (1.-dmc1)*dmc2*Factor	          
            flux(imc1+1,imc2+1) = flux(imc1+1,imc2+1) + (1.-dmc1)*(1.-dmc2)*Factor	          
881       continue
	 enddo
  enddo

  flux  = flux + flux_loss(kind,:,:)
  CALL MPI_ALLREDUCE(flux,fluxt,size(flux),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  fluxt(1,:)        = 2.*fluxt(1,:)
  fluxt(nv,:)       = 2.*fluxt(nv,:)
  fluxt(:,1)        = 2.*fluxt(:,1)
  fluxt(:,na)       = 2.*fluxt(:,na)

  do iplot = 1, 1

     if(mype==0) call check( nf90_put_var(ncid, varid(7+iplot), fluxt, start = istart, &
                              count = icount) )

  enddo

  if(mype==0) call check( nf90_close(ncid) )

  deallocate(flux,fluxt,array_e,array_a)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_fetps_nc1file

! -------------------------------------------------------------------------------------------

subroutine output_fvtps_nc1file(kind)  ! test particles  !add f(vpar,vper) 
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  implicit none
  include 'mpif.h'
  type(TestParticle_type) :: TP
  real,  pointer,dimension(:) :: array_vx,array_vy,array_vper
  real,  pointer,dimension(:,:,:) :: fv,fvt
  real,  pointer,dimension(:,:) :: plot

  integer :: ic,ic0(4),i,iplot,m,n,iflag
  character*4 fldname(4)
  character*3 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(15),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  :: param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,dmc3,dmc4,dmc5,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,angle,vtot,dvper
  real(p2)  :: amass,wx1,wy1,wz1,wx0,wy0,wz0,vv(3)
  real(p2)  :: uu(3),gamma
  integer   :: kind,ii,jj,kk,ip,jp,kp
  integer :: j,k,L,imc1,imc2,imc3,imc4,imc5,index,IERROR
  logical file_exist
  real*8  :: xyz(3),p,q,w

  save ic,ic0
  data ic/0/,ic0/0,0,0,0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
 
 
!  vvmin      = -3000./uspeed
!  vvmax      =  3000./uspeed
!  if(kind == 4) then
!     vvmin      = -10000./uspeed
!     vvmax      =  10000./uspeed
!  endif

  dv(kind)         =  (vvmax(kind)-vvmin(kind))/real(nv-1)
  dvper            =  dv(kind)/2.
  allocate(array_vx(nv),array_vy(nv),array_vper(nv))
  do m=1,nv
     array_vx(m) = vvmin(kind) + dv(kind)*(m-1.)
     array_vy(m) = vvmin(kind) + dv(kind)*(m-1.)
  enddo

  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed



  param      = 0.
  param(1)   = nV
  param(2)   = nV
  param(3)   = 0
  param(4)   = 0
  param(5)   = 0
  param(6)   = 0
  param(7)   = vvmin(kind)
  param(8)   = vvmax(kind)
  param(9)   = 0
  param(10)  = 0
  param(11)  = dv(kind)
  param(12)  = dvper
  param(21)  = utime
  param(22)  = uspeed
  param(23)  = ut0
  param(24)  = cspeed

  fldname = (/"fvxy","fvxz","fvyz","fvpp"/)

  write(ct,'(I3.3)')(istep/ndiagfv - 1)/1000
  ic            = mod((istep/ndiagfv-1),1000) + 1   !max(mod(ic,20000),1)
  ic0(kind)     = ic0(kind) + 1 

  if(kind==1) filename = 'fvP'//ct//'.nc'
  if(kind==2) filename = 'fvHe'//ct//'.nc'
  if(kind==3) filename = 'fvO'//ct//'.nc'
  if(kind==4) filename = 'fvE'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if(.not. file_exist  .or. ic ==0  .or. (ic0(kind)==1 .and. irun ==0)) then
      ic0_fv(kind) = ic - 1   
  endif
  ic = ic - ic0_fv(kind)


  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, ic /)
  if(mype==0 ) then
        icount = (/nv, nv, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )


          call check( nf90_def_dim(ncid, 'nvx', nv, dimid(3)) )
          call check( nf90_def_var(ncid, 'nvx', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'km/s') )

          call check( nf90_def_dim(ncid, 'nvy', nv, dimid(4)) )
          call check( nf90_def_var(ncid, 'nvy', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'km/s') )

          call check( nf90_def_dim(ncid, 'nvp', nv, dimid(7)) )
          call check( nf90_def_var(ncid, 'nvp', NF90_REAL, dimid(7), varid(7)) )
          call check( nf90_put_att(ncid,  varid(7), UNITS, 'km/s') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 's') )

          dimids=(/dimid(3),dimid(4),dimid(5)/)
          do iplot = 1, 4
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(8+iplot)) )
           call check( nf90_put_att(ncid, varid(8+iplot), UNITS, 'psd') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(3), array_vx) )
          call check( nf90_put_var(ncid, varid(4), array_vy) )
          call check( nf90_put_var(ncid, varid(7), array_vper) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 4
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(8+iplot)) )
          enddo
        endif
  endif


  allocate(fv(4,nv,nv),fvt(4,nv,nv),plot(nv,nv))
  fv = 0.
  do m=1,mblocks
	 if(kind ==1) TP = block(m)%TPi
	 if(kind ==2) TP = block(m)%TPHe
     if(kind ==3) TP = block(m)%TPO
     if(kind ==4) TP = block(m)%TPe
     dO  L = 1, TP%mi
         Factor = TP%qv(9,L)
         amass  = TP%amass

         p        = TP%qv(1,L)
         q        = TP%qv(2,L)
         w        = TP%qv(3,L)

	     wx1      = p +1.0 
         ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	     ip       = min(block(m)%nx1,ii+1)
	     wx1      = wx1 - ii
	     wx0      = 1.0 - wx1

	     wy1      = q +1.0 
         jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	     jp       = min(block(m)%ny1,jj+1)
	     wy1      = wy1 - jj
	     wy0      = 1.0 - wy1

         wz1      = w +1.0 
         kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	     kp       = min(block(m)%nz1,kk+1)
	     wz1      = wz1 - kk
	     wz0      = 1.0 - wz1

	 	 bx       = bcar(m)%vector(1,ii,jj,kk)
		 by       = bcar(m)%vector(2,ii,jj,kk)
		 bz       = bcar(m)%vector(3,ii,jj,kk)


         uu(1:3)  = TP%qv(4:6,L)/amass
	     gamma    = sqrt(1.+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
		 vv       = uu/gamma

	     vpar     = (vv(1)*bx+vv(2)*by+vv(3)*bz)/max(sqrt(bx**2+by**2+bz**2),1.e-4)
		 vtot     = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
		 angle    = acos( max(min(abs(vpar)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2]
		 vper     = sqrt(abs(vtot**2-vpar**2))


         if(vv(1)<vvmin(kind) .or. vv(1)>=vvmax(kind) .or. &
            vv(2)<vvmin(kind) .or. vv(2)>=vvmax(kind) .or. &
            vv(3)<vvmin(kind) .or. vv(3)>=vvmax(kind) ) goto 881
            
			imc1  = min(max(int(max(vv(1)-vvmin(kind),0.)/dv(kind)+1.),1),nv-1)
		    dmc1  = imc1-(vv(1) - vvmin(kind))/dv(kind)
            imc2  = min(max(int(max(vv(2)-vvmin(kind),0.)/dv(kind)+1.),1),nv-1)
		    dmc2  = imc2-(vv(2) - vvmin(kind))/dv(kind)
            imc3  = min(max(int(max(vv(3)-vvmin(kind),0.)/dv(kind)+1.),1),nv-1)
		    dmc3  = imc3-(vv(3) - vvmin(kind))/dv(kind)

            imc4  = min(max(int(max(vpar-vvmin(kind),0.)/dv(kind)+1.),1),nv-1)
		    dmc4  = imc4-(vpar - vvmin(kind))/dv(kind)

            imc5  = min(max(int(max(vper,0.)/dvper+1.),1),nv-1)
		    dmc5  = imc5-(vper)/dvper


            fv(1,imc1,imc2)     = fv(1,imc1,imc2) + dmc1*dmc2*Factor	          
            fv(1,imc1,imc2+1)   = fv(1,imc1,imc2+1) + dmc1*(1.-dmc2)*Factor	          
            fv(1,imc1+1,imc2)   = fv(1,imc1+1,imc2) + (1.-dmc1)*dmc2*Factor	          
            fv(1,imc1+1,imc2+1) = fv(1,imc1+1,imc2+1) + (1.-dmc1)*(1.-dmc2)*Factor	          

            fv(3,imc2,imc3)     = fv(3,imc2,imc3) + dmc2*dmc3*Factor	          
            fv(3,imc2,imc3+1)   = fv(3,imc2,imc3+1) + dmc2*(1.-dmc3)*Factor	          
            fv(3,imc2+1,imc3)   = fv(3,imc2+1,imc3) + (1.-dmc2)*dmc3*Factor	          
            fv(3,imc2+1,imc3+1) = fv(3,imc2+1,imc3+1) + (1.-dmc2)*(1.-dmc3)*Factor	          

            fv(2,imc1,imc3)     = fv(2,imc1,imc3) + dmc1*dmc3*Factor	          
            fv(2,imc1,imc3+1)   = fv(2,imc1,imc3+1) + dmc1*(1.-dmc3)*Factor	          
            fv(2,imc1+1,imc3)   = fv(2,imc1+1,imc3) + (1.-dmc1)*dmc3*Factor	          
            fv(2,imc1+1,imc3+1) = fv(2,imc1+1,imc3+1) + (1.-dmc1)*(1.-dmc3)*Factor	          

            fv(4,imc4,imc5)     = fv(4,imc4,imc5) + dmc4*dmc5*Factor	          
            fv(4,imc4,imc5+1)   = fv(4,imc4,imc5+1) + dmc4*(1.-dmc5)*Factor	          
            fv(4,imc4+1,imc5)   = fv(4,imc4+1,imc5) + (1.-dmc4)*dmc5*Factor	          
            fv(4,imc4+1,imc5+1) = fv(4,imc4+1,imc5+1) + (1.-dmc4)*(1.-dmc5)*Factor	          

        
881       continue

	 enddo
  enddo

  CALL MPI_ALLREDUCE(fv,fvt,size(fv),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  fvt(:,1,:)        = 2.*fvt(:,1,:)
  fvt(:,nv,:)       = 2.*fvt(:,nv,:)
  fvt(:,:,1)        = 2.*fvt(:,:,1)
  fvt(:,:,nv)       = 2.*fvt(:,:,nv)

  do iplot = 1, 4
     plot = fvt(iplot,:,:)
     if(mype==0) call check( nf90_put_var(ncid, varid(8+iplot), plot, start = istart, &
                              count = icount) )

  enddo

  if(mype==0) call check( nf90_close(ncid) )

  deallocate(fv,fvt,array_vx,array_vy,array_vper)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_fvtps_nc1file

!-----------------------------------------------------------------
subroutine output_DDtps_nc1file(kind)  ! test particles
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  implicit none
  include 'mpif.h'
  type(TestParticle_type) :: TP
  real,  pointer,dimension(:,:) :: flux,fluxt,flux2,flux2t,den,dent
  integer :: ic,ic0(4),i,iplot,m,n,iflag
  character*3 fldname(2)
  character*3 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  :: param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,angle,vtot
  real(p2)  :: amass,wx1,wy1,wz1,wx0,wy0,wz0,bc(3)
  real(p2)  :: uu(3),gamma,beq0,B,wei  !,vvmin,vvmax,dv
  integer   :: kind,ii,jj,kk,ip,jp,kp,nEk,nEa
  integer :: j,k,L,imc1,imc2,index,IERROR
  logical    file_exist
  real*8    :: xyz(3),p,q,w

  save ic,ic0
  data ic/0/,ic0/0,0,0,0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
 
  if(mod(simtype(2),100) == 66) then
     call output_DDtps_nc1file_66(kind)
     return
  endif
  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed



  param      = 0.
  param(1)   = nEch
  param(2)   = nalpha0
  param(3)   = Emin
  param(4)   = Emax
  param(5)   = amin
  param(6)   = amax
  param(7)   = vvmin(kind)
  param(8)   = vvmax(kind)
  param(9)   = dE
  param(10)  = dalpha
  param(11)  = dv(kind)
  param(21)  = utime
  param(22)  = uspeed
  param(23)  = ut0
  param(24)  = cspeed

  fldname = (/"DD "," E "/)

  write(ct,'(I3.3)')((istep-istep0)/ndiagfe - 1)/100000
  ic      = mod(((istep-istep0)/ndiagfe-1),100000) + 1   !max(mod(ic,20000),1) 
  ic0(kind)     = ic0(kind)+1

  if(kind==1) filename = 'DDP'//ct//'.nc'
  if(kind==2) filename = 'DDHe'//ct//'.nc'
  if(kind==3) filename = 'DDO'//ct//'.nc'
  if(kind==4) filename = 'DDe'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if(.not. file_exist  .or. ic ==0  .or. (ic0(kind)==1 .and. irun ==0)) then
      ic0_dd(kind) = ic - 1   
  endif
  ic = ic - ic0_dd(kind)

  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, ic /)
  if(mype==0 ) then
        icount = (/nEch, nalpha0, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )


          call check( nf90_def_dim(ncid, 'nE', nEch, dimid(3)) )
          call check( nf90_def_var(ncid, 'nE', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'eV') )

          call check( nf90_def_dim(ncid, 'nalpha', nalpha0, dimid(4)) )
          call check( nf90_def_var(ncid, 'nalpha', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, ' ') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )

          dimids=(/dimid(3),dimid(4),dimid(5)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(3), echannel) )
          call check( nf90_put_var(ncid, varid(4), achannel) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif


  allocate(flux(nEch,nalpha0),flux2(nEch,nalpha0),fluxt(nEch,nalpha0),&
           flux2t(nEch,nalpha0),den(nEch,nalpha0),dent(nEch,nalpha0))
  flux  = 0.
  flux2 = 0.
  den   = 0.

  do m=1,mblocks
	 if(kind ==1) TP = block(m)%TPi
	 if(kind ==2) TP = block(m)%TPHe
     if(kind ==3) TP = block(m)%TPO
     if(kind ==4) TP = block(m)%TPe
     dO  L = 1, TP%mi
         amass  = TP%amass

         p        = TP%qv(1,L)
         q        = TP%qv(2,L)
         w        = TP%qv(3,L)

	     wx1      = p +1.0 
         ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	     ip       = min(block(m)%nx1,ii+1)
	     wx1      = wx1 - ii
	     wx0      = 1.0 - wx1

	     wy1      = q +1.0 
         jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	     jp       = min(block(m)%ny1,jj+1)
	     wy1      = wy1 - jj
	     wy0      = 1.0 - wy1

         wz1      = w +1.0 
         kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	     kp       = min(block(m)%nz1,kk+1)
	     wz1      = wz1 - kk
	     wz0      = 1.0 - wz1

	 	 bx       = bcar(m)%vector(1,ii,jj,kk)
		 by       = bcar(m)%vector(2,ii,jj,kk)
		 bz       = bcar(m)%vector(3,ii,jj,kk)
		 B        = sqrt(bx**2+by**2+bz**2)

         uu(1:3)  = TP%qv(4:6,L)/amass
	     gamma    = sqrt(1.+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	     vpar     = (uu(1)*bx+uu(2)*by+uu(3)*bz)/max(sqrt(bx**2+by**2+bz**2),1.e-4)
		 vtot     = sqrt(uu(1)**2+uu(2)**2+uu(3)**2)
		 vper     = sqrt(abs(vtot**2 - vpar**2))
		 angle    = acos( max(min(abs(vpar)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2]
! to get equatorial pitch angle for dipole case

         if(dipole) then
		    xyz           = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
            x             = xyz(1)
            y             = xyz(2)
            z             = xyz(3)
			bc            = dipole3DB(x,y,z)
		    B             = sqrt(bc(1)**2 + bc(2)**2 + bc(3)**2)

		    xyz           = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !get magnetic field at equator
            x             = xyz(1)
            y             = xyz(2)
            z             = xyz(3)
			bc            = dipole3DB(x,y,z)
		    Beq0          = sqrt(bc(1)**2 + bc(2)**2 + bc(3)**2)

			vper          = vper * sqrt(Beq0/max(B,1.e-6))
			vpar          = vpar/max(abs(vpar),1.e-10_p2)*sqrt(abs(vtot**2-vper**2))
!		    angle         = asin( max(min(abs(vper)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2] alpha0
		    angle         = acos( max(min(vpar/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2]

!			if(ic == 1) then
!               TP%qv(7,L) = ee
!               TP%qv(8,L) = angle
!			endif

		 endif
		 nEk              = int(TP%qv(7,L)+ 0.01)   !e
		 nEa              = int(TP%qv(8,L)+ 0.01)
		 wei              = TP%qv(9,L)
         flux(nEk,nEa)    = flux(nEk,nEa) + angle * wei
		 flux2(nEk,nEa)   = flux2(nEk,nEa)+ angle**2 *wei
		 den(nEk,nEa)     = den(nEk,nEa) + wei
	 enddo
  enddo
  CALL MPI_ALLREDUCE(den,dent,size(den),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(flux,fluxt,size(flux),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(flux2,flux2t,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  fluxt  = flux2t/max(dent,1.e-8) - (fluxt/max(dent,1.e-4))**2 

  do iplot = 1, 1

     if(mype==0) call check( nf90_put_var(ncid, varid(7+iplot), fluxt, start = istart, &
                              count = icount) )

  enddo

  if(mype==0) call check( nf90_close(ncid) )

  deallocate(flux,fluxt,flux2t,flux2,den,dent)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_DDtps_nc1file

!-----------------------------------------------------------------
subroutine output_DDtps_nc1file_66(kind)  ! test particles
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  use vector_functions
  implicit none
  include 'mpif.h'
  type(TestParticle_type) :: TP
  type(bvector_type), dimension(:), pointer    :: bvt,dummy
  real,  pointer,dimension(:,:,:) :: flux,fluxt,flux2,flux2t,&
                                     fluxe,fluxet,fluxe2,fluxe2t,den,dent,den0,den0t
  real,  pointer,dimension(:,:,:) :: fluxpar,fluxpart,fluxper,fluxpert
  real,  pointer,dimension(:) :: ech,ach
  integer :: ic,ic0(4),i,iplot,m,n,iflag
  character*3 fldname(6)
  character*3 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(4),icount(4),istart(4),icount1(1),istart1(1)
  real  :: param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,angle,vtot,qx,qy,w1,w2,w3,w4,w11,w12,w21,w22
  real(p2)  :: amass,wx1,wy1,wz1,wx0,wy0,wz0,bc(3),e(3),b(3),bb
  real(p2)  :: uu(3),gamma,beq0,wei,charge,dangle,angle0,ee0,dee
  real(p2)  :: bbpar(3), bper1(3),bper2(3),epar,eper1,eper2,&
	           vper1,vper2,vper12,vper22,evdot

  integer   :: kind,ii,jj,kk,ip,jp,kp,nEk,nEa,mn
  integer :: j,k,L,imc1,imc2,index,IERROR
  logical    file_exist
  real*8    :: xyz(3),p,q,w

  save ic,ic0
  data ic/0/,ic0/0,0,0,0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

      subroutine bconcarB(c,d)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c,d
      end subroutine bconcarB

      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
 	  end subroutine bcovcarE


  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange

  call allocate_bvector(bvt,3)
  do m=1,mblocks
     bvt(m)%vector = bv0(m)%vector !+  bv1(m)%vector
  enddo
  
  call bconcarB(bvt,bcar)
  call bcovcarE(ev1,Ecar)
  call allocate_bvector(bvt,-3)


  allocate(ech(nEch),ach(nalpha0))
  do i = 1,nEch
     ech(i) = emin + (i-1.)*de
  enddo
  do i = 1,nalpha0
     ach(i) = 0. + (i-1.)*dalpha
  enddo

  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed



 
  param      = 0.
  param(1)   = nEch
  param(2)   = nalpha0
  param(3)   = Emin
  param(4)   = Emax
  param(5)   = amin
  param(6)   = amax
  param(7)   = vvmin(kind)
  param(8)   = vvmax(kind)
  param(9)   = dE
  param(10)  = dalpha
  param(11)  = dv(kind)
  param(21)  = utime
  param(22)  = uspeed
  param(23)  = ut0
  param(24)  = cspeed

  fldname = (/"DD ","Den","dA ", "dE ","par","per"/)

  write(ct,'(I3.3)')((istep-istep0)/ndiagfe - 1)/100000
  ic      = mod(((istep-istep0)/ndiagfe-1),100000) + 1   !max(mod(ic,20000),1) 
  ic0(kind)     = ic0(kind)+1

  if(kind==1) filename = 'DDP'//ct//'.nc'
  if(kind==2) filename = 'DDHe'//ct//'.nc'
  if(kind==3) filename = 'DDO'//ct//'.nc'
  if(kind==4) filename = 'DDe'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if(.not. file_exist  .or. ic ==0  .or. (ic0(kind)==1 .and. irun ==0)) then
      ic0_dd(kind) = ic - 1   
  endif
  ic = ic - ic0_dd(kind)

  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, 1, ic /)
  if(mype==0 ) then
        icount = (/nEch, nalpha0,n_TP_region, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'n_region', n_TP_region, dimid(2)) )
          call check( nf90_def_var(ncid, 'n_region', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nE', nEch, dimid(3)) )
          call check( nf90_def_var(ncid, 'nE', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'eV') )


          call check( nf90_def_dim(ncid, 'nalpha', nalpha0, dimid(4)) )
          call check( nf90_def_var(ncid, 'nalpha', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'degree') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )

          dimids=(/dimid(3),dimid(4),dimid(2),dimid(5)/)
          do iplot = 1, 6
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(3), ech) )
          call check( nf90_put_var(ncid, varid(4), ach) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 6
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif

  deallocate(ech,ach)
  allocate(flux(nEch,nalpha0,n_TP_region),flux2(nEch,nalpha0,n_TP_region),&
           fluxt(nEch,nalpha0,n_TP_region),&
           flux2t(nEch,nalpha0,n_TP_region),fluxe(nEch,nalpha0,n_TP_region),&
		   fluxe2(nEch,nalpha0,n_TP_region),fluxet(nEch,nalpha0,n_TP_region),&
           fluxe2t(nEch,nalpha0,n_TP_region),den(nEch,nalpha0,n_TP_region),&
		   dent(nEch,nalpha0,n_TP_region),&
		   den0(nEch,nalpha0,n_TP_region),den0t(nEch,nalpha0,n_TP_region))
  allocate(fluxpar(nEch,nalpha0,n_TP_region),fluxpart(nEch,nalpha0,n_TP_region),&
           fluxper(nEch,nalpha0,n_TP_region),fluxpert(nEch,nalpha0,n_TP_region))

  flux    = 0.
  flux2   = 0.
  fluxe   = 0.
  fluxe2  = 0.
  fluxpar = 0.
  fluxper = 0.
  den0    = 0.
  den     = 0.

  do m=1,mblocks
	 if(kind ==1) TP = block(m)%TPi
	 if(kind ==2) TP = block(m)%TPHe
     if(kind ==3) TP = block(m)%TPO
     if(kind ==4) TP = block(m)%TPe
     dO  L = 1, TP%mi
         amass    = TP%amass
		 mn       = TP%qv(9,L)
		 if(mn <1 .or. mn > n_TP_region) goto 772
		 wei      = 1.

         p        = TP%qv(1,L)
         q        = TP%qv(2,L)
         w        = TP%qv(3,L)

	     wx1      = p +1.0 
         ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	     ip       = min(block(m)%nx1,ii+1)
	     wx1      = wx1 - ii
	     wx0      = 1.0 - wx1

	     wy1      = q +1.0 
         jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	     jp       = min(block(m)%ny1,jj+1)
	     wy1      = wy1 - jj
	     wy0      = 1.0 - wy1

         wz1      = w +1.0 
         kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	     kp       = min(block(m)%nz1,kk+1)
	     wz1      = wz1 - kk
	     wz0      = 1.0 - wz1

         e    = Ecar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			    Ecar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
			    Ecar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
			    Ecar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
			    Ecar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			    Ecar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
			    Ecar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
			    Ecar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
         
		 b    = Bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			    Bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
			    Bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
			    Bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
			    Bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			    Bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
			    Bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
			    Bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
!	 	 bx       = bcar(m)%vector(1,ii,jj,kk)
!		 by       = bcar(m)%vector(2,ii,jj,kk)
!		 bz       = bcar(m)%vector(3,ii,jj,kk)
		 bx       = b(1)
		 by       = b(2)
		 bz       = b(3)
		 BB       = sqrt(bx**2+by**2+bz**2)

         uu(1:3)  = TP%qv(4:6,L)/amass
	     gamma    = sqrt(1.+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
		 ee       = (gamma-1.)*cspeed**2*amass * uT0  ! in eV
		 ee       = log10(ee)
         if(ifrelastic == 0) gamma   = 1.
         uu       = uu/gamma
		 vx       = uu(1)
		 vy       = uu(2)
		 vz       = uu(3)
		 vtot     = sqrt(uu(1)**2+uu(2)**2+uu(3)**2)



         bbpar(1)= bx/bb
		 bbpar(2)= by/bb
		 bbpar(3)= bz/bb
	     if((bb .ne. bb) .or. bb <= 1.e-4) bbpar = zaix
	     bper1  = a_cross_b(bbpar,xaix)
		 if(abs_a(bper1) <= 1.e-4_p2) bper1 = a_cross_b(bbpar,xaix)
		 bper1  = bper1/abs_a(bper1)
		 bper2  = a_cross_b(bbpar,bper1)
		 bper2  = bper2/abs_a(bper2)


         VPAR   = vx*bbpar(1)  + vy*bbpar(2)  + vz*bbpar(3)
         VPer1  = vx*bper1(1)  + vy*bper1(2)  + vz*bper1(3)
         VPer2  = vx*bper2(1)  + vy*bper2(2)  + vz*bper2(3)

         EPAR   = e(1)*bbpar(1)  + e(2)*bbpar(2)  + e(3)*bbpar(3)
         EPer1  = e(1)*bper1(1)  + e(2)*bper1(2)  + e(3)*bper1(3)
         EPer2  = e(1)*bper2(1)  + e(2)*bper2(2)  + e(3)*bper2(3)

		 evdot  = eper1 * vper1 + eper2*vper2
		 vper   = sqrt(abs(vtot**2 - vpar**2))
		 angle  = acos( max(min(vpar/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2]


! to get equatorial pitch angle for dipole case

         if(dipole) then
		    xyz           = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
            x             = xyz(1)
            y             = xyz(2)
            z             = xyz(3)
			bc            = dipole3DB(x,y,z)
		    BB            = sqrt(bc(1)**2 + bc(2)**2 + bc(3)**2)

            bbpar(1)= bc(1)/bb
		    bbpar(2)= bc(2)/bb
		    bbpar(3)= bc(3)/bb
	        if((bb .ne. bb) .or. bb <= 1.e-4) bbpar = zaix
	        bper1  = a_cross_b(bbpar,xaix)
		    if(abs_a(bper1) <= 1.e-4_p2) bper1 = a_cross_b(bbpar,xaix)
		    bper1  = bper1/abs_a(bper1)
		    bper2  = a_cross_b(bbpar,bper1)
		    bper2  = bper2/abs_a(bper2)


            VPAR   = vx*bbpar(1)  + vy*bbpar(2)  + vz*bbpar(3)
            VPer1  = vx*bper1(1)  + vy*bper1(2)  + vz*bper1(3)
            VPer2  = vx*bper2(1)  + vy*bper2(2)  + vz*bper2(3)

            EPAR   = e(1)*bbpar(1)  + e(2)*bbpar(2)  + e(3)*bbpar(3)
            EPer1  = e(1)*bper1(1)  + e(2)*bper1(2)  + e(3)*bper1(3)
            EPer2  = e(1)*bper2(1)  + e(2)*bper2(2)  + e(3)*bper2(3)

		    evdot  = eper1 * vper1 + eper2*vper2
		    vper   = sqrt(abs(vtot**2 - vpar**2))



		    xyz           = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !get magnetic field at equator
            x             = xyz(1)
            y             = xyz(2)
            z             = xyz(3)
			bc            = dipole3DB(x,y,z)
		    Beq0          = sqrt(bc(1)**2 + bc(2)**2 + bc(3)**2)

			vper          = vper * sqrt(Beq0/max(BB,1.e-6))
			vpar          = vpar/max(abs(vpar),1.e-10_p2)*sqrt(max(vtot**2-vper**2,1.e-10_p2))
		    angle         = acos( max(min(abs(vpar)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2]

			if(ic == 1) then
               TP%qv(7,L) = ee
               TP%qv(8,L) = angle
			endif
		 endif
         ee0    = TP%qv(7,L)
		 angle0 = TP%qv(8,L) 

!        if(ee >= emin .and. ee < emax ) then 

		    dee   = 10.**ee - 10.**ee0
			dangle= angle-angle0 

		    qx    = (ee0 -Emin) /de + 1.
		    qy    = (angle0 -0.) /dalpha  + 1.
            if(int(qx) < 1 .or. int(qy) < 1 .or. int(qx) >= nEch .or. int(qy) >= nalpha0) goto 771  

		    i     = min(max(int(qx),1),nEch-1)
		    j     = min(max(int(qy),1),nalpha0-1)
		    w2    = qx - i
		    w1    = 1.- w2
			w4    = qy - j    
			w3    = 1. -w4    

            w11   = w1*w3
            w12   = w1*w4
            w21   = w2*w3
            w22   = w2*w4
            flux(i,j,mn)        = flux(i,j,mn)     + w11* dangle * wei
            flux(i,j+1,mn)      = flux(i,j+1,mn)   + w12* dangle * wei 
            flux(i+1,j,mn)      = flux(i+1,j,mn)   + w21* dangle * wei 
            flux(i+1,j+1,mn)    = flux(i+1,j+1,mn) + w22* dangle * wei 

            flux2(i,j,mn)       = flux2(i,j,mn)     + w11* dangle**2 * wei
            flux2(i,j+1,mn)     = flux2(i,j+1,mn)   + w12* dangle**2 * wei 
            flux2(i+1,j,mn)     = flux2(i+1,j,mn)   + w21* dangle**2 * wei
            flux2(i+1,j+1,mn)   = flux2(i+1,j+1,mn) + w22* dangle**2 * wei

            fluxe2(i,j,mn)       = fluxe2(i,j,mn)     + w11* dee * wei
            fluxe2(i,j+1,mn)     = fluxe2(i,j+1,mn)   + w12* dee * wei 
            fluxe2(i+1,j,mn)     = fluxe2(i+1,j,mn)   + w21* dee * wei
            fluxe2(i+1,j+1,mn)   = fluxe2(i+1,j+1,mn) + w22* dee * wei


            fluxpar(i,j,mn)       = fluxpar(i,j,mn)     + w11* vpar*epar * wei
            fluxpar(i,j+1,mn)     = fluxpar(i,j+1,mn)   + w12* vpar*epar * wei 
            fluxpar(i+1,j,mn)     = fluxpar(i+1,j,mn)   + w21* vpar*epar * wei
            fluxpar(i+1,j+1,mn)   = fluxpar(i+1,j+1,mn) + w22* vpar*epar * wei


            fluxper(i,j,mn)       = fluxper(i,j,mn)     + w11* evdot * wei
            fluxper(i,j+1,mn)     = fluxper(i,j+1,mn)   + w12* evdot * wei 
            fluxper(i+1,j,mn)     = fluxper(i+1,j,mn)   + w21* evdot * wei
            fluxper(i+1,j+1,mn)   = fluxper(i+1,j+1,mn) + w22* evdot * wei

            den0(i,j,mn)          = den0(i,j,mn)     + w11 * wei 
            den0(i,j+1,mn)        = den0(i,j+1,mn)   + w12 * wei 
            den0(i+1,j,mn)        = den0(i+1,j,mn)   + w21 * wei
            den0(i+1,j+1,mn)      = den0(i+1,j+1,mn) + w22 * wei 

! -----------------------------------------------
771         continue
		    qx    = (ee -Emin) /de + 1.
		    qy    = (angle -0.) /dalpha  + 1.
            if(int(qx) < 1 .or. int(qy) < 1 .or. int(qx) >= nEch .or. int(qy) >= nalpha0) goto 772  
		    i     = min(max(int(qx),1),nEch-1)
		    j     = min(max(int(qy),1),nalpha0-1)
		    w2    = qx - i
		    w1    = 1.- w2
			w4    = qy - j    
			w3    = 1. -w4    

            w11   = w1*w3
            w12   = w1*w4
            w21   = w2*w3
            w22   = w2*w4


            den(i,j,mn)         = den(i,j,mn)     + w11 * wei 
            den(i,j+1,mn)       = den(i,j+1,mn)   + w12 * wei 
            den(i+1,j,mn)       = den(i+1,j,mn)   + w21 * wei
            den(i+1,j+1,mn)     = den(i+1,j+1,mn) + w22 * wei 
772         continue
!		 endif
	 enddo  ! end of L
  enddo  ! end of mblock
  CALL MPI_ALLREDUCE(den,dent,size(den),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(den0,den0t,size(den0),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(flux,fluxt,size(flux),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(flux2,flux2t,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  CALL MPI_ALLREDUCE(fluxe2,fluxe2t,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fluxpar,fluxpart,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fluxper,fluxpert,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)


  flux2t    = flux2t/max(den0t,1.e-8) - (fluxt/max(den0t,1.e-4))**2 
  fluxt     = fluxt/max(den0t,1.e-8)
  fluxe2t   = fluxe2t/max(den0t,1.e-8)
  fluxpart  = fluxpart/max(den0t,1.e-8)
  fluxpert  = fluxpert/max(den0t,1.e-8)

  dent(1,:,:)       = 2.* dent(1,:,:)
  dent(nEch,:,:)    = 2.* dent(nEch,:,:)
  dent(:,1,:)       = 2.* dent(:,1,:)
  dent(:,nalpha0,:) = 2.* dent(:,nalpha0,:)

  do iplot = 1, 1

     if(mype==0) call check( nf90_put_var(ncid, varid(7+iplot), flux2t, start = istart, &
                              count = icount) )
!     if(mype==0) call check( nf90_put_var(ncid, varid(8+iplot), flux2t, start = istart, &
!                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(8+iplot), dent, start = istart, &
                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(9+iplot), fluxt, start = istart, &
                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(10+iplot), fluxe2t, start = istart, &
                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(11+iplot), fluxpart, start = istart, &
                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(12+iplot), fluxpert, start = istart, &
                              count = icount) )

  enddo

  if(mype==0) call check( nf90_close(ncid) )  

  deallocate(flux,flux2,fluxt,flux2t,fluxe,fluxe2,fluxet,fluxe2t,den,dent,den0,den0t)
  deallocate(fluxpar,fluxpart,fluxper,fluxpert)



  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_DDtps_nc1file_66


!-----------------------------------------------------------------
subroutine DiagEdistribution_gyro_nc(rstart)
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use vector_functions
  use ionRingbeam_Case, only : vring
  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer     :: bv0_,bfld,bv1_,ev1_
  type(particle_type)    :: particle
  real(p2),pointer,dimension(:,:,:,:,:,:) :: fv,fvt
  real(p2),pointer,dimension(:) :: vpar,vper,csi
  real,pointer,dimension(:,:) :: num,numt
  real,pointer,dimension(:,:,:,:,:) :: plot
  real(p2),pointer,dimension(:,:,:) :: ff,fft

  real(p2) :: vminp,vmaxp,p,q,w,x,y,z
  real(p2) :: dcsi,dvpar,dvper,scale,rstart
  real(8)  :: xyz(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),vb(3),vb1(3),phase,det
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k
  real(p2) :: v1,v2,v3,ee,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,angle
  real(p2) :: bw(3),vw(3),vpar_,vper_ ,weight,dweight,gamma,amass,vpermax,dmc1,dmc2,dmc3, &
              ew(3),ve1(3),flux(5),abw,aew,frac

  integer  :: nv,nvper,ncsi,ic,m,ierror,L,Li,Le,kind,n,nbuff,ic0,iplot,&
              iflag,kinde,imc1,imc2,imc3,mkind,pkind,mn

  character*4 fldname(6)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(6),icount(6),istart(6),icount1(1),istart1(1)
  real  :: param(80)

  integer  :: ind(20),index
  logical file_exist,flagi,flage

  save ic,ic0
  data ic/0/,ic0/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE


     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB


  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange

  if(nr_fv < 1 ) return
! if(ic0_save(3) >=50) return
  kinde      = int(pvdomain(3))
  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed

  param(10)  = kinde
  flagi      =.false.
  flage      =.false.
  ind        = 0
  do m=1,kinde
     param(10+m)           =  int(pvdomain(3+m))	
     do n = 1, kinds
	    if(	int(param(10+m)) == n   ) ind(n)    = m  ! ions
	    if(	int(param(10+m)) == n+10) ind(n+10) = m  ! eles
	 enddo	 
	  
	 if(param(10+m) < 10) flagi = .true.
	 if(param(10+m) > 10) flage = .true.
  enddo

  do m=1,kinde
     if(param(10+m) <10 .and. param(10+m)>0 ) then
	    param(20+m) = sqrt(Tions(int(param(10+m)))) * vthi           !scales at param(21,22...)
		if((vring(1)+vring(2)) * vions(int(param(10+m))) > 0.) param(20+m) = 1.
     endif
     if(param(10+m) >10) then
	    param(20+m) = sqrt(Teles(int(param(10+m)-10)))
     endif
  enddo


  call allocate_bvector(bfld,3)
  call allocate_bvector(bv0_,3)
  call allocate_bvector(bv1_,3)
  call allocate_bvector(ev1_,3)
  do m= 1, mblocks
     bfld(m)%vector     = bv0(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv0_)
  call bconcarB(ev1,ev1_)
  do m= 1, mblocks
     bfld(m)%vector     = bv1(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv1_)
  call allocate_bvector(bfld,-3)

! -----------------------------------------------------------
  nv      =  201
  ncsi    =  101
  nvper   =  6
  vmaxp   =  5.
  vminp   = -5.
  

  dvpar   = (vmaxp-vminp)/real(nv-1)
  dcsi    = (pi*2)/real(ncsi-1)
  
  dvper   = 1.
  vpermax = 6.
  dvper   = vpermax/nvper


  allocate(fv(kinde,nv,ncsi,nvper,nr_fv,5),&
           fvt(kinde,nv,ncsi,nvper,nr_fv,5), &
		   vpar(nv),csi(ncsi),vper(nvper))
  allocate(num(kinde,nr_fv),numt(kinde,nr_fv))

! allocate(ff(n_e,n_alpha,nr_fv),fft(n_e,n_alpha,nr_fv))
! allocate(array_e(n_e),array_alpha(n_alpha),volum(nr_fv))

  do i = 1, nv
     vpar(i) = vminp + (i-1)*dvpar
  enddo

  do i = 1, ncsi
	 csi(i)  = (i-1) * dcsi
  enddo

  do i = 1, nvper
     vper(i) = 0. + (i-0.5)*dvper
  enddo

  fldname = (/"fvci","Je_1","Je_2","Jb_1","Jb_2","dfv "/)
  write(ct,'(I6.6)')(int((istep-1)*dt)/ndiagfv)*ndiagfv
  filename = 'fvcsi'//ct//'.nc'

  inquire(file=filename,exist=file_exist)
  if((.not. file_exist) .or. abs(rstart)<5*dt) then
      ic0_save(3) = 0 
  endif
  ic0_save(3) = ic0_save(3) +1
  ic          = ic0_save(3)

  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, 1, 1, 1, ic /)
  if(mype==0 ) then
        icount = (/kinde,nv, ncsi, nvper, nr_fv, 1 /)

!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'kind', kinde, dimid(1)) )
          call check( nf90_def_var(ncid, 'kind', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nvpar', nv, dimid(2)) )
          call check( nf90_def_var(ncid, 'nvpar', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'ncsi', ncsi, dimid(3)) )
          call check( nf90_def_var(ncid, 'ncsi', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nvper', nvper, dimid(4)) )
          call check( nf90_def_var(ncid, 'nvper', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nregion', nr_fv, dimid(5)) )
          call check( nf90_def_var(ncid, 'nregion', NF90_REAL, dimid(5), varid(5)) )
          call check( nf90_put_att(ncid,  varid(5), UNITS, 'Re') )



          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(6)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(6), varid(6)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(7)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(7), varid(7)) )
          call check( nf90_put_att(ncid,  varid(7), UNITS, 'non') )


          dimids2=(/dimid(5),dimid(4)/)
          call check( nf90_def_var(ncid, 'position', NF90_REAL, dimids2, varid(8)) )
          call check( nf90_put_att(ncid,  varid(8), UNITS, 'non') )


          dimids=(/dimid(1),dimid(2),dimid(3),dimid(4),dimid(5),dimid(6)/)
          do iplot = 1, 5
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(8+iplot)) )
           call check( nf90_put_att(ncid, varid(8+iplot), UNITS, 'non') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(7), param) )
          call check( nf90_put_var(ncid, varid(8), domain_vregion(:,1:6)) )

          call check( nf90_put_var(ncid, varid(2), vpar) )
          call check( nf90_put_var(ncid, varid(3), csi) )
          call check( nf90_put_var(ncid, varid(4), vper) )
          call check( nf90_put_var(ncid, varid(6), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(6)) )
          call check( nf90_put_var(ncid, varid(6), stime*aion, start = istart1) )
          do iplot = 1, 5
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(8+iplot)) )
          enddo
        endif
  endif


  fv      = 0.
 
! --------------------------------------------------------------
  do n=1,nr_fv
     do mn=1,kinds
        if(flage) particle = eles(mn)
        if(flagi) particle = ions(mn)
	    frac    = particle%frac
		amass   = particle%amass
		scale   = particle%vth
        index   = ind(mn + 10)
	    if(index == 0) goto 221
        do m=1,mblocks
           do L = 1,particle%block(m)%mi
!			  if(flage)  scale  = sqrt(Teles(mn))
              p     = particle%block(m)%qv(1,L)
              q     = particle%block(m)%qv(2,L)
              w     = particle%block(m)%qv(3,L)
			  weight  = frac
			  dweight = particle%block(m)%qv(13,L) * particle%block(m)%qv(14,L) * frac

              if(ghybrid) then
			     xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
		         x     = xyz(1)
		         y     = xyz(2)
		         z     = xyz(3)
              else
			     x     = p
				 y     = q
				 z     = w
              endif
			  if(x>= domain_vregion(n,1) .and. x< domain_vregion(n,2) .and. & 
			     y>= domain_vregion(n,3) .and. y< domain_vregion(n,4) .and. & 
			     z>= domain_vregion(n,5) .and. z< domain_vregion(n,6)) then

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1

                 vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 ve1   = ev1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             ev1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 ev1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 ev1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 ev1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             ev1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 ev1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 ev1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 
				 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
				 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2) 

                 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

				 uu       = particle%block(m)%qv(4:6,L)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma  /scale
				 ee       = 1.


				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)
				 abw      = max(sqrt(bw(1)**2+bw(2)**2+bw(3)**2),1.e-10_p2)

				 ew(1)    = x2b(1,1)*ve1(1) + x2b(1,2)*ve1(2)+x2b(1,3)*ve1(3)
				 ew(2)    = x2b(2,1)*ve1(1) + x2b(2,2)*ve1(2)+x2b(2,3)*ve1(3)
				 ew(3)    = x2b(3,1)*ve1(1) + x2b(3,2)*ve1(2)+x2b(3,3)*ve1(3)
				 aew      = max(sqrt(ew(1)**2+ew(2)**2+ew(3)**2),1.e-10_p2)

                 angle    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 angle    = mod(angle+pi2*2,pi2)

				 flux     = (/1._p2, qelectron*vw(1)*ew(1)/aew , &
				                     qelectron*(vw(2)*ew(2)+vw(3)*ew(3))/aew , &
                                     qelectron*vw(1)*bw(1)/abw , &
				                     qelectron*(vw(2)*bw(2)+vw(3)*bw(3))/abw /) * ee/mecell
                 
				 flux(2:5) = flux(2:5)/alpha

                 if(vpar_ < vminp .or. vpar_ >= vmaxp .or. vper_ >= vpermax ) goto 222 


                 imc1     = max(min(int((vpar_ - vminp)/dvpar + 1.),  nv-1),1)
		         dmc1     = imc1-(vpar_ - vminp)/dvpar
                 imc2     = min(int(angle/dcsi + 1.),ncsi-1)
		         dmc2     = imc2 - angle/dcsi
                 imc3     = min(int(vper_/dvper +1.),nvper)



                 fv(index,imc1,imc2,imc3,n,:)     = fv(index,imc1,imc2,imc3,n,:)    + &
				                                    dmc1*dmc2*dweight * flux 	          
                 fv(index,imc1,imc2+1,imc3,n,:)   = fv(index,imc1,imc2+1,imc3,n,:)  + &
				                                    dmc1*(1.-dmc2)*dweight * flux            
                 fv(index,imc1+1,imc2,imc3,n,:)   = fv(index,imc1+1,imc2,imc3,n,:)  + &
				                                    (1.-dmc1)*dmc2*dweight * flux          
                 fv(index,imc1+1,imc2+1,imc3,n,:) = fv(index,imc1+1,imc2+1,imc3,n,:)+ &
				                                    (1.-dmc1)*(1.-dmc2)*dweight * flux  	          

				 num(index,n)        = num(index,n) +  1. *frac

222              continue

			  endif	  			   
		   enddo
        enddo
221	    continue
     enddo
  enddo


  CALL MPI_ALLREDUCE(num,numt,kinde*nr_fv,mpi_real,MPI_SUM,MPI_COMM_WORLD,IERROR)
  if(p2==8) then
     call mpi_allreduce(fv,fvt,size(fv),mpi_double_precision,mpi_sum,MPI_COMM_WORLD, ierror)
  else
     call mpi_allreduce(fv,fvt,size(fv),mpi_real,mpi_sum,MPI_COMM_WORLD, ierror)
  endif

  do n=1,nr_fv
     do k=1,kinde
	    fvt(k,:,:,:,n,:)  = fvt(k,:,:,:,n,:)/domain_vregion(n,7)  !/max(numt(k,n),1.)
     enddo
  enddo

  do iplot = 1, 5
     if(mype==0) then
        allocate (plot(kinde,nv,ncsi,nvper,nr_fv))
		plot =  fvt(:,:,:,:,:,iplot)

        call check( nf90_put_var(ncid, varid(8+iplot), plot, start = istart, &
                              count = icount) )

        deallocate(plot)
     endif
  enddo

  if(mype==0) call check( nf90_close(ncid) )


  DEALLOCATE(vpar, fv,fvt,num,numt,csi,vper)
  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  call allocate_bvector(ev1_,-3)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine DiagEdistribution_gyro_nc

! -------------------------------------------------------------------------------------------

subroutine output_4fe_nc1file()
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use ionRingbeam_Case
  implicit none
  include 'mpif.h'
  real,  pointer,dimension(:,:,:,:) :: flux,fluxt
  integer :: ic,ic0,i,iplot,m,n,iflag
  character*3 fldname(2)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(11),dimids2(2),icount2(2),varid(30), &
             dimids(5),icount(5),istart(5),icount1(1),istart1(1)
  real  :: param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,angle,vtot
  real(p2)  :: xyz(3),p,q,w,amass
  integer :: j,k,L,imc1,imc2,index,IERROR,kinde,ind(20),ieeror,kind
  logical file_exist,flagi,flage

!  save ic,ic0
!  data ic/0/,ic0/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
 
  if(nr_flux == 0) return

  kinde     = int(pedomain(3))

  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed

  
  param(10)  = kinde
  flagi      =.false.
  flage      =.false.
  ind        = 0

  do m=1,kinde
     param(10+m)      =  int(pedomain(3+m))
	 ind(int(param(10+m))) = m
	 if(param(10+m) < 10) flagi = .true.
	 if(param(10+m) > 10) flage = .true.
  enddo

  do m=1,kinde
     if(param(10+m) <10) then
	    param(20+m) = sqrt(Tions(int(param(10+m)))) * vthi
		if((vring(1)+vring(2)) * vions(int(param(10+m))) > 0.) param(20+m) = 1.
     endif
     if(param(10+m) >10) param(20+m) = sqrt(Teles(int(param(10+m)-10)))
  enddo


  fldname = (/"PSD"," E "/)

  write(ct,'(I3.3)')(istep/ndiagfe - 1)/10000
  filename = 'fluxE'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if(.not. file_exist) then
      ic0_save(1) = 0 
  endif

  ic0_save(1) = ic0_save(1) +1
  ic          = ic0_save(1)

  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, 1, 1, ic /)
  if(mype==0 ) then
        icount = (/kinde,nr_flux, ne_flux, na_flux, 1 /)

!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'kind', kinde, dimid(1)) )
          call check( nf90_def_var(ncid, 'kind', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nr', nr_flux, dimid(2)) )
          call check( nf90_def_var(ncid, 'nr', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nE', nE_flux, dimid(3)) )
          call check( nf90_def_var(ncid, 'nE', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nalpha', na_flux, dimid(4)) )
          call check( nf90_def_var(ncid, 'nalpha', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )


          call check( nf90_def_dim(ncid, 'ver', 3, dimid(10)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(10), varid(10)) )
          call check( nf90_put_att(ncid,  varid(10), UNITS, 'non') )


          dimids2=(/dimid(10),dimid(2)/)
          call check( nf90_def_var(ncid, 'position', NF90_REAL, dimids2, varid(7)) )
          call check( nf90_put_att(ncid,  varid(7), UNITS, 'non') )


          dimids=(/dimid(1),dimid(2),dimid(3),dimid(4),dimid(5)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(7), array_r_flux) )
          call check( nf90_put_var(ncid, varid(3), array_e_flux) )
          call check( nf90_put_var(ncid, varid(4), array_alpha_flux) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif


  allocate(flux(kinde,nr_flux,nE_flux,na_flux),fluxt(kinde,nr_flux,nE_flux,na_flux))
  flux = 0.

  if(flagi) then
     do m=1,mblocks
        dO  L = 1, block(m)%mi
            index  =   ind(block(m)%ion(L)%kind)
		    if(index == 0) goto 221
		    kind   = block(m)%ion(L)%kind
            Factor = block(m)%ion(L)%ww * fraci(kind) !* qions(kind)
			if(fullf>=2) Factor = block(m)%ion(L)%w * Factor
            amass = aion*mions(kind)
	        p     = block(m)%ion(L)%p(1)
	        q     = block(m)%ion(L)%p(2) 
	        w     = block(m)%ion(L)%p(3)
			i     = int(p) + 1
			j     = int(q) + 1
			k     = int(w) + 1
            if(ghybrid) then
			     xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
		         x     = xyz(1)
		         y     = xyz(2)
		         z     = xyz(3)
            else
			     x     = p
				 y     = q
				 z     = w
            endif

			do n= 1, nr_flux
			   if(x>= domain_eregion(n,1) .and. x< domain_eregion(n,2) .and. & 
			      y>= domain_eregion(n,3) .and. y< domain_eregion(n,4) .and. & 
			      z>= domain_eregion(n,5) .and. z< domain_eregion(n,6)) then
	 	             bx    = bcar(m)%vector(1,i,j,k)
		             by    = bcar(m)%vector(2,i,j,k)
		             bz    = bcar(m)%vector(3,i,j,k)
	                 vx    = block(m)%ion(L)%v(1)/ amass
	                 vy    = block(m)%ion(L)%v(2)/ amass
	                 vz    = block(m)%ion(L)%v(3)/ amass
		             vpar  = (vx*bx+vy*by+vz*bz)/max(sqrt(bx**2+by**2+bz**2),1.e-4)
		             vtot  = sqrt(vx**2+vy**2+vz**2)
		       
			         angle = acos( max(min(vpar/max(vtot,1.e-6),1.),-1.))   ![0., pi/2]
		             ee    = log10(0.5*amass*max(vtot,1.e-10_p2)**2 * real(ut0) )
                     if(ee<emin_flux .or. ee>emax_flux) goto 881
                        imc1  = min(max(int(max(ee-emin_flux,0.)/de_flux+1.),1),ne_flux-1)
		                dmc1  = imc1-(ee - emin_flux)/de_flux  

                        imc2  = min(max(int(max(angle,0.)/da_flux+1.),1),na_flux-1)
		                dmc2  = imc2-(angle - 0.)/da_flux 
                        flux(index,n,imc1,imc2)     = flux(index,n,imc1,imc2) + &
						                              dmc1*dmc2*Factor *vtot**2*amass	          
                        flux(index,n,imc1,imc2+1)   = flux(index,n,imc1,imc2+1) + &
						                              dmc1*(1.-dmc2)*Factor	*vtot**2*amass          
                        flux(index,n,imc1+1,imc2)   = flux(index,n,imc1+1,imc2) + &
						                              (1.-dmc1)*dmc2*Factor	*vtot**2*amass          
                        flux(index,n,imc1+1,imc2+1) = flux(index,n,imc1+1,imc2+1) + &
						                              (1.-dmc1)*(1.-dmc2)*Factor*vtot**2*amass	          
881                  continue
				endif
			enddo
221         continue

	 enddo
  enddo
  endif

! -----------------------------------------------------------------------

  if(flage) then
     do m=1,mblocks
        dO  L = 1, block(m)%me
            index  =   ind(block(m)%ele(L)%kind + 10)
		    if(index == 0) goto 222
		    kind   = block(m)%ele(L)%kind
            Factor = block(m)%ele(L)%ww * frace(kind) !* qions(kind)
			if(fullf>=2) Factor = block(m)%ele(L)%w * Factor
            amass = aelectron
	        p     = block(m)%ele(L)%p(1)
	        q     = block(m)%ele(L)%p(2) 
	        w     = block(m)%ele(L)%p(3)
			i     = int(p) + 1
			j     = int(q) + 1
			k     = int(w) + 1
            if(ghybrid) then
			     xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
		         x     = xyz(1)
		         y     = xyz(2)
		         z     = xyz(3)
            else
			     x     = p
				 y     = q
				 z     = w
            endif

			do n= 1, nr_flux
			   if(x>= domain_eregion(n,1) .and. x< domain_eregion(n,2) .and. & 
			      y>= domain_eregion(n,3) .and. y< domain_eregion(n,4) .and. & 
			      z>= domain_eregion(n,5) .and. z< domain_eregion(n,6)) then
	 	             bx    = bcar(m)%vector(1,i,j,k)
		             by    = bcar(m)%vector(2,i,j,k)
		             bz    = bcar(m)%vector(3,i,j,k)
	                 vx    = block(m)%ele(L)%v(1)/ amass
	                 vy    = block(m)%ele(L)%v(2)/ amass
	                 vz    = block(m)%ele(L)%v(3)/ amass
		             vpar  = (vx*bx+vy*by+vz*bz)/max(sqrt(bx**2+by**2+bz**2),1.e-4)
		             vtot  = sqrt(vx**2+vy**2+vz**2)
		       
			         angle = acos( max(min(vpar/max(vtot,1.e-6),1.),-1.))   ![0., pi/2]
		             ee    = log10(0.5*amass*max(vtot,1.e-10_p2)**2 * real(ut0) )
                     if(ee<emin_flux .or. ee>emax_flux) goto 882
                        imc1  = min(max(int(max(ee-emin_flux,0.)/de_flux+1.),1),ne_flux-1)
		                dmc1  = imc1-(ee - emin_flux)/de_flux  

                        imc2  = min(max(int(max(angle,0.)/da_flux+1.),1),na_flux-1)
		                dmc2  = imc2-(angle - 0.)/da_flux 
                        flux(index,n,imc1,imc2)     = flux(index,n,imc1,imc2) + dmc1*dmc2*Factor *vtot**2*amass	          
                        flux(index,n,imc1,imc2+1)   = flux(index,n,imc1,imc2+1) + &
						                              dmc1*(1.-dmc2)*Factor	*vtot**2*amass          
                        flux(index,n,imc1+1,imc2)   = flux(index,n,imc1+1,imc2) + &
						                              (1.-dmc1)*dmc2*Factor	*vtot**2*amass          
                        flux(index,n,imc1+1,imc2+1) = flux(index,n,imc1+1,imc2+1) + &
						                              (1.-dmc1)*(1.-dmc2)*Factor*vtot**2*amass	          
882                  continue
				endif
			enddo
222         continue

	     enddo
      enddo
  endif

  CALL MPI_ALLREDUCE(flux,fluxt,size(flux),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  fluxt(:,:,1,:)         = 2.*fluxt(:,:,1,:)
  fluxt(:,:,ne_flux,:)   = 2.*fluxt(:,:,ne_flux,:)
  fluxt(:,:,:,1)         = 2.*fluxt(:,:,:,1)
  fluxt(:,:,:,na_flux)   = 2.*fluxt(:,:,:,na_flux)

  do iplot = 1, 1
     if(mype==0) call check( nf90_put_var(ncid, varid(7+iplot), fluxt, start = istart, &
                              count = icount) )

  enddo

  if(mype==0) call check( nf90_close(ncid) )

  deallocate(flux,fluxt)

 
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_4fe_nc1file

! -------------------------------------------------------------------------------------------

subroutine output_4fv_nc1file()  !f(v||,vper),f(v||,vper1),f(vper1,vper2)
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use vector_functions
  use ionRingbeam_Case

  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer    :: bv0_,bfld,bv1_
  real,pointer,dimension(:,:,:,:,:) :: fv,fvt
  real,pointer,dimension(:,:,:,:) :: fv4,fv4t
! real(p2),pointer,dimension(:) :: vpar,vper,vper1,vper2
! real,pointer,dimension(:,:) :: num,numt
  real,pointer,dimension(:,:,:,:) :: plot
  real,pointer,dimension(:,:,:) :: plot4

  real(p2) :: vminp,vmaxp,p,q,w,x,y,z,vmaxper
  real(p2) :: dcsi,dvpar,dvper,dvper1,scale
  real(8)  :: xyz(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),vb(3),vb1(3),phase,det
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k,s
  real(p2) :: v1,v2,v3,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,ww0,ww1,angle
  real(p2) :: bw(3),vw(3),vpar_,vper_ ,vper1_,vper2_,weight,gamma,dweight, &
              dmc1,dmc2,dmc3,dmc4,ee,factor,amass

  integer  :: nv,ic,m,ierror,L,kind,n,nbuff,ic0,iplot,iflag,kinde,ic0_d
  integer  :: ind(20),index

  character*10 fldname(6)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(15),dimids2(2),icount2(2),varid(30), &
             dimids(5),icount(5),istart(5),icount1(1),istart1(1), &
             dimids4(4),icount4(4),istart4(4),imc1,imc2,imc3,imc4

  real  :: param(80)

  logical file_exist,flagi,flage

  save ic,ic0,ic0_d
  data ic/0/,ic0/0/,ic0_d/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(nr_fv < 1 ) return
 
  call allocate_bvector(bfld,3)
  call allocate_bvector(bv0_,3)
  call allocate_bvector(bv1_,3)
  do m= 1, mblocks
     bfld(m)%vector     = bv0(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv0_)
  do m= 1, mblocks
     bfld(m)%vector     = bv1(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv1_)
  call allocate_bvector(bfld,-3)
  !------------------------
! pvdomain :  1-- number of region; 2--- istep beginning to make record ; 3: number of species to be recorded 4-9 index of 3, the max number is 6

  if(nr_fv == 0) return

  kinde     = int(pvdomain(3))

! --------------------------------------------------------------
  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed



  
  param(10)  = kinde
  flagi      =.false.
  flage      =.false.
  ind        = 0
  do m=1,kinde
     param(10+m)      =  int(pvdomain(3+m))
	 ind(int(param(10+m))) = m
	 if(param(10+m) < 10) flagi = .true.
	 if(param(10+m) > 10) flage = .true.
  enddo

  do m=1,kinde
     if(param(10+m) <10) then
	    param(20+m) = sqrt(Tions(int(param(10+m)))) * vthi
		if((vring(1)+vring(2)) * vions(int(param(10+m))) > 0.) param(20+m) = 1.
     endif
     if(param(10+m) >10) param(20+m) = sqrt(Teles(int(param(10+m)-10)))
  enddo

  fldname = (/"fvpar-vper","fvpar-vpe1","fvpe1-vpe2", &
              "fvpar     ","fvper     ","fvper1    " /)

  write(ct,'(I3.3)')(istep/ndiagfv - 1)/10000

  filename = 'fv'//ct//'.nc'

  inquire(file=filename,exist=file_exist)
  if(.not. file_exist) then
      ic0_save(2) = 0 
  endif

  ic0_save(2) = ic0_save(2) +1
  ic          = ic0_save(2)



  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, 1, 1, ic /)   !(nkind,nr,vpar,vper,t)
  istart4 = (/1, 1, 1, ic /)   !(nkind,nr,vpar,vper,t)
  if(mype==0 ) then
        icount  = (/kinde,nr_fv, nx_fv, ny_fv, 1 /)
        icount4 = (/kinde,nr_fv, nx_fv, 1 /)

!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'kind', kinde, dimid(1)) )
          call check( nf90_def_var(ncid, 'kind', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nr', nr_fv, dimid(2)) )
          call check( nf90_def_var(ncid, 'nr', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nvpar', nx_fv, dimid(3)) )
          call check( nf90_def_var(ncid, 'nvpar', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nvper', nz_fv, dimid(4)) )
          call check( nf90_def_var(ncid, 'nvper', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nvper1', ny_fv, dimid(5)) )
          call check( nf90_def_var(ncid, 'nvper1', NF90_REAL, dimid(5), varid(5)) )
          call check( nf90_put_att(ncid,  varid(5), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(6)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(6), varid(6)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(7)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(7), varid(7)) )
          call check( nf90_put_att(ncid,  varid(7), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'ver', 3, dimid(13)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(13), varid(13)) )
          call check( nf90_put_att(ncid,  varid(13), UNITS, 'non') )

          dimids2 = (/dimid(13),dimid(2)/)
          call check( nf90_def_var(ncid, 'position', NF90_REAL, dimids2, varid(8)) )
          call check( nf90_put_att(ncid,  varid(8), UNITS, 'non') )

          dimids  = (/dimid(1),dimid(2),dimid(3),dimid(4),dimid(6)/)
          do iplot = 1, 3
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(8+iplot)) )
           call check( nf90_put_att(ncid, varid(8+iplot), UNITS, 'flux') )
          enddo

          dimids4 = (/dimid(1),dimid(2),dimid(3),dimid(6)/)
          do iplot = 1, 3
           call check( nf90_def_var(ncid, fldname(iplot+3), NF90_REAL, dimids, varid(11+iplot)) )
           call check( nf90_put_att(ncid, varid(11+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(7), param) )
          call check( nf90_put_var(ncid, varid(8), array_r_fv) )
          call check( nf90_put_var(ncid, varid(3), array_vx_fv) )
          call check( nf90_put_var(ncid, varid(4), array_vz_fv) )
          call check( nf90_put_var(ncid, varid(5), array_vy_fv) )
          call check( nf90_put_var(ncid, varid(6), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(6)) )
          call check( nf90_put_var(ncid, varid(6), stime*aion, start = istart1) )
          do iplot = 1, 3
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(8+iplot)) )
          enddo
          do iplot = 1, 3
             call check( nf90_inq_varid(ncid, fldname(iplot+3), varid(11+iplot)) )
          enddo
        endif
  endif

  allocate(fv(kinde,3,nr_fv,nx_fv,ny_fv),fvt(kinde,3,nr_fv,nx_fv,ny_fv))
  allocate(fv4(kinde,3,nr_fv,nx_fv),fv4t(kinde,3,nr_fv,nx_fv))

  fv    = 0.
  fv4   = 0.
  if(flage) then
     do m=1,mblocks
           do L = 1,block(m)%me
              index   = ind(block(m)%ele(L)%kind + 10)
			  if(index == 0) goto 221
			  
			  scale  = sqrt(Teles(block(m)%ele(L)%kind))

              p      = block(m)%ele(L)%p(1)
              q      = block(m)%ele(L)%p(2)
              w      = block(m)%ele(L)%p(3)
			  weight = block(m)%ele(L)%ww
			  Factor = block(m)%ele(L)%ww * frace(block(m)%ele(L)%kind)
		  	  if(fullf>=2) Factor = block(m)%ele(L)%w * Factor
              if(ghybrid) then
			     xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
		         x     = xyz(1)
		         y     = xyz(2)
		         z     = xyz(3)
              else
			     x     = p
				 y     = q
				 z     = w
              endif
			  do n= 1, nr_fv
			     if(x< domain_vregion(n,1) .or. x>= domain_vregion(n,2) .or. & 
			        y< domain_vregion(n,3) .or. y>= domain_vregion(n,4) .or. & 
			        z< domain_vregion(n,5) .or. z>= domain_vregion(n,6)) goto 222

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1

                 vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
				 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2) 
                 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

				 uu       = block(m)%ele(L)%v(:)/aelectron
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu /gamma /scale  ! scale factor here
!				 vv       = uu /scale 

				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)

				 vpar_    = vw(1)
				 vper_    = max(sqrt(vw(2)**2+vw(3)**2),1.e-6)
				 vper1_   = vw(2)
				 vper2_   = vw(3)
				 ee       = 1.   ! differiential distribution

                 if(vpar_<vmin_fv .or. vpar_>=vmax_fv .or. &
                    vper1_<vmin_fv .or. vper1_>=vmax_fv .or. &
                    vper2_<vmin_fv .or. vper2_>=vmax_fv .or. &
                    vper_ <vmin_fv*0 .or. vper_>=vmax_fv ) goto 222

                        imc1  = min(max(int(max(vpar_-vmin_fv,0.)/dvx_fv+1.),1),nx_fv-1)
		                dmc1  = imc1-(vpar_ - vmin_fv)/dvx_fv

                        imc2  = min(max(int(max(vper_-vmin_fv*0,0.)/dvz_fv+1.),1),nz_fv-1)
		                dmc2  = imc2-(vper_ - vmin_fv*0)/dvz_fv

                        imc3  = min(max(int(max(vper1_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc3  = imc3-(vper1_ - vmin_fv)/dvy_fv

                        imc4  = min(max(int(max(vper2_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc4  = imc4-(vper2_ - vmin_fv)/dvy_fv
!                       vpar - vper ------
                        fv(index,1,n,imc1,imc2)     = fv(index,1,n,imc1,imc2)    + dmc1*dmc2*Factor * ee /vper_	          
                        fv(index,1,n,imc1,imc2+1)   = fv(index,1,n,imc1,imc2+1)  + dmc1*(1.-dmc2)*Factor * ee /vper_	          
                        fv(index,1,n,imc1+1,imc2)   = fv(index,1,n,imc1+1,imc2)  + (1.-dmc1)*dmc2*Factor * ee /vper_	          
                        fv(index,1,n,imc1+1,imc2+1) = fv(index,1,n,imc1+1,imc2+1)+ (1.-dmc1)*(1.-dmc2)*Factor * ee /vper_	          
!                       vpar - vper1 ------

                        fv(index,2,n,imc1,imc3)     = fv(index,2,n,imc1,imc3)    + dmc1*dmc3*Factor * ee	          
                        fv(index,2,n,imc1,imc3+1)   = fv(index,2,n,imc1,imc3+1)  + dmc1*(1.-dmc3)*Factor * ee	          
                        fv(index,2,n,imc1+1,imc3)   = fv(index,2,n,imc1+1,imc3)  + (1.-dmc1)*dmc3*Factor * ee	          
                        fv(index,2,n,imc1+1,imc3+1) = fv(index,2,n,imc1+1,imc3+1)+ (1.-dmc1)*(1.-dmc3)*Factor * ee	          

!                       vper1 - vper2 ------
    
                        fv(index,3,n,imc3,imc4)     = fv(index,3,n,imc3,imc4)    + dmc3*dmc4*Factor * ee	          
                        fv(index,3,n,imc3,imc4+1)   = fv(index,3,n,imc3,imc4+1)  + dmc3*(1.-dmc4)*Factor	 * ee          
                        fv(index,3,n,imc3+1,imc4)   = fv(index,3,n,imc3+1,imc4)  + (1.-dmc3)*dmc4*Factor	 * ee          
                        fv(index,3,n,imc3+1,imc4+1) = fv(index,3,n,imc3+1,imc4+1)+ (1.-dmc3)*(1.-dmc4)*Factor * ee	          

                        fv4(index,1,n,imc1)     = fv4(index,1,n,imc1)    + dmc1*Factor * ee 	          
                        fv4(index,1,n,imc1+1)   = fv4(index,1,n,imc1+1)  + (1.-dmc1)*Factor * ee 	          

                        fv4(index,2,n,imc2)     = fv4(index,2,n,imc2)    + dmc2*Factor * ee /vper_	          
                        fv4(index,2,n,imc2+1)   = fv4(index,2,n,imc2+1)  + (1.-dmc2)*Factor * ee /vper_	          

                        fv4(index,3,n,imc3)     = fv4(index,3,n,imc3)    + dmc3*Factor * ee 	          
                        fv4(index,3,n,imc3+1)   = fv4(index,3,n,imc3+1)  + (1.-dmc3)*Factor * ee 	          

222             continue
			  enddo	  			   
221           continue
		  enddo
     enddo
  endif

! -----------------------------------------------------------------------
  if(flagi) then
     do m=1,mblocks
           do L = 1,block(m)%mi
              index   = ind(block(m)%ion(L)%kind)
			  if(index == 0) goto 321
			  
			  scale  = sqrt(Tions(block(m)%ion(L)%kind)) * vthi
			  kind   = block(m)%ion(L)%kind
			  if(vring(1)*vions(kind)+vring(2)*vions(kind) > 0.) scale = 1.
			  amass  = aion * mions(block(m)%ion(L)%kind)

              p      = block(m)%ion(L)%p(1)
              q      = block(m)%ion(L)%p(2)
              w      = block(m)%ion(L)%p(3)
			  weight = block(m)%ion(L)%ww
			  Factor = block(m)%ion(L)%ww * fraci(block(m)%ion(L)%kind)
			  if(fullf>=2) Factor = block(m)%ion(L)%w * Factor

			  if(df) Factor = 1.0
              if(ghybrid) then
			     xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
		         x     = xyz(1)
		         y     = xyz(2)
		         z     = xyz(3)
              else
			     x     = p
				 y     = q
				 z     = w
              endif
			  do n= 1, nr_fv
			  if(x>= domain_vregion(n,1) .and. x< domain_vregion(n,2) .and. & 
			     y>= domain_vregion(n,3) .and. y< domain_vregion(n,4) .and. & 
			     z>= domain_vregion(n,5) .and. z< domain_vregion(n,6)) then

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1

                 vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2) 
                 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

				 uu       = block(m)%ion(L)%v(:)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu /gamma /scale  ! scale factor here
!				 vv       = uu /scale 

				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = max(sqrt(vw(2)**2+vw(3)**2),1.e-6)
				 vper1_   = vw(2)
				 vper2_   = vw(3)
				 ee       = 1.   ! differiential distribution
                 if(vpar_<vmin_fv .or. vpar_>=vmax_fv .or. &
                    vper1_<vmin_fv .or. vper1_>=vmax_fv .or. &
                    vper2_<vmin_fv .or. vper2_>=vmax_fv .or. &
                    vper_ <vmin_fv*0 .or. vper_>=vmax_fv ) goto 322

                        imc1  = min(max(int(max(vpar_-vmin_fv,0.)/dvx_fv+1.),1),nx_fv-1)
		                dmc1  = imc1-(vpar_ - vmin_fv)/dvx_fv

                        imc2  = min(max(int(max(vper_-vmin_fv*0,0.)/dvz_fv+1.),1),nz_fv-1)
		                dmc2  = imc2-(vper_ - vmin_fv*0)/dvz_fv

                        imc3  = min(max(int(max(vper1_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc3  = imc3-(vper1_ - vmin_fv)/dvy_fv

                        imc4  = min(max(int(max(vper2_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc4  = imc4-(vper2_ - vmin_fv)/dvy_fv
!                       vpar - vper ------
                        fv(index,1,n,imc1,imc2)     = fv(index,1,n,imc1,imc2)    + dmc1*dmc2*Factor * &
						                              ee /max(vper_,1.e-6)	          
                        fv(index,1,n,imc1,imc2+1)   = fv(index,1,n,imc1,imc2+1)  + dmc1*(1.-dmc2)*Factor * &
						                              ee /max(vper_,1.e-6)	          
                        fv(index,1,n,imc1+1,imc2)   = fv(index,1,n,imc1+1,imc2)  + (1.-dmc1)*dmc2*Factor * &
						                              ee /max(vper_,1.e-6)	          
                        fv(index,1,n,imc1+1,imc2+1) = fv(index,1,n,imc1+1,imc2+1)+ (1.-dmc1)*(1.-dmc2)*Factor * &
						                              ee /max(vper_,1.e-6)	          
!                       vpar - vper1 ------

                        fv(index,2,n,imc1,imc3)     = fv(index,2,n,imc1,imc3)    + dmc1*dmc3*Factor * ee	          
                        fv(index,2,n,imc1,imc3+1)   = fv(index,2,n,imc1,imc3+1)  + dmc1*(1.-dmc3)*Factor * ee	          
                        fv(index,2,n,imc1+1,imc3)   = fv(index,2,n,imc1+1,imc3)  + (1.-dmc1)*dmc3*Factor * ee	          
                        fv(index,2,n,imc1+1,imc3+1) = fv(index,2,n,imc1+1,imc3+1)+ (1.-dmc1)*(1.-dmc3)*Factor * ee	          

!                       vper1 - vper2 ------
    
                        fv(index,3,n,imc3,imc4)     = fv(index,3,n,imc3,imc4)    + dmc3*dmc4*Factor * ee	          
                        fv(index,3,n,imc3,imc4+1)   = fv(index,3,n,imc3,imc4+1)  + dmc3*(1.-dmc4)*Factor	 * ee          
                        fv(index,3,n,imc3+1,imc4)   = fv(index,3,n,imc3+1,imc4)  + (1.-dmc3)*dmc4*Factor	 * ee          
                        fv(index,3,n,imc3+1,imc4+1) = fv(index,3,n,imc3+1,imc4+1)+ (1.-dmc3)*(1.-dmc4)*Factor * ee	          

                        fv4(index,1,n,imc1)     = fv4(index,1,n,imc1)    + dmc1*Factor * ee 	          
                        fv4(index,1,n,imc1+1)   = fv4(index,1,n,imc1+1)  + (1.-dmc1)*Factor * ee 	          

                        fv4(index,2,n,imc2)     = fv4(index,2,n,imc2)    + dmc2*Factor * ee /vper_	          
                        fv4(index,2,n,imc2+1)   = fv4(index,2,n,imc2+1)  + (1.-dmc2)*Factor * ee /vper_	          

                        fv4(index,3,n,imc3)     = fv4(index,3,n,imc3)    + dmc3*Factor * ee 	          
                        fv4(index,3,n,imc3+1)   = fv4(index,3,n,imc3+1)  + (1.-dmc3)*Factor * ee 	          

322              continue
!                endif

			  endif
			  enddo	  			   
321           continue
		  enddo
     enddo
  endif



  CALL MPI_ALLREDUCE(fv,fvt,size(fv),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fv4,fv4t,size(fv4),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  do iplot = 1, 3
     if(mype==0) then
        allocate(plot(kinde,nr_fv,nx_fv,ny_fv))
		plot = fvt(:,iplot,:,:,:)
        call check( nf90_put_var(ncid, varid(8+iplot), plot, start = istart, &
                              count = icount) )
        deallocate(plot)
     endif
  enddo

  do iplot = 1, 3
     if(mype==0) then
        allocate(plot4(kinde,nr_fv,nx_fv))
		plot4 = fv4t(:,iplot,:,:)
        call check( nf90_put_var(ncid, varid(11+iplot), plot4, start = istart4, &
                              count = icount4) )
        deallocate(plot4)
     endif
  enddo

  if(mype==0) call check( nf90_close(ncid) )


  DEALLOCATE(fv,fvt,fv4,fv4t)
  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_4fv_nc1file


! -------------------------------------------------------------------------------------------
!f(v||,vper),f(v||,vper1),f(vper1,vper2)
!
subroutine output_4fv_nc1file66(rstart)  
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use vector_functions
  use ionRingbeam_Case

  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer    :: bv0_,bfld,bv1_,ev1_
  type(particle_type)    :: particle
  real,pointer,dimension(:,:,:,:,:) :: fv,fvt
  real,pointer,dimension(:,:,:,:) :: fv4,fv4t
! real(p2),pointer,dimension(:) :: vpar,vper,vper1,vper2
! real,pointer,dimension(:,:) :: num,numt
  real,pointer,dimension(:,:,:,:) :: plot
  real,pointer,dimension(:,:,:) :: plot4

  real(p2) :: vminp,vmaxp,p,q,w,x,y,z,vmaxper,rstart
  real(p2) :: dcsi,dvpar,dvper,dvper1,scale
  real(8)  :: xyz(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),vb(3),vb1(3),ve1(3),phase,det
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k,s
  real(p2) :: v1,v2,v3,aew,abw,vdb1,vdb2,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,ww0,ww1,angle
  real(p2) :: bw(3),vw(3),ew(3),vpar_,vper_ ,vper1_,vper2_,weight,gamma,dweight, &
              dmc1,dmc2,dmc3,dmc4,ee,factor,dfactor,amass,vde1,vde2,vdb,vdEB(4),frac

  integer  :: nv,ic,m,ierror,L,kind,n,nbuff,ic0,iplot,iflag,kinde,ic0_d
  integer  :: ind(20),index,mkind,pkind,mn

  character*10 fldname(10)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(18),dimids2(2),icount2(2),varid(30), &
             dimids(5),icount(5),istart(5),icount1(1),istart1(1), &
             dimids4(4),icount4(4),istart4(4),imc1,imc2,imc3,imc4

  real  :: param(80)

  logical file_exist,flagi,flage

  save ic,ic0,ic0_d
  data ic/0/,ic0/0/,ic0_d/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(nr_fv < 1 ) return
! if(ic0_save(2) >= 50) return 
 
  call allocate_bvector(bfld,3)
  call allocate_bvector(bv0_,3)
  call allocate_bvector(bv1_,3)
  call allocate_bvector(ev1_,3)
  do m= 1, mblocks
     bfld(m)%vector     = bv0(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv0_)
  do m= 1, mblocks
     bfld(m)%vector     = bv1(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv1_)
  call bcovcarE(ev1, ev1_)
  call allocate_bvector(bfld,-3)
  !------------------------
! pvdomain :  1-- number of region; 2--- istep beginning to make record ; 3: number of species to be recorded 4-9 index of 3, the max number is 6

  
  kinde     = int(pvdomain(3))

! --------------------------------------------------------------
  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed



  
  param(10)  = kinde
  flagi      =.false.
  flage      =.false.
  ind        = 0
  do m=1,kinde
     param(10+m)      =  int(pvdomain(3+m))
	 ind(int(param(10+m))) = m
	 if(param(10+m) < 10) flagi = .true.
	 if(param(10+m) > 10) flage = .true.
  enddo

  do m=1,kinde
     if(param(10+m) <10) then
	    param(20+m) = sqrt(Tions(int(param(10+m)))) * vthi
		if((vring(1)+vring(2)) * vions(int(param(10+m))) > 0.) param(20+m) = 1.
     endif
     if(param(10+m) >10) param(20+m) = sqrt(Teles(int(param(10+m)-10)))
  enddo

  fldname = (/"fvpar-vper","fvpar-vpe1","fvpe1-vpe2", &
              "fvpar     ","fvper     ","fvper1    " ,&
              "JdEpar    ","JdEper    ","JdBpar    ","JdBper    " /)

!  write(ct,'(I3.3)')(istep/ndiagfv - 1)/10000
!  select several time step to get the distribution data e.g.50 

! write(ct,'(I6.6)')(int((istep-1)*dt)/10)*10
  write(ct,'(I6.6)')(int((istep-1)*dt)/ndiagfv)*ndiagfv

  filename = 'fv'//ct//'.nc'



  inquire(file=filename,exist=file_exist)
  if((.not. file_exist) .or. abs(rstart)< 5* dt ) then
      ic0_save(2) = 0 
  endif

  ic0_save(2) = ic0_save(2) +1
  ic          = ic0_save(2)


  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, 1, 1, ic /)   !(nkind,nr,vpar,vper,t)
  istart4 = (/1, 1, 1, ic /)   !(nkind,nr,vpar,vper,t)
  if(mype==0 ) then
        icount  = (/kinde,nr_fv, nx_fv, ny_fv, 1 /)
        icount4 = (/kinde,nr_fv, nx_fv, 1 /)

!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'kind', kinde, dimid(1)) )
          call check( nf90_def_var(ncid, 'kind', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nr', nr_fv, dimid(2)) )
          call check( nf90_def_var(ncid, 'nr', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nvpar', nx_fv, dimid(3)) )
          call check( nf90_def_var(ncid, 'nvpar', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nvper', nz_fv, dimid(4)) )
          call check( nf90_def_var(ncid, 'nvper', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nvper1', ny_fv, dimid(5)) )
          call check( nf90_def_var(ncid, 'nvper1', NF90_REAL, dimid(5), varid(5)) )
          call check( nf90_put_att(ncid,  varid(5), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(6)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(6), varid(6)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(7)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(7), varid(7)) )
          call check( nf90_put_att(ncid,  varid(7), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'ver', 3, dimid(13)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(13), varid(13)) )
          call check( nf90_put_att(ncid,  varid(13), UNITS, 'non') )

          dimids2 = (/dimid(13),dimid(2)/)
          call check( nf90_def_var(ncid, 'position', NF90_REAL, dimids2, varid(8)) )
          call check( nf90_put_att(ncid,  varid(8), UNITS, 'non') )

          dimids  = (/dimid(1),dimid(2),dimid(3),dimid(4),dimid(6)/)
          do iplot = 1, 3
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(8+iplot)) )
           call check( nf90_put_att(ncid, varid(8+iplot), UNITS, 'flux') )
          enddo

          dimids4 = (/dimid(1),dimid(2),dimid(3),dimid(6)/)
          do iplot = 1, 3
           call check( nf90_def_var(ncid, fldname(iplot+3), NF90_REAL, dimids4, varid(11+iplot)) )
           call check( nf90_put_att(ncid, varid(11+iplot), UNITS, 'flux') )
          enddo

          dimids  = (/dimid(1),dimid(2),dimid(3),dimid(4),dimid(6)/)
          do iplot = 1, 4
           call check( nf90_def_var(ncid, fldname(iplot+6), NF90_REAL, dimids, varid(14+iplot)) )
           call check( nf90_put_att(ncid, varid(14+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(7), param) )
          call check( nf90_put_var(ncid, varid(8), array_r_fv) )
          call check( nf90_put_var(ncid, varid(3), array_vx_fv) )
          call check( nf90_put_var(ncid, varid(4), array_vz_fv) )
          call check( nf90_put_var(ncid, varid(5), array_vy_fv) )
          call check( nf90_put_var(ncid, varid(6), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(6)) )
          call check( nf90_put_var(ncid, varid(6), stime*aion, start = istart1) )
          do iplot = 1, 3
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(8+iplot)) )
          enddo
          do iplot = 1, 3
             call check( nf90_inq_varid(ncid, fldname(iplot+3), varid(11+iplot)) )
          enddo
          do iplot = 1, 4
             call check( nf90_inq_varid(ncid, fldname(iplot+6), varid(14+iplot)) )
          enddo
        endif
  endif

  allocate(fv(kinde,7,nr_fv,nx_fv,ny_fv),fvt(kinde,7,nr_fv,nx_fv,ny_fv))
  allocate(fv4(kinde,3,nr_fv,nx_fv),fv4t(kinde,3,nr_fv,nx_fv))

  fv    = 0.
  fv4   = 0.
  
  do mn=1,kinds
     if(flage) particle = eles(mn)
     if(flagi) particle = ions(mn)
    
     frac   = particle%frac
     amass  = particle%amass
     index  = ind(mn + 10)
	 if(index == 0) goto 221
     do m=1,mblocks
           do L = 1,particle%block(m)%mi
			  
			  if(pkind == 2) scale  = sqrt(Teles(mn))

              p      = particle%block(m)%qv(1,L)
              q      = particle%block(m)%qv(2,L)
              w      = particle%block(m)%qv(3,L)
			  weight = particle%block(m)%qv(14,L)
! dfactor maybe to df/f			  
			  dFactor = particle%block(m)%qv(13,L) * particle%block(m)%qv(14,L) * frac
! factor must handle full-f	e.g. to check v*E/v*B		  
			  Factor  = frac
              
			  if(ghybrid) then
			     xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
		         x     = xyz(1)
		         y     = xyz(2)
		         z     = xyz(3)
              else
			     x     = p
				 y     = q
				 z     = w
              endif
			  do n= 1, nr_fv
			     if(x< domain_vregion(n,1) .or. x>= domain_vregion(n,2) .or. & 
			        y< domain_vregion(n,3) .or. y>= domain_vregion(n,4) .or. & 
			        z< domain_vregion(n,5) .or. z>= domain_vregion(n,6)) goto 222

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1

                 vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 ve1   = ev1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             ev1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 ev1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 ev1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 ev1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             ev1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 ev1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 ev1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
				 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2) 
                 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

				 uu       = particle%block(m)%qv(4:6,L)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.

				 vv       = uu /gamma     ! scale factor here
				 uu       = uu /scale 

				 vw(1)    = x2b(1,1)*uu(1) + x2b(1,2)*uu(2)+x2b(1,3)*uu(3)
				 vw(2)    = x2b(2,1)*uu(1) + x2b(2,2)*uu(2)+x2b(2,3)*uu(3)
				 vw(3)    = x2b(3,1)*uu(1) + x2b(3,2)*uu(2)+x2b(3,3)*uu(3)

				 vpar_    = vw(1)
				 vper_    = max(sqrt(vw(2)**2+vw(3)**2),1.e-6)
				 vper1_   = vw(2)
				 vper2_   = vw(3)
				 ee       = 1.   ! differiential distribution


				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 
				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)

				 ew(1)    = x2b(1,1)*ve1(1) + x2b(1,2)*ve1(2)+x2b(1,3)*ve1(3)
				 ew(2)    = x2b(2,1)*ve1(1) + x2b(2,2)*ve1(2)+x2b(2,3)*ve1(3)
				 ew(3)    = x2b(3,1)*ve1(1) + x2b(3,2)*ve1(2)+x2b(3,3)*ve1(3)
				 abw      = max(sqrt(bw(1)**2+bw(2)**2+bw(3)**2),1.e-8_p2)
				 aew      = max(sqrt(ew(1)**2+ew(2)**2+ew(3)**2),1.e-8_p2)

				 vde1     = ew(1)*vw(1) / aew   ! /abw                    !!!
				 vde2     = (ew(2)*vw(2) + ew(3)*vw(3)) / aew   !/abw
				 vdb1     = bw(1)*vw(1) / abw    !/abw
				 vdb2     = (bw(2)*vw(2) + bw(3)*vw(3)) / abw   !/abw
				  
				 vdEB     = qelectron * (/vde1,vde2,vdb1,vdb2/)/mecell/domain_vregion(n,7)/alpha 

                 if(vpar_<vmin_fv .or. vpar_>=vmax_fv .or. &
                    vper1_<vmin_fv .or. vper1_>=vmax_fv .or. &
                    vper2_<vmin_fv .or. vper2_>=vmax_fv .or. &
                    vper_ <vmin_fv*0 .or. vper_>=vmax_fv ) goto 222

                        imc1  = min(max(int(max(vpar_-vmin_fv,0.)/dvx_fv+1.),1),nx_fv-1)
		                dmc1  = imc1-(vpar_ - vmin_fv)/dvx_fv

                        imc2  = min(max(int(max(vper_-vmin_fv*0,0.)/dvz_fv+1.),1),nz_fv-1)
		                dmc2  = imc2-(vper_ - vmin_fv*0)/dvz_fv

                        imc3  = min(max(int(max(vper1_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc3  = imc3-(vper1_ - vmin_fv)/dvy_fv

                        imc4  = min(max(int(max(vper2_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc4  = imc4-(vper2_ - vmin_fv)/dvy_fv
!                       vpar - vper ------
                        fv(index,1,n,imc1,imc2)     = fv(index,1,n,imc1,imc2)    + dmc1*dmc2*dFactor * ee !/vper_	          
                        fv(index,1,n,imc1,imc2+1)   = fv(index,1,n,imc1,imc2+1)  + dmc1*(1.-dmc2)*dFactor * ee !/vper_	          
                        fv(index,1,n,imc1+1,imc2)   = fv(index,1,n,imc1+1,imc2)  + (1.-dmc1)*dmc2*dFactor * ee !/vper_	          
                        fv(index,1,n,imc1+1,imc2+1) = fv(index,1,n,imc1+1,imc2+1)+ (1.-dmc1)*(1.-dmc2)*dFactor * ee !/vper_	          
!                       vpar - vper1 ------

                        fv(index,2,n,imc1,imc3)     = fv(index,2,n,imc1,imc3)    + dmc1*dmc3*dFactor * ee	          
                        fv(index,2,n,imc1,imc3+1)   = fv(index,2,n,imc1,imc3+1)  + dmc1*(1.-dmc3)*dFactor * ee	          
                        fv(index,2,n,imc1+1,imc3)   = fv(index,2,n,imc1+1,imc3)  + (1.-dmc1)*dmc3*dFactor * ee	          
                        fv(index,2,n,imc1+1,imc3+1) = fv(index,2,n,imc1+1,imc3+1)+ (1.-dmc1)*(1.-dmc3)*dFactor * ee	          

!                       vper1 - vper2 ------
    
                        fv(index,3,n,imc3,imc4)     = fv(index,3,n,imc3,imc4)    + dmc3*dmc4*dFactor * ee	          
                        fv(index,3,n,imc3,imc4+1)   = fv(index,3,n,imc3,imc4+1)  + dmc3*(1.-dmc4)*dFactor	 * ee          
                        fv(index,3,n,imc3+1,imc4)   = fv(index,3,n,imc3+1,imc4)  + (1.-dmc3)*dmc4*dFactor	 * ee          
                        fv(index,3,n,imc3+1,imc4+1) = fv(index,3,n,imc3+1,imc4+1)+ (1.-dmc3)*(1.-dmc4)*dFactor * ee	          



!                       VdE(vpar,vper) et al  ------
                        fv(index,4:7,n,imc1,imc2)     = fv(index,4:7,n,imc1,imc2)    + dmc1*dmc2*Factor * ee  * vdeb        
                        fv(index,4:7,n,imc1,imc2+1)   = fv(index,4:7,n,imc1,imc2+1)  + dmc1*(1.-dmc2)*Factor * ee *vdeb	          
                        fv(index,4:7,n,imc1+1,imc2)   = fv(index,4:7,n,imc1+1,imc2)  + (1.-dmc1)*dmc2*Factor * ee *vdEB	          
                        fv(index,4:7,n,imc1+1,imc2+1) = fv(index,4:7,n,imc1+1,imc2+1)+ &
						                                (1.-dmc1)*(1.-dmc2)*Factor * ee *vdEB	          



                        fv4(index,1,n,imc1)     = fv4(index,1,n,imc1)    + dmc1*Factor * ee 	          
                        fv4(index,1,n,imc1+1)   = fv4(index,1,n,imc1+1)  + (1.-dmc1)*Factor * ee 	          

                        fv4(index,2,n,imc2)     = fv4(index,2,n,imc2)    + dmc2*Factor * ee /vper_	          
                        fv4(index,2,n,imc2+1)   = fv4(index,2,n,imc2+1)  + (1.-dmc2)*Factor * ee /vper_	          

                        fv4(index,3,n,imc3)     = fv4(index,3,n,imc3)    + dmc3*Factor * ee 	          
                        fv4(index,3,n,imc3+1)   = fv4(index,3,n,imc3+1)  + (1.-dmc3)*Factor * ee 	          

222             continue
			  enddo	  			   
		  enddo
     enddo
221  continue
  enddo

  CALL MPI_ALLREDUCE(fv,fvt,size(fv),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fv4,fv4t,size(fv4),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  do iplot = 1, 3   !f(vpar,vper)
     if(mype==0) then
        allocate(plot(kinde,nr_fv,nx_fv,ny_fv))
		plot = fvt(:,iplot,:,:,:)
        call check( nf90_put_var(ncid, varid(8+iplot), plot, start = istart, &
                              count = icount) )
        deallocate(plot)
     endif
  enddo

  do iplot = 1, 3   !f(vpar)
     if(mype==0) then
        allocate(plot4(kinde,nr_fv,nx_fv))
		plot4 = fv4t(:,iplot,:,:)
        call check( nf90_put_var(ncid, varid(11+iplot), plot4, start = istart4, &
                              count = icount4) )
        deallocate(plot4)
     endif
  enddo

  do iplot = 1, 4  !Je,Jb(vpar,vper)
     if(mype==0) then
        allocate(plot(kinde,nr_fv,nx_fv,ny_fv))
		plot = fvt(:,iplot+3,:,:,:)
        call check( nf90_put_var(ncid, varid(14+iplot), plot, start = istart, &
                              count = icount) )
        deallocate(plot)
     endif
  enddo


  if(mype==0) call check( nf90_close(ncid) )


  DEALLOCATE(fv,fvt,fv4,fv4t)
  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  call allocate_bvector(ev1_,-3)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_4fv_nc1file66


! -------------------------------------------------------------------------------------------
!f(vper1,vper2,Echannel)
!
subroutine output_4fvNG_nc1file66(rstart)  
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use vector_functions
  use ionRingbeam_Case

  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer     :: bv0_,bfld,bv1_
  type(particle_type)   :: particle
  real,pointer,dimension(:,:,:,:,:,:) :: fv6,fv6t
  real,pointer,dimension(:,:,:,:,:) :: plot

  real(p2) :: vminp,vmaxp,p,q,w,x,y,z,vmaxper,rstart
  real(p2) :: dcsi,dvpar,dvper,dvper1,scale
  real(8)  :: xyz(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),vb(3),vb1(3),ve1(3),&
              phase,det,wbper1(3),wbper2(3)
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k,s
  real(p2) :: v1,v2,v3,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,ww0,ww1,angle
  real(p2) :: bw(3),vw(3),ew(3),vpar_,vper_ ,vper1_,vper2_,wper1_,wper2_,weight,gamma,dweight, &
              dmc1,dmc2,dmc3,dmc4,dmc5,ee,factor,dfactor,amass,vde1,vde2,vdb,vdEB(3),ee_,vtot,frac

  integer  :: nv,ic,m,ierror,L,kind,n,nbuff,ic0,iplot,iflag,kinde,ic0_d
  integer  :: ind(20),index,mkind,pkind,mn

  character*10 fldname(9)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(18),dimids2(2),icount2(2),varid(30), &
             dimids(6),icount(6),istart(6),icount1(1),istart1(1), &
             dimids4(4),icount4(4),istart4(4),imc1,imc2,imc3,imc4,imc5

  real  :: param(80)

  logical file_exist,flagi,flage

  save ic,ic0,ic0_d
  data ic/0/,ic0/0/,ic0_d/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(nr_fv < 1 ) return
! if(ic0_save(2) >= 50) return 
 
  call allocate_bvector(bfld,3)
  call allocate_bvector(bv0_,3)
  call allocate_bvector(bv1_,3)
  do m= 1, mblocks
     bfld(m)%vector     = bv0(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv0_)
  call bconcarB(bv1,bv1_)
  call allocate_bvector(bfld,-3)
  !------------------------
! pvdomain :  1-- number of region; 2--- istep beginning to make record ; 3: number of species to be recorded 4-9 index of 3, the max number is 6

  
  kinde     = int(pvdomain(3))

! --------------------------------------------------------------
  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed



  
  param(10)  = kinde
  flagi      =.false.
  flage      =.false.
  ind        = 0
  do m=1,kinde
     param(10+m)      =  int(pvdomain(3+m))
	 ind(int(param(10+m))) = m
	 if(param(10+m) < 10) flagi = .true.
	 if(param(10+m) > 10) flage = .true.
  enddo

  do m=1,kinde
     if(param(10+m) <10) then
	    param(20+m) = sqrt(Tions(int(param(10+m)))) * vthi
		if((vring(1)+vring(2)) * vions(int(param(10+m))) > 0.) param(20+m) = 1.
     endif
     if(param(10+m) >10) param(20+m) = sqrt(Teles(int(param(10+m)-10)))
  enddo

  fldname = (/"f(v1,v2,E)","fvpar-vpe1","fvpe1-vpe2", &
              "fvpar     ","fvper     ","fvper1    " ,&
              "vdEpar    ","vdEper    ","vdB       " /)

!  write(ct,'(I3.3)')(istep/ndiagfv - 1)/10000
!  select several time step to get the distribution data e.g.50 

! write(ct,'(I6.6)')(int((istep-1)*dt)/10)*10
  write(ct,'(I6.6)')(int((istep-1)*dt)/ndiagfv)*ndiagfv

  filename = 'fvNG'//ct//'.nc'



  inquire(file=filename,exist=file_exist)
  if((.not. file_exist) .or. abs(rstart)< 5*dt ) then
      ic0_save(8) = 0 
  endif

  ic0_save(8) = ic0_save(8) +1
  ic          = ic0_save(8)


  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, 1, 1, 1, ic /)   !(nkind,nr,vpar,vper,nE,t)
  if(mype==0 ) then
        icount  = (/kinde,nr_fv, nx_fv, ny_fv, ne_channel,1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'kind', kinde, dimid(1)) )
          call check( nf90_def_var(ncid, 'kind', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nr', nr_fv, dimid(2)) )
          call check( nf90_def_var(ncid, 'nr', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nvper1', nx_fv, dimid(3)) )
          call check( nf90_def_var(ncid, 'nvper1', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nvper2', ny_fv, dimid(4)) )
          call check( nf90_def_var(ncid, 'nvper2', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nEchannel', nE_channel, dimid(5)) )
          call check( nf90_def_var(ncid, 'nEchannel', NF90_REAL, dimid(5), varid(5)) )
          call check( nf90_put_att(ncid,  varid(5), UNITS, 'Re') )


          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(6)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(6), varid(6)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(7)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(7), varid(7)) )
          call check( nf90_put_att(ncid,  varid(7), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'ver', 3, dimid(13)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(13), varid(13)) )
          call check( nf90_put_att(ncid,  varid(13), UNITS, 'non') )


          dimids2 = (/dimid(13),dimid(2)/)
          call check( nf90_def_var(ncid, 'position', NF90_REAL, dimids2, varid(8)) )
          call check( nf90_put_att(ncid,  varid(8), UNITS, 'non') )

          dimids  = (/dimid(1),dimid(2),dimid(3),dimid(4),dimid(5),dimid(6)/)
          do iplot = 1, 2
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(8+iplot)) )
           call check( nf90_put_att(ncid, varid(8+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(7), param) )
          call check( nf90_put_var(ncid, varid(8), array_r_fv) )
          call check( nf90_put_var(ncid, varid(3), array_vx_fv) )
          call check( nf90_put_var(ncid, varid(4), array_vz_fv) )
          call check( nf90_put_var(ncid, varid(5), array_E_channel) )

          call check( nf90_put_var(ncid, varid(6), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(6)) )
          call check( nf90_put_var(ncid, varid(6), stime*aion, start = istart1) )
          do iplot = 1, 2
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(8+iplot)) )
          enddo
        endif
  endif

  allocate(fv6(kinde,2,nr_fv,nx_fv,ny_fv,ne_channel),fv6t(kinde,2,nr_fv,nx_fv,ny_fv,ne_channel))

  fv6   = 0.

  do mn=1,kinds
     if(flage) particle = eles(mn)
     if(flagi) particle = ions(mn)
     frac   = particle%frac
     amass  = particle%amass
     index  = ind(mn + 10)
	 scale  = particle%vth
	 if(index == 0) goto 221
     do m=1,mblocks
           do L = 1,particle%block(m)%mi
			  

              p      = particle%block(m)%qv(1,L)
              q      = particle%block(m)%qv(2,L)
              w      = particle%block(m)%qv(3,L)
			  weight = particle%block(m)%qv(14,L)
! dfactor maybe to df/f			  
			  dFactor = particle%block(m)%qv(13,L) * particle%block(m)%qv(14,L) * frac
! factor must handle full-f	e.g. to check v*E/v*B		  
			  Factor  = frac
              
			  if(ghybrid) then
			     xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
		         x     = xyz(1)
		         y     = xyz(2)
		         z     = xyz(3)
              else
			     x     = p
				 y     = q
				 z     = w
              endif
			  do n= 1, nr_fv
			     if(x< domain_vregion(n,1) .or. x>= domain_vregion(n,2) .or. & 
			        y< domain_vregion(n,3) .or. y>= domain_vregion(n,4) .or. & 
			        z< domain_vregion(n,5) .or. z>= domain_vregion(n,6)) goto 222

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1

                 vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  


                 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
				 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2) 
                 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

				 uu       = particle%block(m)%qv(4:6,L)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu /gamma     
				 uu       = uu /scale    ! using uu here ..........................
                 vtot     = sqrt(uu(1)**2+uu(2)**2+uu(3)**2)
		         ee_      = log10(0.5*amass*max(vtot,1.e-10_p2)**2 * real(ut0) )



				 vw(1)    = x2b(1,1)*uu(1) + x2b(1,2)*uu(2)+x2b(1,3)*uu(3)
				 vw(2)    = x2b(2,1)*uu(1) + x2b(2,2)*uu(2)+x2b(2,3)*uu(3)
				 vw(3)    = x2b(3,1)*uu(1) + x2b(3,2)*uu(2)+x2b(3,3)*uu(3)

				 vpar_    = vw(1)
				 vper1_   = vw(2)
				 vper2_   = vw(3)
				 ee       = 1.   ! differiential distribution


                 wbper2   = a_cross_b(vb,vb1)
                 wbper1   = a_cross_b(wbper2,vb)
                 wbper1   = wbper1/max(sqrt(wbper1(1)**2+wbper1(2)**2+wbper1(3)**2),1.e-8_p2) 
                 wbper2   = wbper2/max(sqrt(wbper2(1)**2+wbper2(2)**2+wbper2(3)**2),1.e-8_p2) 
                 wper1_   = wbper1(1)*uu(1) + wbper1(2)*uu(2) + wbper1(3)*uu(3)  
                 wper2_   = wbper2(1)*uu(1) + wbper2(2)*uu(2) + wbper2(3)*uu(3)  



                 if(ee_ <emin_flux .or. &
                    vper1_<vmin_fv .or. vper1_>=vmax_fv .or. &
                    vper2_<vmin_fv .or. vper2_>=vmax_fv ) goto 222

                        imc1  = min(max(int(max(vper1_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc1  = imc1-(vper1_ - vmin_fv)/dvy_fv

                        imc2  = min(max(int(max(vper2_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc2  = imc2-(vper2_ - vmin_fv)/dvy_fv
						
						imc3  = min(max(int(max(ee_-emin_flux,0.)/de_channel+1.),1),ne_channel)


                        fv6(index,1,n,imc1,imc2,imc3)     = fv6(index,1,n,imc1,imc2,imc3)    + &
						                                    dmc1*dmc2*dFactor * ee !/vper_	          
                        fv6(index,1,n,imc1,imc2+1,imc3)   = fv6(index,1,n,imc1,imc2+1,imc3)  + &
						                                    dmc1*(1.-dmc2)*dFactor * ee !/vper_	          
                        fv6(index,1,n,imc1+1,imc2,imc3)   = fv6(index,1,n,imc1+1,imc2,imc3)  + &
						                                    (1.-dmc1)*dmc2*dFactor * ee !/vper_	          
                        fv6(index,1,n,imc1+1,imc2+1,imc3) = fv6(index,1,n,imc1+1,imc2+1,imc3)+ &
						                                    (1.-dmc1)*(1.-dmc2)*dFactor * ee !/vper_	          

                        imc1  = min(max(int(max(wper1_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc1  = imc1-(wper1_ - vmin_fv)/dvy_fv

                        imc2  = min(max(int(max(wper2_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc2  = imc2-(wper2_ - vmin_fv)/dvy_fv

                        fv6(index,2,n,imc1,imc2,imc3)     = fv6(index,2,n,imc1,imc2,imc3)    + &
						                                    dmc1*dmc2*dFactor * ee !/vper_	          
                        fv6(index,2,n,imc1,imc2+1,imc3)   = fv6(index,2,n,imc1,imc2+1,imc3)  + &
						                                    dmc1*(1.-dmc2)*dFactor * ee !/vper_	          
                        fv6(index,2,n,imc1+1,imc2,imc3)   = fv6(index,2,n,imc1+1,imc2,imc3)  + &
						                                    (1.-dmc1)*dmc2*dFactor * ee !/vper_	          
                        fv6(index,2,n,imc1+1,imc2+1,imc3) = fv6(index,2,n,imc1+1,imc2+1,imc3)+ &
						                                    (1.-dmc1)*(1.-dmc2)*dFactor * ee !/vper_	          


222             continue
			  enddo	  			   
		  enddo
     enddo
221  continue
  enddo




  CALL MPI_ALLREDUCE(fv6,fv6t,size(fv6),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  do iplot = 1, 2
     if(mype==0) then
        allocate(plot(kinde,nr_fv,nx_fv,ny_fv,ne_channel))
		plot = fv6t(:,iplot,:,:,:,:)
        call check( nf90_put_var(ncid, varid(8+iplot), plot, start = istart, &
                              count = icount) )
        deallocate(plot)
     endif
  enddo



  if(mype==0) call check( nf90_close(ncid) )


  DEALLOCATE(fv6,fv6t)
  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_4fvNG_nc1file66

! -------------------------------------------------------------------------------
subroutine output_DDEA_nc1file(kind)  ! 
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  use vector_functions
  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer    :: bvt,dummy
  real,  pointer,dimension(:,:) :: flux,fluxt,flux2,flux2t,&
                                   fluxe,fluxet,fluxe2,fluxe2t,den,dent,den0,den0t
  real,  pointer,dimension(:,:) :: fluxpar,fluxpart,fluxper,fluxpert
  real,  pointer,dimension(:) :: ech,ach
  integer :: ic,ic0(4),i,iplot,m,n,iflag
  character*3 fldname(6)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  :: param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,angle,vtot,qx,qy,w1,w2,w3,w4,w11,w12,w21,w22
  real(p2)  :: amass,wx1,wy1,wz1,wx0,wy0,wz0,bc(3),e(3),b(3),bb
  real(p2)  :: uu(3),gamma,beq0,wei,charge,dangle,angle0,ee0,dee

  real(p2)  :: bbpar(3), bper1(3),bper2(3),epar,eper1,eper2,&
	           vper1,vper2,vper12,vper22,evdot

  integer   :: kind,ii,jj,kk,ip,jp,kp,nEk,nEa,mn
  integer :: j,k,L,imc1,imc2,index,IERROR
  logical    file_exist
  real*8    :: xyz(3),p,q,w

  save ic,ic0
  data ic/0/,ic0/0,0,0,0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

      subroutine bconcarB(c,d)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c,d
      end subroutine bconcarB

      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
 	  end subroutine bcovcarE

  end interface
  if(kind == 0) return


  call allocate_bvector(bvt,3)
  do m=1,mblocks
     bvt(m)%vector = bv0(m)%vector +  bv1(m)%vector
  enddo
  
  call bconcarB(bvt,bcar)
  call bcovcarE(ev1,Ecar)
  call allocate_bvector(bvt,-3)



! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(iftestParticle(1) == 0) then
     nEch    =  81
     nalpha0 =  101
     Emin    =  1.         !alog10 (E in eV)
     Emax    =  6.
  endif
  dE         =  (Emax-Emin)/real(nEch-1)
  dalpha     =  180./real(nalpha0-1)  ! assuming [0,pi]
  allocate(ech(nEch),ach(nalpha0))
  do i = 1,nEch
     ech(i) = emin + (i-1.)*de
  enddo
  do i = 1,nalpha0
     ach(i) = 0. + (i-1.)*dalpha
  enddo

  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed



 
  param      = 0.
  param(1)   = nEch
  param(2)   = nalpha0
  param(3)   = Emin
  param(4)   = Emax
  param(5)   = amin
  param(6)   = amax
  param(7)   = vvmin(kind)
  param(8)   = vvmax(kind)
  param(9)   = dE
  param(10)  = dalpha
  param(11)  = dv(kind)
  param(21)  = utime
  param(22)  = uspeed
  param(23)  = ut0
  param(24)  = cspeed

  fldname = (/"DD ","Den","EA ","XI ","EV1","EV2"/)

  write(ct,'(I3.3)')((istep-istep0)/ndiagfe - 1)/100000
  ic      = mod(((istep-istep0)/ndiagfe-1),100000) + 1   !max(mod(ic,20000),1) 
  ic0(kind)     = ic0(kind)+1

  if(kind==1) filename = 'DD_EA0_i'//ct//'.nc'
  if(kind==2) filename = 'DD_EA0_e'//ct//'.nc'

  inquire(file=filename,exist=file_exist)
  if(.not. file_exist  .or. ic ==0  .or. (ic0(kind)==1 .and. irun ==0)) then
      ic0_dd(kind) = ic - 1   
  endif
  ic = ic - ic0_dd(kind)

  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, ic /)
  if(mype==0 ) then
        icount = (/nEch, nalpha0, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )


          call check( nf90_def_dim(ncid, 'nE', nEch, dimid(3)) )
          call check( nf90_def_var(ncid, 'nE', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'eV') )

          call check( nf90_def_dim(ncid, 'nalpha', nalpha0, dimid(4)) )
          call check( nf90_def_var(ncid, 'nalpha', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, ' ') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )

          dimids=(/dimid(3),dimid(4),dimid(5)/)
          do iplot = 1, 6
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(3), ech) )
          call check( nf90_put_var(ncid, varid(4), ach) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 6
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif

  deallocate(ech,ach)
  allocate(flux(nEch,nalpha0),flux2(nEch,nalpha0),fluxt(nEch,nalpha0),&
           flux2t(nEch,nalpha0),fluxe(nEch,nalpha0),fluxe2(nEch,nalpha0),fluxet(nEch,nalpha0),&
           fluxe2t(nEch,nalpha0),den(nEch,nalpha0),dent(nEch,nalpha0),den0(nEch,nalpha0),den0t(nEch,nalpha0))
  allocate(fluxpar(nEch,nalpha0),fluxpart(nEch,nalpha0),fluxper(nEch,nalpha0),fluxpert(nEch,nalpha0))
  flux    = 0.
  flux2   = 0.
  fluxe   = 0.
  fluxe2  = 0.
  fluxpar = 0.
  fluxper = 0.
  den0    = 0.
  den     = 0.

  if(kind == 2) then 
  do m=1,mblocks
     dO  L = 1, block(m)%me
         mn       = block(m)%ele(L)%kind
	     charge   = qelectron 
	     amass    = aelectron * meles(mn)
		 wei      = block(m)%ele(L)%w * frace(mn)

         p        = block(m)%ele(L)%p(1)
         q        = block(m)%ele(L)%p(2)
         w        = block(m)%ele(L)%p(3)

	     wx1      = p +1.0 
         ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	     ip       = min(block(m)%nx1,ii+1)
	     wx1      = wx1 - ii
	     wx0      = 1.0 - wx1

	     wy1      = q +1.0 
         jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	     jp       = min(block(m)%ny1,jj+1)
	     wy1      = wy1 - jj
	     wy0      = 1.0 - wy1

         wz1      = w +1.0 
         kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	     kp       = min(block(m)%nz1,kk+1)
	     wz1      = wz1 - kk
	     wz0      = 1.0 - wz1

         e    = Ecar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			    Ecar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
			    Ecar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
			    Ecar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
			    Ecar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			    Ecar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
			    Ecar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
			    Ecar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
         
		 b    = Bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			    Bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
			    Bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
			    Bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
			    Bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			    Bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
			    Bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
			    Bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

!        b    = Bcar(m)%vector(:,ii,jj,kk)
		 bx       = b(1)
		 by       = b(2)
		 bz       = b(3)
		 BB       = sqrt(bx**2+by**2+bz**2)

         uu(1:3)  = block(m)%ele(L)%v/amass
	     gamma    = sqrt(1.+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
		 ee       = (gamma-1.)*cspeed**2*amass * uT0  
		 ee       = log10(ee)

         if(ifrelastic == 0) gamma   = 1.
         uu       = uu/gamma
		 vx       = uu(1)
		 vy       = uu(2)
		 vz       = uu(3)
		 vtot     = sqrt(uu(1)**2+uu(2)**2+uu(3)**2)
         bbpar(1)= bx/bb
		 bbpar(2)= by/bb
		 bbpar(3)= bz/bb
	     if((bb .ne. bb) .or. bb <= 1.e-4) bbpar = zaix
	     bper1  = a_cross_b(bbpar,xaix)
		 if(abs_a(bper1) <= 1.e-4_p2) bper1 = a_cross_b(bbpar,xaix)
		 bper1  = bper1/abs_a(bper1)
		 bper2  = a_cross_b(bbpar,bper1)
		 bper2  = bper2/abs_a(bper2)


         VPAR   = vx*bbpar(1)  + vy*bbpar(2)  + vz*bbpar(3)
         VPer1  = vx*bper1(1)  + vy*bper1(2)  + vz*bper1(3)
         VPer2  = vx*bper2(1)  + vy*bper2(2)  + vz*bper2(3)

         EPAR   = e(1)*bbpar(1)  + e(2)*bbpar(2)  + e(3)*bbpar(3)
         EPer1  = e(1)*bper1(1)  + e(2)*bper1(2)  + e(3)*bper1(3)
         EPer2  = e(1)*bper2(1)  + e(2)*bper2(2)  + e(3)*bper2(3)

		 evdot  = eper1 * vper1 + eper2*vper2
		 vper   = sqrt(abs(vtot**2 - vpar**2))
		 angle  = acos( max(min((vpar)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2]
! to get equatorial pitch angle for dipole case

         if(dipole) then
		    xyz           = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
            x             = xyz(1)
            y             = xyz(2)
            z             = xyz(3)
			bc            = dipole3DB(x,y,z)
		    BB            = sqrt(bc(1)**2 + bc(2)**2 + bc(3)**2)

            bbpar(1)= bc(1)/bb
		    bbpar(2)= bc(2)/bb
		    bbpar(3)= bc(3)/bb
	        if((bb .ne. bb) .or. bb <= 1.e-4) bbpar = zaix
	        bper1  = a_cross_b(bbpar,xaix)
		    if(abs_a(bper1) <= 1.e-4_p2) bper1 = a_cross_b(bbpar,xaix)
		    bper1  = bper1/abs_a(bper1)
		    bper2  = a_cross_b(bbpar,bper1)
		    bper2  = bper2/abs_a(bper2)


            VPAR   = vx*bbpar(1)  + vy*bbpar(2)  + vz*bbpar(3)
            VPer1  = vx*bper1(1)  + vy*bper1(2)  + vz*bper1(3)
            VPer2  = vx*bper2(1)  + vy*bper2(2)  + vz*bper2(3)

            EPAR   = e(1)*bbpar(1)  + e(2)*bbpar(2)  + e(3)*bbpar(3)
            EPer1  = e(1)*bper1(1)  + e(2)*bper1(2)  + e(3)*bper1(3)
            EPer2  = e(1)*bper2(1)  + e(2)*bper2(2)  + e(3)*bper2(3)

		    evdot  = eper1 * vper1 + eper2*vper2
		    vper   = sqrt(abs(vtot**2 - vpar**2))

		    xyz           = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !get magnetic field at equator
            x             = xyz(1)
            y             = xyz(2)
            z             = xyz(3)
			bc            = dipole3DB(x,y,z)
		    Beq0          = sqrt(bc(1)**2 + bc(2)**2 + bc(3)**2)

			vper          = vper * sqrt(Beq0/max(BB,1.e-6))
			vpar          = vpar/max(abs(vpar),1.e-10_p2)*sqrt(abs(vtot**2-vper**2))
!		    angle         = asin( max(min(abs(vper)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2] alpha0
		    angle         = acos( max(min(abs(vpar)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., pi/2]
            
			if(ic == 1) then
!               block(m)%ele(L)%e0 = ee        ! need to be modified here at late time
!               block(m)%ele(L)%w0 = angle
			endif
		 endif
!		 angle0 = block(m)%ele(L)%w0  
!		 ee0    = block(m)%ele(L)%e0  

         if(ee >= emin .and. ee < emax .and. mn == 1) then 

         if(p < idomain_TP_region(1,1)-2  .or. p > idomain_TP_region(1,2)+2 .or. &
            q < idomain_TP_region(1,3)-50 .or. q > idomain_TP_region(1,4)+50 .or. &
            w < idomain_TP_region(1,5) .or. w > idomain_TP_region(1,6) ) goto 777


		    dee   = 10.**ee - 10.**ee0
			dangle= angle-angle0 
		    qy    = (angle0 -0.)/dalpha  + 1.
		    qx    = (ee0 -Emin) /de + 1.
            if(int(qx) < 1 .or. int(qy) < 1 .or. int(qx) >= nEch .or. int(qy) >= nalpha0) goto 771  

		    i     = min(int(qx),nEch-1)
		    j     = min(int(qy),nalpha0-1)
		    w2    = qx - i
		    w1    = 1.- w2
			w4    = qy - j    
			w3    = 1. -w4    

            w11   = w1*w3
            w12   = w1*w4
            w21   = w2*w3
            w22   = w2*w4
            flux(i,j)        = flux(i,j)     + w11* dangle * wei
            flux(i,j+1)      = flux(i,j+1)   + w12* dangle * wei 
            flux(i+1,j)      = flux(i+1,j)   + w21* dangle * wei 
            flux(i+1,j+1)    = flux(i+1,j+1) + w22* dangle * wei 

            flux2(i,j)       = flux2(i,j)     + w11* dangle**2 * wei
            flux2(i,j+1)     = flux2(i,j+1)   + w12* dangle**2 * wei 
            flux2(i+1,j)     = flux2(i+1,j)   + w21* dangle**2 * wei
            flux2(i+1,j+1)   = flux2(i+1,j+1) + w22* dangle**2 * wei

            fluxe(i,j)        = fluxe(i,j)     + w11* dee * wei
            fluxe(i,j+1)      = fluxe(i,j+1)   + w12* dee * wei 
            fluxe(i+1,j)      = fluxe(i+1,j)   + w21* dee * wei 
            fluxe(i+1,j+1)    = fluxe(i+1,j+1) + w22* dee * wei 

            fluxe2(i,j)       = fluxe2(i,j)     + w11* dee**2 * wei
            fluxe2(i,j+1)     = fluxe2(i,j+1)   + w12* dee**2 * wei 
            fluxe2(i+1,j)     = fluxe2(i+1,j)   + w21* dee**2 * wei
            fluxe2(i+1,j+1)   = fluxe2(i+1,j+1) + w22* dee**2 * wei

            fluxpar(i,j)       = fluxpar(i,j)     + w11* vpar*epar * wei
            fluxpar(i,j+1)     = fluxpar(i,j+1)   + w12* vpar*epar * wei 
            fluxpar(i+1,j)     = fluxpar(i+1,j)   + w21* vpar*epar * wei
            fluxpar(i+1,j+1)   = fluxpar(i+1,j+1) + w22* vpar*epar * wei


            fluxper(i,j)       = fluxper(i,j)     + w11* evdot * wei
            fluxper(i,j+1)     = fluxper(i,j+1)   + w12* evdot * wei 
            fluxper(i+1,j)     = fluxper(i+1,j)   + w21* evdot * wei
            fluxper(i+1,j+1)   = fluxper(i+1,j+1) + w22* evdot * wei

            den0(i,j)         = den0(i,j)     + w11 * wei 
            den0(i,j+1)       = den0(i,j+1)   + w12 * wei 
            den0(i+1,j)       = den0(i+1,j)   + w21 * wei
            den0(i+1,j+1)     = den0(i+1,j+1) + w22 * wei 

771  continue
		    qy    = (angle -0.)/dalpha  + 1.
		    qx    = (ee -Emin) /de + 1.
            if(int(qx) < 1 .or. int(qy) < 1 .or. int(qx) >= nEch .or. int(qy) >= nalpha0) goto 772  

		    i     = min(int(qx),nEch-1)
		    j     = min(int(qy),nalpha0-1)
		    w2    = qx - i
		    w1    = 1.- w2
			w4    = qy - j    
			w3    = 1. -w4    

            w11   = w1*w3
            w12   = w1*w4
            w21   = w2*w3
            w22   = w2*w4

            den(i,j)         = den(i,j)     + w11 * wei 
            den(i,j+1)       = den(i,j+1)   + w12 * wei 
            den(i+1,j)       = den(i+1,j)   + w21 * wei
            den(i+1,j+1)     = den(i+1,j+1) + w22 * wei 
772  continue
		 endif
777      continue
	 enddo
  enddo
  endif   ! for electrons


  CALL MPI_ALLREDUCE(den,dent,size(den),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(den0,den0t,size(den0),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(flux,fluxt,size(flux),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(flux2,flux2t,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  CALL MPI_ALLREDUCE(fluxe,fluxet,size(flux),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fluxe2,fluxe2t,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  CALL MPI_ALLREDUCE(fluxpar,fluxpart,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fluxper,fluxpert,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  flux2t   = flux2t/max(den0t,1.e-4)  - (fluxt/max(den0t,1.e-4))**2 
  fluxt    = fluxt/max(den0t,1.e-4)  
!  fluxe2t  = fluxe2t/max(den0t,1.e-4) - (fluxet/max(den0t,1.e-4))**2 
  fluxe2t  = fluxet/max(den0t,1.e-4) !/ max(sqrt(abs(fluxe2t)),1.e-6_p2)

  fluxpart  = fluxpart/max(den0t,1.e-4)
  fluxpert  = fluxpert/max(den0t,1.e-4)

  do iplot = 1, 1

     if(mype==0) call check( nf90_put_var(ncid, varid(7+iplot), flux2t, start = istart, &
                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(8+iplot), dent, start = istart, &
                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(9+iplot), fluxt, start = istart, &
                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(10+iplot), fluxe2t, start = istart, &
                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(11+iplot), fluxpart, start = istart, &
                              count = icount) )
     if(mype==0) call check( nf90_put_var(ncid, varid(12+iplot), fluxpert, start = istart, &
                              count = icount) )

  enddo

  if(mype==0) call check( nf90_close(ncid) )

  deallocate(flux,fluxt,flux2t,flux2,fluxe,fluxet,fluxe2t,fluxe2,den,dent,den0,den0t, &
             fluxpar,fluxpart,fluxper,fluxpert)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_DDEA_nc1file


! ------------------------------------------------------------------------------------
subroutine output_Orbit_tps_nc1file(kind)  ! test particles
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  implicit none
  include 'mpif.h'
  type(TestParticle_type) :: TP

  real,pointer,dimension(:,:) :: send,recv,dummy,save,plot
  integer,pointer,dimension(:) :: num_recv,request
  type(bvector_type), dimension(:), pointer    :: bv0_,bv1_

  integer :: ic,ic0(4),i,iplot,m,n,iflag,mn
  character*3 fldname(2)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename

  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  :: param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,vtot,angle
  real(p2)  :: xyz(3),p,q,w,amass,u(3),weight,charge
  integer :: j,k,L,imc1,imc2,index,IERROR,kind,nsum,status(mpi_status_size),&
             iend,ind,mark,nbuff,num,ier
  logical file_exist

  logical  :: gyro
  real(p2) :: wx0,wx1,wy0,wy1,wz0,wz1,vb(3),vb1(3),vbper1(3),vbper2(3),x2b(3,3), &
              uu(3),vv(3),gamma,vw(3),vpar_,vper_,bw(3),csi
  integer  :: ii,jj,kk,ip,jp,kp


  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
     function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	 end function dipole3DB
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(kind == 0) return
 
  nsum       = mark_tps(kind) 
  ind        = 9   !ind_tps

  gyro           =  int(iftraceParticle(10)) > 0 
  if(gyro) then
     ind        = 10
     call allocate_bvector(bv0_,3)
     call allocate_bvector(bv1_,3)
     call bconcarB(bv0,bv0_)
     call bconcarB(bv1,bv1_)
  endif



  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed




  fldname = (/"PSD"," E "/)


  if(kind==1) filename = 'Orbit_tps_P.nc'
  if(kind==2) filename = 'Orbit_tps_He.nc'
  if(kind==3) filename = 'Orbit_tps_O.nc'
  if(kind==4) filename = 'Orbit_tps_e.nc'


  inquire(file=filename,exist=file_exist)
  if(.not. file_exist) then
      ic0_save(16+kind) = 0 
  endif
  ic0_save(16+kind) = ic0_save(16+kind) + 1
  ic                = ic0_save(16+kind)



  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, ic /)

  if(mype==0 ) then
        icount = (/ind, nsum, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'index', ind, dimid(2)) )
          call check( nf90_def_var(ncid, 'index', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'num', nsum, dimid(3)) )
          call check( nf90_def_var(ncid, 'num', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )



          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )


          dimids=(/dimid(2),dimid(3),dimid(5)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif

     nbuff = 2000 
     allocate(send(ind, nbuff),num_recv(numberpe),request(numberpe))
  num   = 0
  do m=1,mblocks
	 if(kind ==1) TP = block(m)%TPi
	 if(kind ==2) TP = block(m)%TPHe
     if(kind ==3) TP = block(m)%TPO
     if(kind ==4) TP = block(m)%TPe
     dO  L = 1, TP%mi
         amass    = TP%amass
		 mn       = TP%qv(9,L)    !original regions
		 mark     = TP%qv(10,L)
         
		 if(mark > 0) then
		    num  = num + 1
		    if(num >= nbuff) then
                    allocate(dummy(ind,num-1))
					dummy = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(ind,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
		    endif

			     p   = TP%qv(1,L)
			     q   = TP%qv(2,L)
			     w   = TP%qv(3,L)

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1
                 vb    = bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
				 vv      = TP%qv(4:6,L)/amass
				 vpar    = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)
				 angle   = acos(vpar/max(sqrt(vv(1)**2+vv(2)**2+vv(3)**2),1.e-6_p2))

		         xyz     = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)

			     send(1,num)   = mark
			     send(2,num)   = TP%qv(9,L)
			     send(3,num)   = xyz(1)
			     send(4,num)   = xyz(2)
			     send(5,num)   = xyz(3)
			     send(6,num)   = TP%qv(4,L)/amass
			     send(7,num)   = TP%qv(5,L)/amass
			     send(8,num)   = TP%qv(6,L)/amass
			     send(9,num)   = angle *180./pi

		 endif
	 enddo
  enddo

     call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)
     if(mype > 0 .and. num > 0) then
	       call mpi_isend(send(:,1:num),ind*num, mpi_real,0, &
		                           mype, MPI_COMM_WORLD, request(mype),ier)

     endif
     if(mype==0) then
	    allocate(save(ind,nbuff)) !nbuff
	    iplot = 0
        do i=1,num
		   iplot = iplot +1
		   if(iplot >= nbuff) then
              allocate(dummy(ind,iplot-1))
			  dummy = save(:,1:iplot-1)
			  deallocate(save)
			  nbuff = nbuff*1.5
			  allocate(save(ind,nbuff))
			  save(:,1:iplot-1) = dummy
			  deallocate(dummy)
		   endif
           save(:,iplot) = send(:,i)
        enddo
     endif

     do L=1,numberpe-1
	    if(num_recv(L+1) > 0) then
	       if(mype==0) then
		      allocate(recv(ind,num_recv(L+1)))
              call mpi_recv(recv,ind*num_recv(L+1), mpi_real,L, &
		                           L, MPI_COMM_WORLD, status, ier)
		  
		      do i=1,num_recv(L+1)
		         iplot = iplot +1
		         if(iplot >= nbuff) then
                    allocate(dummy(ind,iplot-1))
			        dummy = save(:,1:iplot-1)
			        deallocate(save)
			        nbuff = nbuff*1.5
			        allocate(save(ind,nbuff))
			        save(:,1:iplot-1) = dummy
			        deallocate(dummy)
		         endif
                 save(:,iplot) = recv(:,i)
              enddo
			  deallocate(recv)
           elseif(mype==L)then
		      call MPI_Wait(request(L), status, ier)
		   endif
		endif
	 enddo 

  if(mype==0) then
     allocate(plot(ind, nsum))
     plot = 0.
	 do L = 1,iplot
        mark = save(1,L)
		plot(:,mark) = save(:,L) 
	 enddo
  endif

  do iplot = 1, 1

     if(mype==0) call check( nf90_put_var(ncid, varid(7+iplot), plot, start = istart, &
                              count = icount) )

  enddo

  if(mype ==0) call check( nf90_close(ncid) )
  if(mype ==0) deallocate(plot)

  deallocate(send,request,num_recv)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_Orbit_tps_nc1file


! -------------------------------------------------------------------------------------------

subroutine output_Orbit_trs_nc1file(kind)  ! tracing particles
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  use trace_particles
  use vector_functions
  implicit none
  include 'mpif.h'

  real,pointer,dimension(:,:) :: send,recv,dummy,save,plot
  type(bvector_type), dimension(:), pointer    :: bv0_,bv1_
  type(particle_type)   :: particle

  integer :: ic,ic0(4),i,iplot,m,n,iflag,mn
  character*3 fldname(2)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename

  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  :: param(80),btot,vdotEpar,vdotEper
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper(3),ee,vtot,angle,ve(3),epar_,eper_(3)
  real(p2)  :: xyz(3),p,q,w,amass,u(3),weight,charge
  integer :: j,k,L,imc1,imc2,index,IERROR,kind,numt,status(mpi_status_size),&
             iend,ind,nbuff,num,num0,ier
  integer(p8) :: mark
  logical file_exist

  logical  :: gyro
  real(p2) :: wx0,wx1,wy0,wy1,wz0,wz1,vb(3),vb1(3),vbper1(3),vbper2(3),x2b(3,3), &
              uu(3),vv(3),gamma,vw(3),vpar_,vper_,bw(3),csi
  integer  :: ii,jj,kk,ip,jp,kp,mkind,pkind
  integer,pointer,dimension(:) ::  nsum,nsum0


  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
     function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	 end function dipole3DB
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(kind == 0) return

  if(mod(simtype(2),100) == 66 ) then
     call output_Orbit_trs_nc1file_gyro(kind)
	 return
  endif
 


  numt       = mark_trs(kind) 
  if(numt==0) return

  ind        = 10

  gyro           =  int(iftraceParticle(10)) > 0 
  if(gyro) then
     ind        = 11
     call allocate_bvector(bv0_,3)
     call allocate_bvector(bv1_,3)
     call bconcarB(bv0,bv0_)
     call bconcarB(bv1,bv1_)
  endif


  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = cspeed
  param(5)   = ulength
  param(6)   = uBfield
  param(7)   = uefield
  param(8)   = va
  param(9)   = alpha
  param(10)  = kinds
  param(11:kinds+10)      = mions(1:kinds)
  if(kind == 1) param(21) = aion
  if(kind == 2) param(21) = aelectron


  fldname = (/"qv "," E "/)


  if(kind==1) filename = 'Orbit_trs_ion.nc'
  if(kind==2) filename = 'Orbit_trs_ele.nc'


  inquire(file=filename,exist=file_exist)
  if(.not. file_exist) then
      ic0_save(20+kind) = 0 
  endif
  ic0_save(20+kind) = ic0_save(20+kind) + 1
  ic                = ic0_save(20+kind)



  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, ic /)

  if(mype==0 ) then
        icount = (/ind, numt, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'index', ind, dimid(2)) )
          call check( nf90_def_var(ncid, 'index', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'num', numt, dimid(3)) )
          call check( nf90_def_var(ncid, 'num', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )



          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )


          dimids=(/dimid(2),dimid(3),dimid(5)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif

  nbuff = numt 
  allocate(send(ind, nbuff))
  num   = 0

  do mn = 1,kinds
     if(kind == 1) particle = ions(mn)
     if(kind == 2) particle = eles(mn)
     amass = particle%amass
       do m=1,mblocks
          do  L = 1, particle%block(m)%mi
		      mark     = particle%block(m)%mark(3,L)    
		      if(mark > 0) then
		         num  = num + 1
		         if(num >= nbuff) then
                    allocate(dummy(ind,num-1))
					dummy = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(ind,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
		         endif
			     p   = particle%block(m)%qv(1,L)
			     q   = particle%block(m)%qv(2,L)
			     w   = particle%block(m)%qv(3,L)

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1
				 xyz      = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
                 if(dipole .and. gyro) then
				  	 vb = dipole3DB(xyz(1),xyz(2),xyz(3)) 
                     vb1= bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
				          bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
					      bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                  else
                     vb = bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						  bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						  bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  


                 endif
                 ve     = Ecar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              Ecar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
				 		  Ecar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  Ecar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  Ecar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              Ecar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						  Ecar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						  Ecar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1 
						   
		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
				 btot    = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
				 
				 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2
				 uu       = particle%block(m)%qv(4:6,L)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma
				  
				 vtot    = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
				 vpar    = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)
				 vper    = vv - vpar*vb
				 angle   = acos(vpar/max(vtot,1.e-6_p2))
                 Epar_   = vb(1)*ve(1) + vb(2)*vE(2) + vb(3)*vE(3)
				 Eper_   = vE - Epar_*vb
				 vdotEpar = Epar_*vpar
				 vdotEper = Eper_(1)*vper(1) + Eper_(2)*vper(2) + Eper_(3)*vper(3) 



			     send(1,num)   = mark
			     send(2,num)   = xyz(1)
			     send(3,num)   = xyz(2)
			     send(4,num)   = xyz(3)
			     send(5,num)   = vtot
			     send(6,num)   = angle *180./pi
			     send(7,num)   = vdotepar
			     send(8,num)   = vdoteper
			     send(9,num)   = Btot
			     send(10,num)  = mn

                 if(gyro) then
				    vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				    vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				    vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				    vpar_    = vw(1)
				    vper_    = sqrt(vw(2)**2+vw(3)**2)

				    bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				    bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				    bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)

                    csi      = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                    csi      = mod(csi + pi2, pi2)
                    send(11,num)  = csi
                 endif

		       endif
	        enddo
          enddo
     enddo

     CALL MPI_ALLREDUCE(num, num0, 1, MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
     allocate(nsum(numberpe),nsum0(numberpe))
	 allocate(save(ind,num0))


     nsum           = 0
     nsum(mype + 1) = num
     CALL MPI_ALLREDUCE(nsum, nsum0,size(nsum),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
!    call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)

     if(mype/=0 .and. num >0) then
        call mpi_send(send(:,1:num), ind*num, mpi_real, &
  	                        0, mype+100, MPI_COMM_WORLD,IERROR)
     endif

     if(mype==0) then
        save(:,1:num) = send(:,1:num)
        do m=1,numberpe-1
	       if(nsum0(m+1) > 0) then
              num   = sum(nsum0(1:m))
		      allocate(recv(ind,nsum0(m+1)))
		      call mpi_recv(recv,size(recv),mpi_real, &
			                       m, m+100, MPI_COMM_WORLD,status,IERROR)
              save(:,num+1:num+nsum0(m+1))  = recv
			  deallocate(recv)
           endif
        enddo
        
		allocate(plot(ind, num0))
	    do m = 1,num0
           mark = int(save(1,m)/10)
		   plot(:,mark) = save(:,m) ! two informations in save(1,L) index*10 + mn
	    enddo

        do iplot = 1, 1
           call check( nf90_put_var(ncid, varid(7+iplot), plot, start = istart, &
                              count = icount) )
        enddo

        call check( nf90_close(ncid) )
        deallocate(plot)
	 endif 

  deallocate(send,save)
  deallocate(nsum,nsum0)


  if(gyro) then
     call allocate_bvector(bv0_,-3)
     call allocate_bvector(bv1_,-3)
  endif

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_Orbit_trs_nc1file

! -------------------------------------------------------------------------------------------

subroutine output_Orbit_trs_nc1file_gyro(kind)  ! racing particles
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  use trace_particles
  use vector_functions
  implicit none
  include 'mpif.h'

  real,pointer,dimension(:,:) :: send,recv,dummy,save,plot,TMsav
  type(bvector_type), dimension(:), pointer    :: bv0_,bv1_
  type(particle_type)  :: particle

  integer :: ic,ic0(4),i,iplot,m,n,iflag,mn
  character*3 fldname(2)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename

  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  :: param(80),btot,vdotEpar,vdotEper
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper(3),ee,vtot,angle,ve(3),epar_,eper_(3)
  real(p2)  :: xyz(3),p,q,w,amass,u(3),weight,charge,csi_,vdotB,vdotE,wper1_,wper2_
  integer :: j,k,L,imc1,imc2,index,IERROR,kind,numt,status(mpi_status_size),&
             iend,ind,nbuff,num,num0,ier,mkind,pkind
  integer(p8) :: mark
  logical file_exist

  logical  :: gyro
  real(p2) :: wx0,wx1,wy0,wy1,wz0,wz1,vb(3),vb1(3),vbper1(3),vbper2(3),&
              wbper1(3),wbper2(3),x2b(3,3), &
              uu(3),vv(3),gamma,vw(3),vpar_,vper_,bw(3),ew(3),csi, &
			  theta_v,theta_w

  integer  :: ii,jj,kk,ip,jp,kp
  integer,pointer,dimension(:) ::  nsum,nsum0
  integer,pointer,dimension(:) ::  Shell,indi


  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
     function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	 end function dipole3DB
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(kind == 0) return


  numt       = mark_trs(kind)
  if(numt > 10000) then
     write(*,*)'wrong particle numbers' , numt
	 return
  endif
  if(numt==0) return

  ind        = 10

  gyro       =  int(iftraceParticle(10)) > 0 

  if(gyro) then
     ind        = 16
     call allocate_bvector(bv0_,3)
     call allocate_bvector(bv1_,3)
     call bconcarB(bv0,bv0_)
     call bconcarB(bv1,bv1_)
  endif


  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = cspeed
  param(5)   = ulength
  param(6)   = uBfield
  param(7)   = uefield
  param(8)   = va
  param(9)   = alpha
  param(10)  = kinds
  param(11:kinds+10)      = mions(1:kinds)
  if(kind == 1) param(21) = aion
  if(kind == 2) param(21) = aelectron


  fldname = (/"qv "," E "/)


  if(kind==1) filename = 'Orbit_trs_ion.nc'
  if(kind==2) filename = 'Orbit_trs_ele.nc'


  inquire(file=filename,exist=file_exist)
  if(.not. file_exist) then
      ic0_save(20+kind) = 0 
  endif
  ic0_save(20+kind) = ic0_save(20+kind) + 1
  ic                = ic0_save(20+kind)

  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, ic /)

  if(mype==0 ) then
        icount = (/ind, numt, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'index', ind, dimid(2)) )
          call check( nf90_def_var(ncid, 'index', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'num', numt, dimid(3)) )
          call check( nf90_def_var(ncid, 'num', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )



          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )


          dimids=(/dimid(2),dimid(3),dimid(5)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif

  nbuff = numt 
  allocate(send(ind, nbuff))
  num   = 0

  do mn=1,kinds
     if(kind == 1) particle =ions(mn)
     if(kind == 2) particle =eles(mn)

     amass = particle%amass
       do m=1,mblocks
          dO  L = 1, particle%block(m)%mi
		      mark     = particle%block(m)%mark(3,L)    
		      if(mark > 0 .and. mark < 200000 ) then
		         num  = num + 1
		         if(num >= nbuff) then
                    allocate(dummy(ind,num-1))
					dummy = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(ind,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
		         endif
			     p   = particle%block(m)%qv(1,L)
			     q   = particle%block(m)%qv(2,L)
			     w   = particle%block(m)%qv(3,L)

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1
				 xyz      = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)

                 vb1      = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
				          bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
					      bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
				  	 vb = dipole3DB(xyz(1),xyz(2),xyz(3)) 
                  else
                     vb = bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						  bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						  bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 endif
                 ve     = Ecar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              Ecar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
				 		  Ecar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  Ecar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  Ecar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              Ecar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						  Ecar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						  Ecar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1 
						   
		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
				 btot    = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
				 
				 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

                 wbper2   = a_cross_b(vb,vb1)
                 wbper1   = a_cross_b(wbper2,vb)
                 wbper1   = wbper1/max(sqrt(wbper1(1)**2+wbper1(2)**2+wbper1(3)**2),1.e-8_p2) 
                 wbper2   = wbper2/max(sqrt(wbper2(1)**2+wbper2(2)**2+wbper2(3)**2),1.e-8_p2) 

! using v ------------------------------------------------------------
				 uu       = particle%block(m)%qv(4:6,L)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma
				  
				 vtot    = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
				 vpar    = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)

				 vper    = vv - vpar*vb
				 angle   = acos(vpar/max(vtot,1.e-6_p2))

                 Epar_   = vb(1)*ve(1) + vb(2)*vE(2) + vb(3)*vE(3)
				 Eper_   = vE - Epar_*vb
				 vdotEpar = Epar_*vpar
				 vdotEper = Eper_(1)*vper(1) + Eper_(2)*vper(2) + Eper_(3)*vper(3) 

				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)


				 ew(1)    = x2b(1,1)*ve(1) + x2b(1,2)*ve(2)+x2b(1,3)*ve(3)
				 ew(2)    = x2b(2,1)*ve(1) + x2b(2,2)*ve(2)+x2b(2,3)*ve(3)
				 ew(3)    = x2b(3,1)*ve(1) + x2b(3,2)*ve(2)+x2b(3,3)*ve(3)


                 csi_    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 csi_    = mod(csi_ + pi2*2  + pi ,pi2)

				 vdotB   = vw(1)*bw(1) + vw(2)*bw(2) + vw(3)*bw(3) 
				 vdotE   = vw(1)*ew(1) + vw(2)*ew(2) + vw(3)*ew(3) 

                 wper1_   = wbper1(1)*vv(1) + wbper1(2)*vv(2) + wbper1(3)*vv(3)  
                 wper2_   = wbper2(1)*vv(1) + wbper2(2)*vv(2) + wbper2(3)*vv(3)
				   
                 theta_w  = atan2theta(wper2_,wper1_) 
                 theta_v  = atan2theta(vw(3),vw(2)) 



			     send(1,num)   = mark
			     send(2,num)   = xyz(1)
			     send(3,num)   = xyz(2)
			     send(4,num)   = xyz(3)
			     send(5,num)   = vtot
			     send(6,num)   = angle *180./pi
			     send(7,num)   = Btot
			     send(8,num)   = bw(1)      
			     send(9,num)   = sqrt(bw(2)**2+bw(3)**2)      
			     send(10,num)  = ew(1)      
			     send(11,num)  = sqrt(ew(2)**2+ew(3)**2)      
			     send(12,num)  = vdotE
			     send(13,num)  = vdotB

                 send(14,num)  = csi_
                 send(15,num)  = theta_v
                 send(16,num)  = theta_w

		       endif
	        enddo
         enddo
   enddo

     CALL MPI_ALLREDUCE(num, num0, 1, MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
     allocate(nsum(numberpe),nsum0(numberpe))
	 allocate(save(ind,num0))


     nsum           = 0
     nsum(mype + 1) = num
     CALL MPI_ALLREDUCE(nsum, nsum0,size(nsum),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
!    call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)

     if(mype/=0 .and. num >0) then
        call mpi_send(send(:,1:num), ind*num, mpi_real, &
  	                        0, mype+100, MPI_COMM_WORLD,IERROR)
     endif

     if(mype==0) then
        save(:,1:num) = send(:,1:num)
        do m=1,numberpe-1
	       if(nsum0(m+1) > 0) then
              num   = sum(nsum0(1:m))
		      allocate(recv(ind,nsum0(m+1)))
		      call mpi_recv(recv,size(recv),mpi_real, &
			                       m, m+100, MPI_COMM_WORLD,status,IERROR)
              save(:,num+1:num+nsum0(m+1))  = recv
			  deallocate(recv)
           endif
        enddo
        
		allocate(plot(ind, num0))
	    do m = 1,num0
           mark = int(save(1,m)/10)
		   plot(:,mark) = save(:,m) ! two informations in save(1,L) index*10 + mn
	    enddo

!  ----  shell apply ---------------------------
!        allocate (shell(num0),indi(num0),TMsav(ind,num0))
!		TMsav  = plot
!		shell  = plot(1,:)
!		call shellIs(num0,shell,1,indi)
!		do m=1,num0
!		   plot(:,m) = TMsav(:,indi(m))
!        enddo
!		deallocate(shell,indi,TMsav)
! -------------------------------------------
        do iplot = 1, 1
           call check( nf90_put_var(ncid, varid(7+iplot), plot, start = istart, &
                              count = icount) )
        enddo

        call check( nf90_close(ncid) )
        deallocate(plot)
	 endif 

  deallocate(send,save)
  deallocate(nsum,nsum0)


  if(gyro) then
     call allocate_bvector(bv0_,-3)
     call allocate_bvector(bv1_,-3)
  endif

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_Orbit_trs_nc1file_gyro
! -------------------------------------------------------------------------------------------


subroutine output_Orbit_Rtrs_nc1file(kind)  ! racing particles
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  use trace_particles
  use vector_functions
  implicit none
  include 'mpif.h'

  real,pointer,dimension(:,:) :: send,recv,dummy,save,plot,TMsav
  type(bvector_type), dimension(:), pointer    :: bv0_,bv1_
  type(particle_type)   :: particle

  integer :: ic,ic0(4),i,iplot,m,n,iflag,mn
  character*3 fldname(2)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename

  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  :: param(80),btot,vdotEpar,vdotEper
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper(3),ee,vtot,angle,ve(3),epar_,eper_(3)
  real(p2)  :: xyz(3),p,q,w,amass,u(3),weight,charge,csi_,vdotB,vdotE,wper1_,wper2_
  integer :: j,k,L,imc1,imc2,index,IERROR,kind,numt,status(mpi_status_size),&
             iend,ind,nbuff,num,num0,ier,pkind,mkind
  integer(p8) :: mark
  logical file_exist

  logical  :: gyro
  real(p2) :: wx0,wx1,wy0,wy1,wz0,wz1,vb(3),vb1(3),vbper1(3),vbper2(3),&
              wbper1(3),wbper2(3),x2b(3,3), &
              uu(3),vv(3),gamma,vw(3),vpar_,vper_,bw(3),ew(3),csi, &
			  theta_v,theta_w

  integer  :: ii,jj,kk,ip,jp,kp
  integer,pointer,dimension(:) ::  nsum,nsum0
  integer,pointer,dimension(:) ::  Shell,indi


  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
     function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	 end function dipole3DB
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(kind == 0) return


  numt       = mark_Rtrs(kind)
  if(iftraceParticle(9) == 13) numt  = mark_Rtrs(kind+5) 

  if(numt >= 20000) then
     write(*,*)'wrong particle numbers' , numt
	 return
  endif
  if(numt==0) return

  ind        = 10

  gyro       =  int(iftraceParticle(10)) > 0 

  if(gyro) then
     ind        = 16
     call allocate_bvector(bv0_,3)
     call allocate_bvector(bv1_,3)
     call bconcarB(bv0,bv0_)
     call bconcarB(bv1,bv1_)
  endif


  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = cspeed
  param(5)   = ulength
  param(6)   = uBfield
  param(7)   = uefield
  param(8)   = va
  param(9)   = alpha
  param(10)  = kinds
  param(11:kinds+10)      = mions(1:kinds)
  if(kind == 1) param(21) = aion
  if(kind == 2) param(21) = aelectron


  fldname = (/"qv "," E "/)


  if(kind==1) filename = 'Orbit_Rtrs_ion.nc'
  if(kind==2) filename = 'Orbit_Rtrs_ele.nc'


  inquire(file=filename,exist=file_exist)
  if(.not. file_exist) then
      ic0_save(25+kind) = 0 
  endif
  ic0_save(25+kind) = ic0_save(25+kind) + 1
  ic                = ic0_save(25+kind)

  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, ic /)

  if(mype==0 ) then
        icount = (/ind, numt, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'index', ind, dimid(2)) )
          call check( nf90_def_var(ncid, 'index', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'num', numt, dimid(3)) )
          call check( nf90_def_var(ncid, 'num', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )



          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )


          dimids=(/dimid(2),dimid(3),dimid(5)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif

  nbuff = numt 
  allocate(send(ind, nbuff))
  num   = 0


  do mn=1,kinds
     if(kind == 1) particle = ions(mn)
     if(kind == 2) particle = eles(mn)
     amass = particle%amass
       do m=1,mblocks
          do  L = 1, particle%block(m)%mi
		      mark     = particle%block(m)%mark(3,L)    
		      if(mark > 0 .and. mark < 200000 ) then
		         num  = num + 1
		         if(num >= nbuff) then
                    allocate(dummy(ind,num-1))
					dummy = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(ind,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
		         endif
			     p   = particle%block(m)%qv(1,L)
			     q   = particle%block(m)%qv(2,L)
			     w   = particle%block(m)%qv(3,L)

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1
				 xyz      = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)

                 vb1      = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
				          bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
					      bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
				  	 vb = dipole3DB(xyz(1),xyz(2),xyz(3)) 
                  else
                     vb = bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						  bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						  bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 endif
                 ve     = Ecar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              Ecar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
				 		  Ecar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  Ecar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  Ecar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              Ecar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						  Ecar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						  Ecar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1 
						   
		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
				 btot    = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
				 
				 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

                 wbper2   = a_cross_b(vb,vb1)
                 wbper1   = a_cross_b(wbper2,vb)
                 wbper1   = wbper1/max(sqrt(wbper1(1)**2+wbper1(2)**2+wbper1(3)**2),1.e-8_p2) 
                 wbper2   = wbper2/max(sqrt(wbper2(1)**2+wbper2(2)**2+wbper2(3)**2),1.e-8_p2) 

! using v ------------------------------------------------------------
				 uu       = particle%block(m)%qv(4:6,L)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma
				  
				 vtot    = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
				 vpar    = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)

				 vper    = vv - vpar*vb
				 angle   = acos(vpar/max(vtot,1.e-6_p2))

                 Epar_   = vb(1)*ve(1) + vb(2)*vE(2) + vb(3)*vE(3)
				 Eper_   = vE - Epar_*vb
				 vdotEpar = Epar_*vpar
				 vdotEper = Eper_(1)*vper(1) + Eper_(2)*vper(2) + Eper_(3)*vper(3) 

				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)


				 ew(1)    = x2b(1,1)*ve(1) + x2b(1,2)*ve(2)+x2b(1,3)*ve(3)
				 ew(2)    = x2b(2,1)*ve(1) + x2b(2,2)*ve(2)+x2b(2,3)*ve(3)
				 ew(3)    = x2b(3,1)*ve(1) + x2b(3,2)*ve(2)+x2b(3,3)*ve(3)


                 csi_    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 csi_    = mod(csi_ + pi2*2  + pi ,pi2)

				 vdotB   = vw(1)*bw(1) + vw(2)*bw(2) + vw(3)*bw(3) 
				 vdotE   = vw(1)*ew(1) + vw(2)*ew(2) + vw(3)*ew(3) 

                 wper1_   = wbper1(1)*vv(1) + wbper1(2)*vv(2) + wbper1(3)*vv(3)  
                 wper2_   = wbper2(1)*vv(1) + wbper2(2)*vv(2) + wbper2(3)*vv(3)
				   
                 theta_w  = atan2theta(wper2_,wper1_) 
                 theta_v  = atan2theta(vw(3),vw(2)) 



			     send(1,num)   = mark
			     send(2,num)   = xyz(1)
			     send(3,num)   = xyz(2)
			     send(4,num)   = xyz(3)
			     send(5,num)   = vtot
			     send(6,num)   = angle *180./pi
			     send(7,num)   = Btot
			     send(8,num)   = bw(1)      
			     send(9,num)   = sqrt(bw(2)**2+bw(3)**2)      
			     send(10,num)  = ew(1)      
			     send(11,num)  = sqrt(ew(2)**2+ew(3)**2)      
			     send(12,num)  = vdotE
			     send(13,num)  = vdotB

                 send(14,num)  = csi_
                 send(15,num)  = theta_v
                 send(16,num)  = theta_w

		       endif
	        enddo
         enddo
     enddo
     
	 CALL MPI_ALLREDUCE(num, num0, 1, MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
     allocate(nsum(numberpe),nsum0(numberpe))
	 allocate(save(ind,num0))


     nsum           = 0
     nsum(mype + 1) = num
     CALL MPI_ALLREDUCE(nsum, nsum0,size(nsum),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
!    call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)

     if(mype/=0 .and. num >0) then
        call mpi_send(send(:,1:num), ind*num, mpi_real, &
  	                        0, mype+100, MPI_COMM_WORLD,IERROR)
     endif

     if(mype==0) then
        save(:,1:num) = send(:,1:num)
        do m=1,numberpe-1
	       if(nsum0(m+1) > 0) then
              num   = sum(nsum0(1:m))
		      allocate(recv(ind,nsum0(m+1)))
		      call mpi_recv(recv,size(recv),mpi_real, &
			                       m, m+100, MPI_COMM_WORLD,status,IERROR)
              save(:,num+1:num+nsum0(m+1))  = recv
			  deallocate(recv)
           endif
        enddo
        
		allocate(plot(ind, num0))
	    do m = 1,num0
           mark = save(1,m)/10
		   plot(:,mark) = save(:,m) ! two informations in save(1,L) index*10 + mn
	    enddo

!  ----  shell apply ---------------------------
!        allocate (shell(num0),indi(num0),TMsav(ind,num0))
!		TMsav  = plot
!		shell  = plot(1,:)
!		call shellIs(num0,shell,1,indi)
!		do m=1,num0
!		   plot(:,m) = TMsav(:,indi(m))
!        enddo
!		deallocate(shell,indi,TMsav)
! -------------------------------------------
        do iplot = 1, 1
           call check( nf90_put_var(ncid, varid(7+iplot), plot, start = istart, &
                              count = icount) )
        enddo

        call check( nf90_close(ncid) )
        deallocate(plot)
	 endif 

  deallocate(send,save)
  deallocate(nsum,nsum0)


  if(gyro) then
     call allocate_bvector(bv0_,-3)
     call allocate_bvector(bv1_,-3)
  endif

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_Orbit_Rtrs_nc1file

! --------------------------------------------------------------------------------
! fv(vpar,vper),fv(vper1,vper2),fv(dvper1,dvper2),,fv(vpar,csi,vper),den(csi,q),  je,jb(q)
subroutine output_fvcsi_Trs_nc(kind,rstart)  
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use vector_functions
  use ionRingbeam_Case
  use trace_particles

  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer    :: bv0_,bfld,bv1_,ev1_
  real,pointer,dimension(:,:,:,:) :: fv1,fv2,fv3,fv4,fv1t,fv2t,fv3t,fv4t
  real,pointer,dimension(:,:,:)  :: plot1
  real,pointer,dimension(:,:,:,:) :: plot2
  real,pointer,dimension(:) :: vpar,vper,cxi

  real(p2) :: vminp,vmaxp,p,q,w,x,y,z,vmaxper,rstart
  real(p2) :: dcsi,dvpar,dvper,dvper1,scale
  real(8)  :: xyz(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),&
              vb(3),vb1(3),ve1(3),phase,det,vb_eq(3),B_eq,B_
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k,s,nregion,nre
  real(p2) :: v1,v2,v3,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,ww0,ww1,angle
  real(p2) :: bw(3),vw(3),ew(3),vpar_,vper_ ,vper1_,vper2_,weight,gamma,dweight, &
              dmc1,dmc2,dmc3,dmc4,dmc5,dmc6,dmc7,ee,factor,dfactor,amass,vpermax,&
			  abw,aew

  real(p2) :: vde_,vdb_, wbper1(3),wbper2(3),wx2b(3,3),wb2x(3,3),wper1_,wper2_,qflux(5), &
              upar_,uper1_,uper2_,uper_,upar_eq,uper_eq
  integer  :: nv,ic,m,ierror,L,n,nbuff,ic0,iplot,iflag,kinde,ic0_d,num,numt
  integer  :: ind(20),index,kind,mark,ncsi,nsum,nvper,imc1,imc2,imc3,imc4,imc5,imc6,imc7

  character*10 fldname(10)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(30),varid(30), &
             dimids1(4),icount1(4),istart1(4),icount0(1),istart0(1), &
             dimids2(5),icount2(5),istart2(5),&
			 dimids3(4),icount3(4),istart3(4),&
			 dimids4(5),icount4(5),istart4(5)


  real  :: param(80)

  logical file_exist,flagi,flage

  save ic,ic0,ic0_d
  save dimid,varid,ncid,dimids1,icount1,istart1,icount0,istart0,&
       dimids2,icount2,istart2,dimids3,icount3,istart3

  data ic/0/,ic0/0/,ic0_d/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(kind == 0) return
  nsum       = mark_trs(kind) 
  if(nsum == 0) return
  nregion    = iftraceParticle(10 + kind )
  if(nregion == 0) return
! ---------------------------------------------------
  call allocate_bvector(bfld,3)
  call allocate_bvector(bv0_,3)
  call allocate_bvector(bv1_,3)
  call allocate_bvector(ev1_,3)
  do m= 1, mblocks
     bfld(m)%vector     = bv0(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv0_)
  do m= 1, mblocks
     bfld(m)%vector     = bv1(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv1_)
  call bcovcarE(ev1, ev1_)
  call allocate_bvector(bfld,-3)
  !------------------------
! pvdomain :  1-- number of region; 2--- istep beginning to make record ; 3: number of species to be recorded 4-9 index of 3, the max number is 6

  nv         =  101  ! 101
  ncsi       =  61  ! 101
  nvper      =  6
  vmaxp      =  5.
  vminp      = -5.
  

  dvpar      = (vmaxp-vminp)/real(nv-1)
  dcsi       = (pi*2)/real(ncsi-1)
  
  dvper      = 1.
  vpermax    = 6.
  dvper      = vpermax/nvper

  allocate(vpar(nv),cxi(ncsi),vper(nvper))

  do i = 1, nv
     vpar(i) = vminp + (i-1)*dvpar
  enddo

  do i = 1, ncsi
	 cxi(i)  = (i-1) * dcsi
  enddo

  do i = 1, nvper
     vper(i) = 0. + (i-0.5)*dvper
  enddo

! --------------------------------------------------------------
  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed
  param(8)   = va
  param(9)   = alpha
  param(10)  = kinds
  param(11:kinds+10)      = mions(1:kinds)
  if(kind == 1) param(21) = aion
  if(kind == 2) param(21) = aelectron

  fldname = (/"fvpar-vper","fver1-vpe2","fvdpe1vpe2", &  ! fv(3,vpar,vper,t)
              "fvparCsi  ","nCsiQ     ","JE1csiQ   " ,&  ! fv(vpar,csi,vper,t)
              "JE2csiQ   ","JB1csiQ   ","JB2csiQ   " ,&
			  "f(si,nv,Q)"/)   !fv(3,csi,nq,t)


  write(ct,'(I6.6)')(int((istep-1)*dt)/1)*1   ! 2023???

  filename = 'TracingF'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if((.not. file_exist) .or. abs(rstart)< 5* dt ) then
      ic0_save(29) = 0 
  endif

  ic0_save(29)= ic0_save(29) +1
  ic          = ic0_save(29)

! ------------------------------------------------------------------------------

  icount0 = (/1/)
  istart0 = (/ic/)
  istart1 = (/1, 1, 1,  ic   /)   !(nkind,nr,vpar,vper,t,nregion)
  istart2 = (/1, 1, 1 , 1, ic/)   !(nkind,nr,vpar,vper,t,nregion)
  istart3 = (/1, 1, 1 , ic   /)   !(nkind,nr,vpar,vper,t,nregion)
  istart4 = (/1, 1, 1,  1, ic/)   !(nkind,nr,vpar,vper,t,nregion)

  if(mype==0 ) then
        icount1 = (/nx_fv, ny_fv,   nregion,   1 /)
        icount2 = (/ncsi,  nv,      nvper,     nregion,  1 /)
        icount3 = (/ncsi,  nalongQ, nregion,   1/)
        icount4 = (/ncsi,  nv,      nalongQ,   nregion , 1/)

!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'kind', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'kind', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nupar', nx_fv, dimid(2)) )
          call check( nf90_def_var(ncid, 'nupar', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nuper1', ny_fv, dimid(3)) )
          call check( nf90_def_var(ncid, 'nuper1', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'ncsi', ncsi, dimid(4)) )
          call check( nf90_def_var(ncid, 'ncsi', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nv', nv, dimid(5)) )
          call check( nf90_def_var(ncid, 'nv', NF90_REAL, dimid(5), varid(5)) )
          call check( nf90_put_att(ncid,  varid(5), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'CHvper', nvper, dimid(6)) )
          call check( nf90_def_var(ncid, 'CHvper', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nalongQ', nalongQ, dimid(7)) )
          call check( nf90_def_var(ncid, 'nalongQ', NF90_REAL, dimid(7), varid(7)) )
          call check( nf90_put_att(ncid,  varid(7), UNITS, 'Re') )


          call check( nf90_def_dim(ncid, 'nregion', nregion, dimid(26)) )
          call check( nf90_def_var(ncid, 'nregion', NF90_REAL, dimid(26), varid(26)) )
          call check( nf90_put_att(ncid,  varid(26), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(8)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(8), varid(8)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(9)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(9), varid(9)) )
          call check( nf90_put_att(ncid,  varid(9), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nuper', nz_fv, dimid(10)) )
          call check( nf90_def_var(ncid, 'nuper', NF90_REAL, dimid(10), varid(10)) )
          call check( nf90_put_att(ncid,  varid(10), UNITS, 'Re') )



! -----------------------------------------------

          dimids1  = (/dimid(2),dimid(3),dimid(26),dimid(8)/)
          do iplot = 1, 3
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids1, varid(10+iplot)) )
           call check( nf90_put_att(ncid, varid(10+iplot), UNITS, 'flux') )
          enddo

          dimids2 = (/dimid(4),dimid(5),dimid(6),dimid(26),dimid(8)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot+3), NF90_REAL, dimids2, varid(13+iplot)) )
           call check( nf90_put_att(ncid, varid(13+iplot), UNITS, 'flux') )
          enddo

          dimids3 = (/dimid(4),dimid(7),dimid(26),dimid(8)/)
          do iplot = 1, 5
           call check( nf90_def_var(ncid, fldname(iplot+4), NF90_REAL, dimids3, varid(14+iplot)) )
           call check( nf90_put_att(ncid, varid(14+iplot), UNITS, 'flux') )
          enddo

          dimids4 = (/dimid(4),dimid(5),dimid(7),dimid(26),dimid(8)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot+9), NF90_REAL, dimids4, varid(19+iplot)) )
           call check( nf90_put_att(ncid, varid(19+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(2), array_vx_fv) )
          call check( nf90_put_var(ncid, varid(3), array_vy_fv) )
          call check( nf90_put_var(ncid, varid(4), cxi) )
          call check( nf90_put_var(ncid, varid(5), vpar) )
          call check( nf90_put_var(ncid, varid(6), vper) )
          call check( nf90_put_var(ncid, varid(7), array_alongQ) )
          call check( nf90_put_var(ncid, varid(9), param) )
          call check( nf90_put_var(ncid, varid(10),array_vz_fv) )

          call check( nf90_put_var(ncid, varid(8), stime*aion, start = istart0) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(8)) )
          call check( nf90_put_var(ncid, varid(8), stime*aion, start = istart0) )
          
		  do iplot = 1, 3
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(10+iplot)) )
          enddo
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot+3), varid(13+iplot)) )
          enddo
          do iplot = 1, 5
             call check( nf90_inq_varid(ncid, fldname(iplot+4), varid(14+iplot)) )
          enddo
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot+9), varid(19+iplot)) )
          enddo
       
	    endif
  endif
! ---------------------------------------------------------------------------

  allocate(fv1(3,nx_fv,ny_fv,nregion) , fv1t(3,nx_fv,ny_fv,nregion))
  allocate(fv2(ncsi,nv,nvper,nregion),  fv2t(ncsi,nv,nvper,nregion))
  allocate(fv3(5,ncsi,nalongQ,nregion), fv3t(5,ncsi,nalongQ,nregion))
  allocate(fv4(ncsi,nv,nalongQ,nregion),fv4t(ncsi,nv,nalongQ,nregion))

  fv1   = 0.
  fv2   = 0.
  fv3   = 0.
  fv4   = 0.
  num   = 0

  if(kind == 2) then
     do m=1,mblocks
           do L = 1,block(m)%me
		      mark     = block(m)%ele(L)%mark(3)    
		      if(mark > 0) then
			     amass = meles(block(m)%ele(L)%kind)*aelectron
				 nre    = mod(mark,10)
			     p   = block(m)%ele(L)%p(1)
			     q   = block(m)%ele(L)%p(2)
			     w   = block(m)%ele(L)%p(3)
			     weight = block(m)%ele(L)%w * block(m)%ele(L)%ww  * frace(block(m)%ele(L)%kind)

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1

                 vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 ve1   = ev1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             ev1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 ev1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 ev1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 ev1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             ev1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 ev1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 ev1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
					 B_   = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
		             xyz  = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb_eq= dipole3DB(x,y,z)
					 B_eq = sqrt(vb_eq(1)**2+vb_eq(2)**2+vb_eq(3)**2)
				 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2) 
                 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

				 uu       = block(m)%ele(L)%v(:)/aelectron
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.

				 vv       = uu /gamma     ! scale factor here
!				 uu       = uu /scale 
! convert to equatorial distribution by keeping miu=cont. ---------------------------
				 vw(1)    = x2b(1,1)*uu(1) + x2b(1,2)*uu(2)+x2b(1,3)*uu(3)
				 vw(2)    = x2b(2,1)*uu(1) + x2b(2,2)*uu(2)+x2b(2,3)*uu(3)
				 vw(3)    = x2b(3,1)*uu(1) + x2b(3,2)*uu(2)+x2b(3,3)*uu(3)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)
				 
				 ew(1)    = x2b(1,1)*ve1(1) + x2b(1,2)*ve1(2)+x2b(1,3)*ve1(3)
				 ew(2)    = x2b(2,1)*ve1(1) + x2b(2,2)*ve1(2)+x2b(2,3)*ve1(3)
				 ew(3)    = x2b(3,1)*ve1(1) + x2b(3,2)*ve1(2)+x2b(3,3)*ve1(3)
				 abw      = max(sqrt(bw(1)**2+bw(2)**2+bw(3)**2),1.e-8_p2)
				 aew      = max(sqrt(ew(1)**2+ew(2)**2+ew(3)**2),1.e-8_p2)

                 angle    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 angle    = mod(angle+pi2*2 + pi*0,pi2)    ! southward ????

! for diagnosising f(upar,uper) --------------------------
				 upar_    = vw(1)
				 uper_    = max(sqrt(vw(2)**2+vw(3)**2),1.e-8)

                 uper_eq  = uper_ !uper_*sqrt(B_eq/B_)
				 upar_eq  = upar_ !upar_/max(abs(upar_),1.e-8_p2) &
				           !*sqrt(max(upar_**2+uper_**2-uper_eq**2,1.e-8_p2))


                 uper1_   = vbper1(1)*uu(1) + vbper1(2)*uu(2) + vbper1(3)*uu(3)  
                 uper2_   = vbper2(1)*uu(1) + vbper2(2)*uu(2) + vbper2(3)*uu(3)  

				 ee       = 1.   ! differiential distribution
				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

! ----------------------------------------------------------------------------------


                 vde_     = a_dot_b(vv, ve1) /max(sqrt(ve1(1)**2+ve1(2)**2+ve1(3)**2),1.e-8_p2)
                 vdb_     = a_dot_b(vv, vb1) /max(sqrt(vb1(1)**2+vb1(2)**2+vb1(3)**2),1.e-8_p2)
				 qflux    = (/1._p2, qelectron*vw(1)*ew(1)/aew /abw,  &
				                     qelectron*(vw(2)*ew(2)+vw(3)*ew(3))/aew /abw, &
							         qelectron*vw(1)*bw(1)/abw/abw,  &
				                     qelectron*(vw(2)*bw(2)+vw(3)*bw(3))/abw/abw  &
				 
				            /)

                 wbper2   = a_cross_b(vb,vb1)
                 wbper1   = a_cross_b(wbper2,vb)
                 wbper1   = wbper1/max(sqrt(wbper1(1)**2+wbper1(2)**2+wbper1(3)**2),1.e-8_p2) 
                 wbper2   = wbper2/max(sqrt(wbper2(1)**2+wbper2(2)**2+wbper2(3)**2),1.e-8_p2) 
                 wper1_   = wbper1(1)*uu(1) + wbper1(2)*uu(2) + wbper1(3)*uu(3)  
                 wper2_   = wbper2(1)*uu(1) + wbper2(2)*uu(2) + wbper2(3)*uu(3)  



!                 if(vpar_ <vmin_fv .or. vpar_ >=vmax_fv .or. &
!                    vper1_<vmin_fv .or. vper1_>=vmax_fv .or. &
!                    vper2_<vmin_fv .or. vper2_>=vmax_fv .or. &
!                    vper_ <vmin_fv*0 .or. vper_>=vmax_fv ) goto 221

                        imc1  = min(max(int(max(upar_eq-vmin_fv,0.)/dvx_fv+1.),1),nx_fv-1)
		                dmc1  = imc1-(upar_eq - vmin_fv)/dvx_fv

                        imc2  = min(max(int(max(uper_eq-vmin_fv*0,0.)/dvz_fv+1.),1),nz_fv-1)
		                dmc2  = imc2-(uper_eq - vmin_fv*0)/dvz_fv

                        imc3  = min(max(int(max(uper1_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc3  = imc3-(uper1_ - vmin_fv)/dvy_fv

                        imc4  = min(max(int(max(uper2_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc4  = imc4-(uper2_ - vmin_fv)/dvy_fv

                        imc5  = min(max(int(max(wper1_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc5  = imc5-(wper1_ - vmin_fv)/dvy_fv

                        imc6  = min(max(int(max(wper2_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc6  = imc6-(wper2_ - vmin_fv)/dvy_fv

!                       vpar - vper ------
                        fv1(1,imc1,imc2,nre)     = fv1(1,imc1,imc2,nre)    + dmc1*dmc2*weight * ee !/vper_	          
                        fv1(1,imc1,imc2+1,nre)   = fv1(1,imc1,imc2+1,nre)  + dmc1*(1.-dmc2)*weight * ee !/vper_	          
                        fv1(1,imc1+1,imc2,nre)   = fv1(1,imc1+1,imc2,nre)  + (1.-dmc1)*dmc2*weight * ee !/vper_	          
                        fv1(1,imc1+1,imc2+1,nre) = fv1(1,imc1+1,imc2+1,nre)+ (1.-dmc1)*(1.-dmc2)*weight * ee !/vper_	          
!                       vper1 - vper2 ------

                        fv1(2,imc3,imc4,nre)     = fv1(2,imc3,imc4,nre)    + dmc3*dmc4*weight * ee	          
                        fv1(2,imc3,imc4+1,nre)   = fv1(2,imc3,imc4+1,nre)  + dmc3*(1.-dmc4)*weight	 * ee          
                        fv1(2,imc3+1,imc4,nre)   = fv1(2,imc3+1,imc4,nre)  + (1.-dmc3)*dmc4*weight	 * ee          
                        fv1(2,imc3+1,imc4+1,nre) = fv1(2,imc3+1,imc4+1,nre)+ (1.-dmc3)*(1.-dmc4)*weight * ee	          
!                       vper1 - vper2 ------based on wave frame ---------------------------------------------

                        fv1(3,imc5,imc6,nre)     = fv1(3,imc5,imc6,nre)    + dmc5*dmc6*weight * ee	          
                        fv1(3,imc5,imc6+1,nre)   = fv1(3,imc5,imc6+1,nre)  + dmc5*(1.-dmc6)*weight	 * ee          
                        fv1(3,imc5+1,imc6,nre)   = fv1(3,imc5+1,imc6,nre)  + (1.-dmc5)*dmc6*weight	 * ee          
                        fv1(3,imc5+1,imc6+1,nre) = fv1(3,imc5+1,imc6+1,nre)+ (1.-dmc5)*(1.-dmc6)*weight * ee	          
!   ----------------------------------------------------------------------------------------------------------
                        imc1  = min(int(angle/dcsi + 1.),ncsi-1)
		                dmc1  = imc1 - angle/dcsi

                        imc2  = max(min(int((vpar_ - vminp)/dvpar + 1.),  nv-1),1)
		                dmc2  = imc2-(vpar_ - vminp)/dvpar

                        imc3  = min(int(vper_/dvper +1.),nvper)

                        fv2(imc1,imc2,imc3,nre)     = fv2(imc1,imc2,imc3,nre)    + dmc1*dmc2*weight * ee !/vper_	          
                        fv2(imc1,imc2+1,imc3,nre)   = fv2(imc1,imc2+1,imc3,nre)  + dmc1*(1.-dmc2)*weight * ee !/vper_	          
                        fv2(imc1+1,imc2,imc3,nre)   = fv2(imc1+1,imc2,imc3,nre)  + (1.-dmc1)*dmc2*weight * ee !/vper_	          
                        fv2(imc1+1,imc2+1,imc3,nre) = fv2(imc1+1,imc2+1,imc3,nre)+ (1.-dmc1)*(1.-dmc2)*weight * ee !/vper_	          

!   ----------------------------------------------------------------------------------------------------------
                        imc4  = min(int(q/dalongQ + 1.),nalongQ-1)
		                dmc4  = imc4 - q/dalongQ
                        fv3(:,imc1,imc4,nre)     = fv3(:,imc1,imc4,nre)    + dmc1*dmc4*weight * ee *qflux	          
                        fv3(:,imc1,imc4+1,nre)   = fv3(:,imc1,imc4+1,nre)  + dmc1*(1.-dmc4)*weight * ee  *qflux          
                        fv3(:,imc1+1,imc4,nre)   = fv3(:,imc1+1,imc4,nre)  + (1.-dmc1)*dmc4*weight * ee  *qflux          
                        fv3(:,imc1+1,imc4+1,nre) = fv3(:,imc1+1,imc4+1,nre)+ (1.-dmc1)*(1.-dmc4)*weight * ee  *qflux	          
!   ----------------------------------------------------------------------------------------------------------
                        fv4(imc1,imc2,imc4,nre)     = fv4(imc1,imc2,imc4,nre)    + dmc1*dmc2*weight * ee * dmc4	          
                        fv4(imc1,imc2+1,imc4,nre)   = fv4(imc1,imc2+1,imc4,nre)  + dmc1*(1.-dmc2)*weight * ee * dmc4	          
                        fv4(imc1+1,imc2,imc4,nre)   = fv4(imc1+1,imc2,imc4,nre)  + (1.-dmc1)*dmc2*weight * ee * dmc4	          
                        fv4(imc1+1,imc2+1,imc4,nre) = fv4(imc1+1,imc2+1,imc4,nre)+ &
						                              (1.-dmc1)*(1.-dmc2)*weight * ee * dmc4	          

                        fv4(imc1,imc2,imc4+1,nre)     = fv4(imc1,imc2,imc4+1,nre)    + dmc1*dmc2*weight * ee * (1.-dmc4)          
                        fv4(imc1,imc2+1,imc4+1,nre)   = fv4(imc1,imc2+1,imc4+1,nre)  + &
						                                dmc1*(1.-dmc2)*weight * ee * (1.-dmc4)	          
                        fv4(imc1+1,imc2,imc4+1,nre)   = fv4(imc1+1,imc2,imc4+1,nre)  + &
						                                (1.-dmc1)*dmc2*weight * ee * (1.-dmc4)	          
                        fv4(imc1+1,imc2+1,imc4+1,nre) = fv4(imc1+1,imc2+1,imc4+1,nre)+ &
						                                (1.-dmc1)*(1.-dmc2)*weight * ee * (1.-dmc4)	          


                        num   = num + 1
!222              continue
			  endif	  			   
221           continue
		  enddo
     enddo
  endif



  CALL MPI_ALLREDUCE(fv1,fv1t,size(fv1),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fv2,fv2t,size(fv2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fv3,fv3t,size(fv3),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fv4,fv4t,size(fv4),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(num,numt,1,MPI_integer,MPI_SUM,mpi_comm_world,IERROR)


  do iplot = 1, 3
     if(mype==0) then
        allocate(plot1(nx_fv,ny_fv,nregion))
		plot1 = fv1t(iplot,:,:,:)
        call check( nf90_put_var(ncid, varid(10+iplot), plot1, start = istart1, &
                              count = icount1) )
        deallocate(plot1)
     endif
  enddo

  do iplot = 1, 1
     if(mype==0) then
        allocate(plot2(ncsi,nv,nvper,nregion))
		plot2 = fv2t(:,:,:,:)
        call check( nf90_put_var(ncid, varid(13+iplot), plot2, start = istart2, &
                              count = icount2) )
        deallocate(plot2)
     endif
  enddo

  do iplot = 1, 5
     if(mype==0) then
        allocate(plot1(ncsi,nalongQ,nregion))
		plot1 = fv3t(iplot,:,:,:)
        call check( nf90_put_var(ncid, varid(14+iplot), plot1, start = istart3, &
                              count = icount3) )
        deallocate(plot1)
     endif
  enddo

  do iplot = 1, 1
     if(mype==0) then
        allocate(plot2(ncsi,nv,nalongQ,nregion))
		plot2 = fv4t(:,:,:,:)
        call check( nf90_put_var(ncid, varid(19+iplot), plot2, start = istart4, &
                              count = icount4) )
        deallocate(plot2)
     endif
  enddo


  if(mype==0) call check( nf90_close(ncid) )


  DEALLOCATE(fv1,fv1t,fv2,fv2t,fv3,fv3t,fv4,fv4t)
  DEALLOCATE(vpar,cxi,vper)

  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  call allocate_bvector(ev1_,-3)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_fvcsi_Trs_nc


! --------------------------------------------------------------------------------
! fv(vpar,vper),fv(vper1,vper2),fv(dvper1,dvper2),,fv(vpar,csi,vper),den(csi,q),  je,jb(q)
subroutine output_fvcsi_ReTrs_nc(kind,rstart)  
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use vector_functions
  use ionRingbeam_Case
  use trace_particles

  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer    :: bv0_,bfld,bv1_,ev1_
  real,pointer,dimension(:,:,:,:) :: fv1,fv2,fv3,fv4,fv1t,fv2t,fv3t,fv4t
  real,pointer,dimension(:,:,:)  :: plot1
  real,pointer,dimension(:,:,:,:) :: plot2
  real,pointer,dimension(:) :: vpar,vper,cxi

  real(p2) :: vminp,vmaxp,p,q,w,x,y,z,vmaxper,rstart
  real(p2) :: dcsi,dvpar,dvper,dvper1,scale
  real(8)  :: xyz(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),&
              vb(3),vb1(3),ve1(3),phase,det,vb_eq(3),B_eq,B_
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k,s,nregion,nre
  real(p2) :: v1,v2,v3,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,ww0,ww1,angle
  real(p2) :: bw(3),vw(3),ew(3),vpar_,vper_ ,vper1_,vper2_,weight,gamma,dweight, &
              dmc1,dmc2,dmc3,dmc4,dmc5,dmc6,dmc7,ee,factor,dfactor,amass,vpermax,&
			  abw,aew

  real(p2) :: vde_,vdb_, wbper1(3),wbper2(3),wx2b(3,3),wb2x(3,3),wper1_,wper2_,qflux(5), &
              upar_,uper1_,uper2_,uper_,upar_eq,uper_eq
  integer  :: nv,ic,m,ierror,L,n,nbuff,ic0,iplot,iflag,kinde,ic0_d,num,numt
  integer  :: ind(20),index,kind,mark,ncsi,nsum,nvper,imc1,imc2,imc3,imc4,imc5,imc6,imc7

  character*10 fldname(10)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(30),varid(30), &
             dimids1(4),icount1(4),istart1(4),icount0(1),istart0(1), &
             dimids2(5),icount2(5),istart2(5),&
			 dimids3(4),icount3(4),istart3(4),&
			 dimids4(5),icount4(5),istart4(5)
  real  :: param(80)

  logical file_exist,flagi,flage

  save ic,ic0,ic0_d
  save dimid,varid,ncid,dimids1,icount1,istart1,icount0,istart0,&
       dimids2,icount2,istart2,dimids3,icount3,istart3

  data ic/0/,ic0/0/,ic0_d/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
  if(kind == 0) return
  nsum       = mark_Rtrs(kind) 
  if(nsum == 0) return
  nregion    = iftraceParticle(10 + kind )
  if(nregion == 0) return
! ---------------------------------------------------
  call allocate_bvector(bfld,3)
  call allocate_bvector(bv0_,3)
  call allocate_bvector(bv1_,3)
  call allocate_bvector(ev1_,3)
  do m= 1, mblocks
     bfld(m)%vector     = bv0(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv0_)
  do m= 1, mblocks
     bfld(m)%vector     = bv1(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv1_)
  call bcovcarE(ev1, ev1_)
  call allocate_bvector(bfld,-3)
  !------------------------
! pvdomain :  1-- number of region; 2--- istep beginning to make record ; 3: number of species to be recorded 4-9 index of 3, the max number is 6

  nv         =  101   ! 101
  ncsi       =  61    ! 101
  nvper      =  6
  vmaxp      =  5.
  vminp      = -5.
  

  dvpar      = (vmaxp-vminp)/real(nv-1)
  dcsi       = (pi*2)/real(ncsi-1)
  
  dvper      = 1.
  vpermax    = 6.
  dvper      = vpermax/nvper

  allocate(vpar(nv),cxi(ncsi),vper(nvper))

  do i = 1, nv
     vpar(i) = vminp + (i-1)*dvpar
  enddo

  do i = 1, ncsi
	 cxi(i)  = (i-1) * dcsi
  enddo

  do i = 1, nvper
     vper(i) = 0. + (i-0.5)*dvper
  enddo

! --------------------------------------------------------------
  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed
  param(8)   = va
  param(9)   = alpha
  param(10)  = kinds
  param(11:kinds+10)      = mions(1:kinds)
  if(kind == 1) param(21) = aion
  if(kind == 2) param(21) = aelectron

  fldname = (/"fvpar-vper","fver1-vpe2","fvdpe1vpe2", &  ! fv(3,vpar,vper,t)
              "fvparCsi  ","nCsiQ     ","JE1csiQ   " ,&  ! fv(vpar,csi,vper,t)
              "JE2csiQ   ","JB1csiQ   ","JB2csiQ   " ,&
			  "f(si,nv,Q)"/)   !fv(3,csi,nq,t)


! write(ct,'(I6.6)')(int((istep-1)*dt)/10)*10
  write(ct,'(I6.6)')(int((istep-1)*dt)/1)*1

  filename = 'ReTracingF'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if((.not. file_exist) .or. abs(rstart)< 5* dt ) then
      ic0_save(30+kind) = 0 
  endif

  ic0_save(30+kind)= ic0_save(30+kind) +1
  ic               = ic0_save(30+kind)

! ------------------------------------------------------------------------------

  icount0 = (/1/)
  istart0 = (/ic/)
  istart1 = (/1, 1, 1,  ic   /)   !(nkind,nr,vpar,vper,t,nregion)
  istart2 = (/1, 1, 1 , 1, ic/)   !(nkind,nr,vpar,vper,t,nregion)
  istart3 = (/1, 1, 1 , ic   /)   !(nkind,nr,vpar,vper,t,nregion)
  istart4 = (/1, 1, 1,  1, ic/)   !(nkind,nr,vpar,vper,t,nregion)
  if(mype==0 ) then
        icount1 = (/nx_fv, ny_fv,   nregion,   1 /)
        icount2 = (/ncsi,  nv,      nvper,     nregion,  1 /)
        icount3 = (/ncsi,  nalongQ, nregion,   1/)
        icount4 = (/ncsi,  nv,      nalongQ,   nregion , 1/)

!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'kind', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'kind', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nupar', nx_fv, dimid(2)) )
          call check( nf90_def_var(ncid, 'nupar', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nuper1', ny_fv, dimid(3)) )
          call check( nf90_def_var(ncid, 'nuper1', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'ncsi', ncsi, dimid(4)) )
          call check( nf90_def_var(ncid, 'ncsi', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nv', nv, dimid(5)) )
          call check( nf90_def_var(ncid, 'nv', NF90_REAL, dimid(5), varid(5)) )
          call check( nf90_put_att(ncid,  varid(5), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'CHvper', nvper, dimid(6)) )
          call check( nf90_def_var(ncid, 'CHvper', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nalongQ', nalongQ, dimid(7)) )
          call check( nf90_def_var(ncid, 'nalongQ', NF90_REAL, dimid(7), varid(7)) )
          call check( nf90_put_att(ncid,  varid(7), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'nregion', nregion, dimid(26)) )
          call check( nf90_def_var(ncid, 'nregion', NF90_REAL, dimid(26), varid(26)) )
          call check( nf90_put_att(ncid,  varid(26), UNITS, 'Re') )



          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(8)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(8), varid(8)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(9)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(9), varid(9)) )
          call check( nf90_put_att(ncid,  varid(9), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nuper', nz_fv, dimid(10)) )
          call check( nf90_def_var(ncid, 'nuper', NF90_REAL, dimid(10), varid(10)) )
          call check( nf90_put_att(ncid,  varid(10), UNITS, 'Re') )



! -----------------------------------------------

          dimids1  = (/dimid(2),dimid(3),dimid(26),dimid(8)/)
          do iplot = 1, 3
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids1, varid(10+iplot)) )
           call check( nf90_put_att(ncid, varid(10+iplot), UNITS, 'flux') )
          enddo

          dimids2 = (/dimid(4),dimid(5),dimid(6),dimid(26),dimid(8)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot+3), NF90_REAL, dimids2, varid(13+iplot)) )
           call check( nf90_put_att(ncid, varid(13+iplot), UNITS, 'flux') )
          enddo

          dimids3 = (/dimid(4),dimid(7),dimid(26),dimid(8)/)
          do iplot = 1, 5
           call check( nf90_def_var(ncid, fldname(iplot+4), NF90_REAL, dimids3, varid(14+iplot)) )
           call check( nf90_put_att(ncid, varid(14+iplot), UNITS, 'flux') )
          enddo

          dimids4 = (/dimid(4),dimid(5),dimid(7),dimid(26),dimid(8)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot+9), NF90_REAL, dimids4, varid(19+iplot)) )
           call check( nf90_put_att(ncid, varid(19+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(2), array_vx_fv) )
          call check( nf90_put_var(ncid, varid(3), array_vy_fv) )
          call check( nf90_put_var(ncid, varid(4), cxi) )
          call check( nf90_put_var(ncid, varid(5), vpar) )
          call check( nf90_put_var(ncid, varid(6), vper) )
          call check( nf90_put_var(ncid, varid(7), array_alongQ) )
          call check( nf90_put_var(ncid, varid(9), param) )
          call check( nf90_put_var(ncid, varid(10),array_vz_fv) )

          call check( nf90_put_var(ncid, varid(8), stime*aion, start = istart0) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(8)) )
          call check( nf90_put_var(ncid, varid(8), stime*aion, start = istart0) )
          
		  do iplot = 1, 3
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(10+iplot)) )
          enddo
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot+3), varid(13+iplot)) )
          enddo
          do iplot = 1, 5
             call check( nf90_inq_varid(ncid, fldname(iplot+4), varid(14+iplot)) )
          enddo
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot+9), varid(19+iplot)) )
          enddo
       
	    endif
  endif
! ---------------------------------------------------------------------------

  allocate(fv1(3,nx_fv,ny_fv,nregion) , fv1t(3,nx_fv,ny_fv,nregion))
  allocate(fv2(ncsi,nv,nvper,nregion),  fv2t(ncsi,nv,nvper,nregion))
  allocate(fv3(5,ncsi,nalongQ,nregion), fv3t(5,ncsi,nalongQ,nregion))
  allocate(fv4(ncsi,nv,nalongQ,nregion),fv4t(ncsi,nv,nalongQ,nregion))

  fv1   = 0.
  fv2   = 0.
  fv3   = 0.
  fv4   = 0.
  num   = 0

  if(kind == 2) then
     do m=1,mblocks
           do L = 1,block(m)%me
		      mark     = block(m)%ele(L)%mark(3)    
		      if(mark > 0) then
			     amass = meles(block(m)%ele(L)%kind)*aelectron
				 scale  = sqrt(Teles(kind))
				 nre    = mod(mark,10)

			     p   = block(m)%ele(L)%p(1)
			     q   = block(m)%ele(L)%p(2)
			     w   = block(m)%ele(L)%p(3)
			     weight = block(m)%ele(L)%w * block(m)%ele(L)%ww  * frace(block(m)%ele(L)%kind)

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1

                 vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 ve1   = ev1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             ev1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 ev1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 ev1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 ev1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             ev1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 ev1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 ev1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
					 B_   = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
		             xyz  = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb_eq= dipole3DB(x,y,z)
					 B_eq = sqrt(vb_eq(1)**2+vb_eq(2)**2+vb_eq(3)**2)
				 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2) 
                 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

				 uu       = block(m)%ele(L)%v(:)/aelectron
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.

				 vv       = uu /gamma     ! scale factor here
				 uu       = uu /scale 
! convert to equatorial distribution by keeping miu=cont. ---------------------------
				 vw(1)    = x2b(1,1)*uu(1) + x2b(1,2)*uu(2)+x2b(1,3)*uu(3)
				 vw(2)    = x2b(2,1)*uu(1) + x2b(2,2)*uu(2)+x2b(2,3)*uu(3)
				 vw(3)    = x2b(3,1)*uu(1) + x2b(3,2)*uu(2)+x2b(3,3)*uu(3)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)
				 
				 ew(1)    = x2b(1,1)*ve1(1) + x2b(1,2)*ve1(2)+x2b(1,3)*ve1(3)
				 ew(2)    = x2b(2,1)*ve1(1) + x2b(2,2)*ve1(2)+x2b(2,3)*ve1(3)
				 ew(3)    = x2b(3,1)*ve1(1) + x2b(3,2)*ve1(2)+x2b(3,3)*ve1(3)
				 abw      = max(sqrt(bw(1)**2+bw(2)**2+bw(3)**2),1.e-8_p2)
				 aew      = max(sqrt(ew(1)**2+ew(2)**2+ew(3)**2),1.e-8_p2)

                 angle    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 angle    = mod(angle+pi2*2 + pi*0,pi2)    ! southward ????

! for diagnosising f(upar,uper) --------------------------
				 upar_    = vw(1)
				 uper_    = max(sqrt(vw(2)**2+vw(3)**2),1.e-8)

                 uper_eq  = uper_ !uper_*sqrt(B_eq/B_)
				 upar_eq  = upar_ !upar_/max(abs(upar_),1.e-8_p2) &
				           !*sqrt(max(upar_**2+uper_**2-uper_eq**2,1.e-8_p2))


                 uper1_   = vbper1(1)*uu(1) + vbper1(2)*uu(2) + vbper1(3)*uu(3)  
                 uper2_   = vbper2(1)*uu(1) + vbper2(2)*uu(2) + vbper2(3)*uu(3)  

				 ee       = 1.   ! differiential distribution
				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

! ----------------------------------------------------------------------------------


                 vde_     = a_dot_b(vv, ve1) /max(sqrt(ve1(1)**2+ve1(2)**2+ve1(3)**2),1.e-8_p2)
                 vdb_     = a_dot_b(vv, vb1) /max(sqrt(vb1(1)**2+vb1(2)**2+vb1(3)**2),1.e-8_p2)
				 qflux    = (/1._p2, qelectron*vw(1)*ew(1)/aew /abw,  &
				                     qelectron*(vw(2)*ew(2)+vw(3)*ew(3))/aew/abw, &
							         qelectron*vw(1)*bw(1)/abw /abw,  &
				                     qelectron*(vw(2)*bw(2)+vw(3)*bw(3))/abw/abw  &
				 
				            /)

                 wbper2   = a_cross_b(vb,vb1)
                 wbper1   = a_cross_b(wbper2,vb)
                 wbper1   = wbper1/max(sqrt(wbper1(1)**2+wbper1(2)**2+wbper1(3)**2),1.e-8_p2) 
                 wbper2   = wbper2/max(sqrt(wbper2(1)**2+wbper2(2)**2+wbper2(3)**2),1.e-8_p2) 
                 wper1_   = wbper1(1)*uu(1) + wbper1(2)*uu(2) + wbper1(3)*uu(3)  
                 wper2_   = wbper2(1)*uu(1) + wbper2(2)*uu(2) + wbper2(3)*uu(3)  



!                 if(vpar_ <vmin_fv .or. vpar_ >=vmax_fv .or. &
!                    vper1_<vmin_fv .or. vper1_>=vmax_fv .or. &
!                    vper2_<vmin_fv .or. vper2_>=vmax_fv .or. &
!                    vper_ <vmin_fv*0 .or. vper_>=vmax_fv ) goto 221

                        imc1  = min(max(int(max(upar_eq-vmin_fv,0.)/dvx_fv+1.),1),nx_fv-1)
		                dmc1  = imc1-(upar_eq - vmin_fv)/dvx_fv

                        imc2  = min(max(int(max(uper_eq-vmin_fv*0,0.)/dvz_fv+1.),1),nz_fv-1)
		                dmc2  = imc2-(uper_eq - vmin_fv*0)/dvz_fv

                        imc3  = min(max(int(max(uper1_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc3  = imc3-(uper1_ - vmin_fv)/dvy_fv

                        imc4  = min(max(int(max(uper2_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc4  = imc4-(uper2_ - vmin_fv)/dvy_fv

                        imc5  = min(max(int(max(wper1_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc5  = imc5-(wper1_ - vmin_fv)/dvy_fv

                        imc6  = min(max(int(max(wper2_-vmin_fv,0.)/dvy_fv+1.),1),ny_fv-1)
		                dmc6  = imc6-(wper2_ - vmin_fv)/dvy_fv

!                       vpar - vper ------
                        fv1(1,imc1,imc2,nre)     = fv1(1,imc1,imc2,nre)    + dmc1*dmc2*weight * ee !/vper_	          
                        fv1(1,imc1,imc2+1,nre)   = fv1(1,imc1,imc2+1,nre)  + dmc1*(1.-dmc2)*weight * ee !/vper_	          
                        fv1(1,imc1+1,imc2,nre)   = fv1(1,imc1+1,imc2,nre)  + (1.-dmc1)*dmc2*weight * ee !/vper_	          
                        fv1(1,imc1+1,imc2+1,nre) = fv1(1,imc1+1,imc2+1,nre)+ (1.-dmc1)*(1.-dmc2)*weight * ee !/vper_	          
!                       vper1 - vper2 ------

                        fv1(2,imc3,imc4,nre)     = fv1(2,imc3,imc4,nre)    + dmc3*dmc4*weight * ee	          
                        fv1(2,imc3,imc4+1,nre)   = fv1(2,imc3,imc4+1,nre)  + dmc3*(1.-dmc4)*weight	 * ee          
                        fv1(2,imc3+1,imc4,nre)   = fv1(2,imc3+1,imc4,nre)  + (1.-dmc3)*dmc4*weight	 * ee          
                        fv1(2,imc3+1,imc4+1,nre) = fv1(2,imc3+1,imc4+1,nre)+ (1.-dmc3)*(1.-dmc4)*weight * ee	          
!                       vper1 - vper2 ------based on wave frame ---------------------------------------------

                        fv1(3,imc5,imc6,nre)     = fv1(3,imc5,imc6,nre)    + dmc5*dmc6*weight * ee	          
                        fv1(3,imc5,imc6+1,nre)   = fv1(3,imc5,imc6+1,nre)  + dmc5*(1.-dmc6)*weight	 * ee          
                        fv1(3,imc5+1,imc6,nre)   = fv1(3,imc5+1,imc6,nre)  + (1.-dmc5)*dmc6*weight	 * ee          
                        fv1(3,imc5+1,imc6+1,nre) = fv1(3,imc5+1,imc6+1,nre)+ (1.-dmc5)*(1.-dmc6)*weight * ee	          
!   ----------------------------------------------------------------------------------------------------------
                        imc1  = min(int(angle/dcsi + 1.),ncsi-1)
		                dmc1  = imc1 - angle/dcsi

                        imc2  = max(min(int((vpar_ - vminp)/dvpar + 1.),  nv-1),1)
		                dmc2  = imc2-(vpar_ - vminp)/dvpar

                        imc3  = min(int(vper_/dvper +1.),nvper)

                        fv2(imc1,imc2,imc3,nre)     = fv2(imc1,imc2,imc3,nre)    + dmc1*dmc2*weight * ee !/vper_	          
                        fv2(imc1,imc2+1,imc3,nre)   = fv2(imc1,imc2+1,imc3,nre)  + dmc1*(1.-dmc2)*weight * ee !/vper_	          
                        fv2(imc1+1,imc2,imc3,nre)   = fv2(imc1+1,imc2,imc3,nre)  + (1.-dmc1)*dmc2*weight * ee !/vper_	          
                        fv2(imc1+1,imc2+1,imc3,nre) = fv2(imc1+1,imc2+1,imc3,nre)+ (1.-dmc1)*(1.-dmc2)*weight * ee !/vper_	          

!   ----------------------------------------------------------------------------------------------------------
                        imc4  = min(int(q/dalongQ + 1.),nalongQ-1)
		                dmc4  = imc4 - q/dalongQ
                        fv3(:,imc1,imc4,nre)     = fv3(:,imc1,imc4,nre)    + dmc1*dmc4*weight * ee *qflux	          
                        fv3(:,imc1,imc4+1,nre)   = fv3(:,imc1,imc4+1,nre)  + dmc1*(1.-dmc4)*weight * ee  *qflux          
                        fv3(:,imc1+1,imc4,nre)   = fv3(:,imc1+1,imc4,nre)  + (1.-dmc1)*dmc4*weight * ee  *qflux          
                        fv3(:,imc1+1,imc4+1,nre) = fv3(:,imc1+1,imc4+1,nre)+ (1.-dmc1)*(1.-dmc4)*weight * ee  *qflux	          
!   ----------------------------------------------------------------------------------------------------------
                        fv4(imc1,imc2,imc4,nre)     = fv4(imc1,imc2,imc4,nre)    + dmc1*dmc2*weight * ee * dmc4	          
                        fv4(imc1,imc2+1,imc4,nre)   = fv4(imc1,imc2+1,imc4,nre)  + dmc1*(1.-dmc2)*weight * ee * dmc4	          
                        fv4(imc1+1,imc2,imc4,nre)   = fv4(imc1+1,imc2,imc4,nre)  + (1.-dmc1)*dmc2*weight * ee * dmc4	          
                        fv4(imc1+1,imc2+1,imc4,nre) = fv4(imc1+1,imc2+1,imc4,nre)+ &
						                              (1.-dmc1)*(1.-dmc2)*weight * ee * dmc4	          

                        fv4(imc1,imc2,imc4+1,nre)     = fv4(imc1,imc2,imc4+1,nre)    + dmc1*dmc2*weight * ee * (1.-dmc4)          
                        fv4(imc1,imc2+1,imc4+1,nre)   = fv4(imc1,imc2+1,imc4+1,nre)  + &
						                                dmc1*(1.-dmc2)*weight * ee * (1.-dmc4)	          
                        fv4(imc1+1,imc2,imc4+1,nre)   = fv4(imc1+1,imc2,imc4+1,nre)  + &
						                                (1.-dmc1)*dmc2*weight * ee * (1.-dmc4)	          
                        fv4(imc1+1,imc2+1,imc4+1,nre) = fv4(imc1+1,imc2+1,imc4+1,nre)+ &
						                                (1.-dmc1)*(1.-dmc2)*weight * ee * (1.-dmc4)	          

                        num   = num + 1
!222              continue
			  endif	  			   
221           continue
		  enddo
     enddo
  endif



  CALL MPI_ALLREDUCE(fv1,fv1t,size(fv1),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fv2,fv2t,size(fv2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fv3,fv3t,size(fv3),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fv4,fv4t,size(fv4),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(num,numt,1,MPI_integer,MPI_SUM,mpi_comm_world,IERROR)


  do iplot = 1, 3
     if(mype==0) then
        allocate(plot1(nx_fv,ny_fv,nregion))
		plot1 = fv1t(iplot,:,:,:)
        call check( nf90_put_var(ncid, varid(10+iplot), plot1, start = istart1, &
                              count = icount1) )
        deallocate(plot1)
     endif
  enddo

  do iplot = 1, 1
     if(mype==0) then
        allocate(plot2(ncsi,nv,nvper,nregion))
		plot2 = fv2t(:,:,:,:)
        call check( nf90_put_var(ncid, varid(13+iplot), plot2, start = istart2, &
                              count = icount2) )
        deallocate(plot2)
     endif
  enddo

  do iplot = 1, 5
     if(mype==0) then
        allocate(plot1(ncsi,nalongQ,nregion))
		plot1 = fv3t(iplot,:,:,:)
        call check( nf90_put_var(ncid, varid(14+iplot), plot1, start = istart3, &
                              count = icount3) )
        deallocate(plot1)
     endif
  enddo

  do iplot = 1, 1
     if(mype==0) then
        allocate(plot2(ncsi,nv,nalongQ,nregion))
		plot2 = fv4t(:,:,:,:)
        call check( nf90_put_var(ncid, varid(19+iplot), plot2, start = istart4, &
                              count = icount4) )
        deallocate(plot2)
     endif
  enddo



  if(mype==0) call check( nf90_close(ncid) )


  DEALLOCATE(fv1,fv1t,fv2,fv2t,fv3,fv3t,fv4,fv4t)
  DEALLOCATE(vpar,cxi,vper)

  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  call allocate_bvector(ev1_,-3)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_fvcsi_ReTrs_nc
! -----------------------------------------------------------------------
! -------------------------------------------------------------------------------
subroutine output_DDEA_testparticles1(rstart)  ! 
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  use vector_functions
  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer    :: bvt,dummy,bv1_
  real,  pointer,dimension(:,:,:,:) :: flux2,flux2t
  real,  pointer,dimension(:,:,:) :: plot,den
  real,  pointer,dimension(:) :: ech,ach
  integer :: ic,ic0(4),i,iplot,m,n,iflag
  character*4 fldname(13)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(20),dimids2(2),icount2(2),varid(30), &
             dimids(4),icount(4),istart(4),icount1(1),istart1(1)
  real  :: param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,dmc3,dmc4,&
            bx,by,bz,vx,vy,vz,vpar,vper(3),ee,angle,vtot,qx,qy,w1,w2,w3,w4,w11,w12,w21,w22
  real(p2)  :: amass,wx1,wy1,wz1,wx0,wy0,wz0,bc(3),e(3),b(3),bb,bb0
  real(p2)  :: uu(3),gamma,beq0,wei,charge,dangle,angle0,ee0,dee

  real(p2)  :: bbpar(3), bper1(3),bper2(3),epar,eper1,eper2,&
	           vper1,vper2,vper12,vper22,evdot,rstart,vb(3),vb1(3),vb0(3),flux(13),&
			   vbper1(3),vbper2(3),wbper1(3),wbper2(3),vv(3),x2b(3,3),eper_(3),ve(3), &
			   vdotEper,btot,csi_,epar_,theta_w,vdotEpar,vpar_,vper_,wper1,wper2,&
			   ew(3),bw(3),vw(3),wper1_,wper2_,vdotB,baix(3),pLL,pp,pp0

  integer   :: kind,ii,jj,kk,ip,jp,kp,nEk,nEa,mn
  integer   :: j,k,L,imc1,imc2,imc3,imc4,index,IERROR
  logical      file_exist
  real*8    :: xyz(3),p,q,w,xyz0(3),xyz1(3),pqw(3)

  save ic,ic0
  data ic/0/,ic0/0,0,0,0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

      subroutine bconcarB(c,d)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c,d
      end subroutine bconcarB

      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
 	  end subroutine bcovcarE

  end interface

  ndiagf_tps     = max(int(iftestParticle(191)+0.1),1)
 

  call allocate_bvector(bvt,3)
  call allocate_bvector(bv1_,3)
  do m=1,mblocks
     bvt(m)%vector = bv0(m)%vector +  bv1(m)%vector
  enddo
  
  call bconcarB(bvt,bcar)
  call bcovcarE(ev1,Ecar)
  call bconcarB(bv1,bv1_)
  call allocate_bvector(bvt,-3)

! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange

  allocate(ech(ne_tps),ach(na_tps))
  do i = 1,ne_tps
     ech(i) = Emin_tps + (i-1.)*de
  enddo
  do i = 1,na_tps
     ach(i) = amin_tps + (i-1.)*dalpha
  enddo

  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed

  param      = 0.
  param(1)   = ne_tps
  param(2)   = na_tps
  param(3)   = emin_tps
  param(4)   = emax_tps
  param(5)   = amin_tps
  param(6)   = amax_tps
  param(7)   = vvmin(1)
  param(8)   = vvmax(1)
  param(9)   = dE
  param(10)  = dalpha
  param(11)  = nr_tps
  param(21)  = utime
  param(22)  = uspeed
  param(23)  = ut0
  param(24)  = cspeed
  param(25)  = ulength


  fldname = (/"DE  ","DA  ","DL  ","DDE ","DDA ","DDL ", &
              "VdE1","VdE2","VdB ","Csi ","Den0","Den1","Flux"/)



  write(ct,'(I6.6)')(int((istep-1)*dt)/ndiagf_tps) * ndiagf_tps
  filename = 'DDEA_teste1_'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if((.not. file_exist) .or. abs(rstart)< 5* dt ) then
      ic0_save(40) = 0 
  endif

  ic0_save(40) = ic0_save(40) +1
  ic           = ic0_save(40)


  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, 1,ic /)   !(nr,ne,na,t)
  if(mype==0 ) then
        icount = (/nr_tps,ne_tps,na_tps, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'nr', nr_tps, dimid(1)) )
          call check( nf90_def_var(ncid, 'nr', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )


          call check( nf90_def_dim(ncid, 'nE', ne_tps, dimid(3)) )
          call check( nf90_def_var(ncid, 'nE', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'eV') )

          call check( nf90_def_dim(ncid, 'nalpha', na_tps, dimid(4)) )
          call check( nf90_def_var(ncid, 'nalpha', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, ' ') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )

          dimids=(/dimid(1),dimid(3),dimid(4),dimid(5)/)
          do iplot = 1, 13
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(3), ech) )
          call check( nf90_put_var(ncid, varid(4), ach) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 13
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif

  deallocate(ech,ach)
  allocate(flux2(13,nr_tps,ne_tps,na_tps),flux2t(13,nr_tps,ne_tps,na_tps))
  flux2   = 0.


  do mn = 1,nr_tps
     do m  = 1,mblocks
        dO  L = 1, teste(mn)%block(m)%mi
	        charge   = teste(mn)%charge 
	        amass    = teste(mn)%amass
		    wei      = teste(mn)%block(m)%qv(7,L)

            p        = teste(mn)%block(m)%qv(1,L)
            q        = teste(mn)%block(m)%qv(2,L)
            w        = teste(mn)%block(m)%qv(3,L)

	        wx1      = p +1.0 
            ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	        ip       = min(block(m)%nx1,ii+1)
	        wx1      = wx1 - ii
	        wx0      = 1.0 - wx1

	        wy1      = q +1.0 
            jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	        jp       = min(block(m)%ny1,jj+1)
	        wy1      = wy1 - jj
	        wy0      = 1.0 - wy1

            wz1      = w +1.0 
            kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	        kp       = min(block(m)%nz1,kk+1)
	        wz1      = wz1 - kk
	        wz0      = 1.0 - wz1

            ve       = Ecar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			    Ecar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
			    Ecar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
			    Ecar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
			    Ecar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			    Ecar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
			    Ecar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
			    Ecar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
         
            vb1      = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
				          bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
					      bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
            if(dipole) then
		        xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !
                x    = xyz(1)
                y    = xyz(2)
                z    = xyz(3)
			    vb   = dipole3DB(x,y,z)
				bb   = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
		        xyz0 = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !
                x    = xyz0(1)
                y    = xyz0(2)
                z    = xyz0(3)
			    vb0  = dipole3DB(x,y,z)
				bb0  = sqrt(vb0(1)**2+vb0(2)**2+vb0(3)**2)
		    endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
				 btot    = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
				 uu      = teste(mn)%block(m)%qv(4:6,L)/amass
! --------- check D_LL ----------------				 
	             xyz1    = a_cross_b(uu,vb)/bb/charge*amass
	             x       = xyz(1)+xyz1(1)
	             y       = xyz(2)+xyz1(2)
	             z       = xyz(3)+xyz1(3)
	             pqw     = xyz2pqw(0,x,y,z)

	             pLL     = (pqw(1)-pmin)/deltax


! ---------------------------------
            if(dipole) then
		        xyz  = pqw2xyz(0,pmin+pLL*deltax,qmin+q*deltay,wmin+w*deltaz) !
                x    = xyz(1)
                y    = xyz(2)
                z    = xyz(3)
			    vb   = dipole3DB(x,y,z)
				bb   = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
		        xyz0 = pqw2xyz(0,pmin+pLL*deltax,qequator,wmin+w*deltaz) !
                x    = xyz0(1)
                y    = xyz0(2)
                z    = xyz0(3)
			    vb0  = dipole3DB(x,y,z)
				bb0  = sqrt(vb0(1)**2+vb0(2)**2+vb0(3)**2)
		    endif
		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
				 btot    = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)

! ---------------------------------
				 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

                 wbper2   = a_cross_b(vb,vb1)
                 wbper1   = a_cross_b(wbper2,vb)
                 wbper1   = wbper1/max(sqrt(wbper1(1)**2+wbper1(2)**2+wbper1(3)**2),1.e-8_p2) 
                 wbper2   = wbper2/max(sqrt(wbper2(1)**2+wbper2(2)**2+wbper2(3)**2),1.e-8_p2) 
! --------------------------------------------
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma

		         ee       = (gamma-1.)*cspeed**2*amass * uT0  
				 pp       = amass * sqrt(uu(1)**2+uu(2)**2+uu(3)**2)
		         ee       = log10(ee)
		         ee0      = teste(mn)%block(m)%qv(8,L)
				 gamma    = 10.**(ee0)/cspeed**2/amass / uT0 +1.
				 pp0      = sqrt(gamma**2-1.)*cspeed *amass

				 
				  
!				 vtot     = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
!				 vpar     = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)
		         
!				 vper_    = sqrt(vtot**2 -vpar**2)*sqrt(BB0/max(BB,1.e-10_p2))
!				 vpar_    = vpar/abs(max(vpar,1.e-10_p2))* sqrt(abs(vtot**2-vper_**2))
!				 angle    = min(max(vpar/max(vtot,1.e-10_p2),-1.),1.)

				 vtot     = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
				 vpar     = vv(1)*vb(1) + vv(2)*vb(2) + vv(3)*vb(3) 
				 vper     = vv - vpar*vb
				 vper_    = sqrt(max(vtot**2-vpar**2,0._p2)) * sqrt(BB0/max(BB,1.e-10))
			     vpar     = vpar/max(abs(vpar),1.e-10_p2)*sqrt(abs(vtot**2-vper_**2))
				 angle    = vpar/max(vtot,1.e-12)
				 if(ifrad) angle = acos(angle)*180./pi

!				 angle    = acos( max(min((vpar)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., 90]
! check equatorial pitch angle ------------
                 

                 angle0   = teste(mn)%block(m)%qv(9,L)   ![-1,1]

                 Epar_    = vb(1)*ve(1) + vb(2)*vE(2) + vb(3)*vE(3)
				 Eper_    = vE - Epar_*vb

				 vdotEpar = Epar_*vpar
				 vdotEper = Eper_(1)*vper(1) + Eper_(2)*vper(2) + Eper_(3)*vper(3) 
				 vdotB    = vb1(1)*vv(1) + vb1(2)*vv(2) + vb1(3)*vv(3) 

				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)


				 ew(1)    = x2b(1,1)*ve(1) + x2b(1,2)*ve(2)+x2b(1,3)*ve(3)
				 ew(2)    = x2b(2,1)*ve(1) + x2b(2,2)*ve(2)+x2b(2,3)*ve(3)
				 ew(3)    = x2b(3,1)*ve(1) + x2b(3,2)*ve(2)+x2b(3,3)*ve(3)


                 csi_    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 csi_    = mod(csi_ + pi2*2,  pi2)

                 wper1_   = wbper1(1)*vv(1) + wbper1(2)*vv(2) + wbper1(3)*vv(3)  
                 wper2_   = wbper2(1)*vv(1) + wbper2(2)*vv(2) + wbper2(3)*vv(3)
				   
                 theta_w  = atan2theta(wper2_,wper1_) 

                 imc1     = min(max(int(max(ee0-emin_tps,0.)/de +1.),1),ne_tps-1)
		         dmc1     = imc1-(ee0 - emin_tps)/de
                 imc2     = min(max(int(max(angle0-amin_tps,0.)/dalpha +1.),1),na_tps-1)
		         dmc2     = imc2-(angle0 - amin_tps)/dalpha

				 flux(1)  = pp - pp0   !10.**ee - 10.**ee0   !4 DDp

				 flux(2)  = abs(angle) - abs(angle0)         !4 DDa

                 if(ifrad) flux(2) = acos(abs(cos(angle*pi/180.)))*180./pi - &
				                     acos(abs(cos(angle0*pi/180.)))*180./pi
				 flux(3)  = (pLL - teste(mn)%block(m)%qv(10,L))*deltax/ulength
				 flux(4)  = flux(1)**2
				 flux(5)  = flux(2)**2
				 flux(6)  = flux(3)**2
				 flux(7)  = vdotEpar
				 flux(8)  = vdotEper
				 flux(9)  = vdotB
				 flux(10) = csi_
				 flux(11) = 1.
				 flux(12) = 1.
				 flux(13) = wei
                 imc3     = min(max(int(max(ee-emin_tps,0.)/de +1.),1),ne_tps-1)
		         dmc3     = imc3-(ee - emin_tps)/de
                 imc4     = min(max(int(max(angle-amin_tps,0.)/dalpha +1.),1),na_tps-1)
		         dmc4     = imc4-(angle - amin_tps)/dalpha

                 flux2(1:11,mn,imc1,imc2)     = flux2(1:11,mn,imc1,imc2)    + dmc1*dmc2*flux(1:11)	          
                 flux2(1:11,mn,imc1,imc2+1)   = flux2(1:11,mn,imc1,imc2+1)  + dmc1*(1.-dmc2)*flux(1:11)		          
                 flux2(1:11,mn,imc1+1,imc2)   = flux2(1:11,mn,imc1+1,imc2)  + (1.-dmc1)*dmc2*flux(1:11)		          
                 flux2(1:11,mn,imc1+1,imc2+1) = flux2(1:11,mn,imc1+1,imc2+1)+ (1.-dmc1)*(1.-dmc2)*flux(1:11)	          

                 flux2(12:13,mn,imc3,imc4)     = flux2(12:13,mn,imc3,imc4)    + dmc3*dmc4*flux(12:13)	          
                 flux2(12:13,mn,imc3,imc4+1)   = flux2(12:13,mn,imc3,imc4+1)  + dmc3*(1.-dmc4)*flux(12:13)		          
                 flux2(12:13,mn,imc3+1,imc4)   = flux2(12:13,mn,imc3+1,imc4)  + (1.-dmc3)*dmc4*flux(12:13)		          
                 flux2(12:13,mn,imc3+1,imc4+1) = flux2(12:13,mn,imc3+1,imc4+1)+ (1.-dmc3)*(1.-dmc4)*flux(12:13)	          


		 enddo
777      continue
	 enddo
  enddo


  CALL MPI_ALLREDUCE(flux2,flux2t,size(flux2),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  allocate(den(nr_tps,ne_tps,na_tps))
  allocate(plot(nr_tps,ne_tps,na_tps))
  flux2t(:,:,:,na_tps-1) = (flux2t(:,:,:,1) + flux2t(:,:,:,na_tps-1))/2. 
  flux2t(:,:,:,na_tps)   = flux2t(:,:,:,2) 
  flux2t(:,:,:,1)        = flux2t(:,:,:,na_tps-1) 

  den                  = max(flux2t(11,:,:,:),1.e-12)  !/iftestParticle(6)


  do iplot = 1, 13
     if(iplot == 1)  plot = flux2t(1,:,:,:)/den
     if(iplot == 2)  plot = flux2t(2,:,:,:)/den
     if(iplot == 3)  plot = flux2t(3,:,:,:)/den
     if(iplot == 4)  plot = flux2t(4,:,:,:)/den - (flux2t(1,:,:,:)/den)**2  
     if(iplot == 5)  plot = flux2t(5,:,:,:)/den - (flux2t(2,:,:,:)/den)**2
     if(iplot == 6)  plot = flux2t(6,:,:,:)/den - (flux2t(3,:,:,:)/den)**2
     if(iplot == 7)  plot = flux2t(7,:,:,:)/den
     if(iplot == 8)  plot = flux2t(8,:,:,:)/den
     if(iplot == 9)  plot = flux2t(9,:,:,:)/den
     if(iplot == 10) plot = flux2t(10,:,:,:)/den
     if(iplot == 11) plot = flux2t(11,:,:,:)/iftestParticle(6)
     if(iplot == 12) plot = flux2t(12,:,:,:)/iftestParticle(6)
     if(iplot == 13) plot = flux2t(13,:,:,:)/iftestParticle(6)


     if(mype==0) call check( nf90_put_var(ncid, varid(7+iplot), plot, start = istart, &
                              count = icount) )
  enddo

  if(mype==0) call check( nf90_close(ncid) )

  deallocate(flux2t,flux2,den,plot)
  call allocate_bvector(bv1_,-3)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_DDEA_testparticles1

! -------------------------------------------------------------------------
subroutine output_Orbit_testparticle1  ! racing particles
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use test_particles
  use trace_particles
  use vector_functions
  implicit none
  include 'mpif.h'

  real,pointer,dimension(:,:) :: send,recv,dummy,save,plot,TMsav
  type(bvector_type), dimension(:), pointer    :: bv0_,bv1_,ev1_,temp

  integer :: ic,ic0(4),i,iplot,m,n,iflag,mn
  character*3 fldname(2)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename

  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real(p2)  :: param(80),btot,vdotEpar,vdotEper
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper(3),ee,vtot,angle,ve(3),epar_,eper_(3), &
			vb0(3),bb,bb0
  real(p2)  :: xyz(3),p,q,w,amass,u(3),weight,charge,csi_,vdotB,vdotE,wper1_,wper2_
  integer :: j,k,L,imc1,imc2,index,IERROR,kind,numt,status(mpi_status_size),&
             iend,ind,nbuff,num,num0,ier
  integer(p8) :: mark
  logical file_exist

  logical  :: gyro
  real(p2) :: wx0,wx1,wy0,wy1,wz0,wz1,vb(3),vb1(3),vbper1(3),vbper2(3),&
              wbper1(3),wbper2(3),x2b(3,3), &
              uu(3),vv(3),gamma,vw(3),vpar_,vper_,bw(3),ew(3),csi, &
			  theta_v,theta_w

  integer  :: ii,jj,kk,ip,jp,kp
  integer,pointer,dimension(:) ::  nsum,nsum0
  integer,pointer,dimension(:) ::  Shell,indi


  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE
     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
     function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	 end function dipole3DB
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange


  numt       = mark_trs(8)
  if(numt > 10000) then
     write(*,*)'wrong particle numbers' , numt
	 return
  endif
  if(numt==0) return

     ind        = 16
     call allocate_bvector(bv0_,3)
     call allocate_bvector(bv1_,3)
     call allocate_bvector(ev1_,3)
     call allocate_bvector(temp,3)
     call bconcarB(bv0,bv0_)
     call bconcarB(bv1,bv1_)
     do m=1,mblocks
	    temp(m)%vector = ev1(m)%vector
		if(nx <=2 .and. nz <=2 ) temp(m)%vector(2,:,:,:) = 0.
     enddo
     call bcovcarE(temp,ev1_)
     call allocate_bvector(temp,-3)



  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = cspeed
  param(5)   = ulength
  param(6)   = uBfield
  param(7)   = uefield
  param(8)   = va
  param(9)   = alpha
  param(10)  = kinds
  param(11:kinds+10)      = mions(1:kinds)
  if(kind == 1) param(21) = aion
  if(kind == 2) param(21) = aelectron


  fldname = (/"qv "," E "/)


  filename = 'Orbit_testparticles1.nc'

  inquire(file=filename,exist=file_exist)
  if(.not. file_exist) then
      ic0_save(41) = 0 
  endif
  ic0_save(41) = ic0_save(41) + 1
  ic           = ic0_save(41)

  icount1 = (/1/)
  istart1 = (/ic/)
  istart  = (/1, 1, ic /)

  if(mype==0 ) then
        icount = (/ind, numt, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'ver', 3, dimid(1)) )
          call check( nf90_def_var(ncid, 'ver', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'index', ind, dimid(2)) )
          call check( nf90_def_var(ncid, 'index', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'Re') )

          call check( nf90_def_dim(ncid, 'num', numt, dimid(3)) )
          call check( nf90_def_var(ncid, 'num', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'Re') )



          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )


          dimids=(/dimid(2),dimid(3),dimid(5)/)
          do iplot = 1, 1
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids, varid(7+iplot)) )
           call check( nf90_put_att(ncid, varid(7+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart1) )
          do iplot = 1, 1
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(7+iplot)) )
          enddo
        endif
  endif

  nbuff = numt 
  allocate(send(ind, nbuff))
  num   = 0

  do mn = 1,nr_tps
       do m=1,mblocks
          dO  L = 1, teste(mn)%block(m)%mi
		      mark     = teste(mn)%block(m)%mark(3,L)    
		      if(mark > 0 .and. mark < 200000 ) then
		         num  = num + 1
		         if(num >= nbuff) then
                    allocate(dummy(ind,num-1))
					dummy = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(ind,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
		         endif
			     amass = teste(mn)%amass
			     p   = teste(mn)%block(m)%qv(1,L)
			     q   = teste(mn)%block(m)%qv(2,L)
			     w   = teste(mn)%block(m)%qv(3,L)

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1
				 xyz      = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)

                 vb1      = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
				          bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
					      bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 ve       = ev1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              ev1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  ev1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  ev1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  ev1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              ev1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
				          ev1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
					      ev1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
				     xyz      = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz)
				  	 vb0      = dipole3DB(xyz(1),xyz(2),xyz(3)) 
					 bb0      = sqrt(vb0(1)**2+vb0(2)**2+vb0(3)**2)

				     xyz      = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
				  	 vb       = dipole3DB(xyz(1),xyz(2),xyz(3)) 
					 bb       = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)

                  else
                     vb = bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			              bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						  bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						  bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						  bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			              bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						  bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						  bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 endif
						   
		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
				 btot    = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
				 
				 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

                 wbper2   = a_cross_b(vb,vb1)
                 wbper1   = a_cross_b(wbper2,vb)
                 wbper1   = wbper1/max(sqrt(wbper1(1)**2+wbper1(2)**2+wbper1(3)**2),1.e-8_p2) 
                 wbper2   = wbper2/max(sqrt(wbper2(1)**2+wbper2(2)**2+wbper2(3)**2),1.e-8_p2) 

! using v ------------------------------------------------------------
				 uu       = teste(mn)%block(m)%qv(4:6,L)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma
				  
				 vtot    = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
				 vpar    = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)

				 vper    = vv - vpar*vb
				 angle   = acos(vpar/max(vtot,1.e-6_p2))

                 Epar_   = vb(1)*ve(1) + vb(2)*vE(2) + vb(3)*vE(3)
				 Eper_   = vE - Epar_*vb
				 vdotEpar = Epar_*vpar
				 vdotEper = Eper_(1)*vper(1) + Eper_(2)*vper(2) + Eper_(3)*vper(3) 

				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)


				 ew(1)    = x2b(1,1)*ve(1) + x2b(1,2)*ve(2)+x2b(1,3)*ve(3)
				 ew(2)    = x2b(2,1)*ve(1) + x2b(2,2)*ve(2)+x2b(2,3)*ve(3)
				 ew(3)    = x2b(3,1)*ve(1) + x2b(3,2)*ve(2)+x2b(3,3)*ve(3)


                 csi_    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 csi_    = mod(csi_ + pi2*2  + pi ,pi2)
                 vper_   = vper_ * sqrt(bb0/max(bb,1.e-6))
				 vpar_   = sqrt(abs(vtot**2-vper_**2))
				 csi_    = vpar/max(vtot,1.e-6_p2)

				 vdotB   = vw(1)*bw(1) + vw(2)*bw(2) + vw(3)*bw(3) 
				 vdotE   = vw(1)*ew(1) + vw(2)*ew(2) + vw(3)*ew(3) 

                 wper1_   = wbper1(1)*vv(1) + wbper1(2)*vv(2) + wbper1(3)*vv(3)  
                 wper2_   = wbper2(1)*vv(1) + wbper2(2)*vv(2) + wbper2(3)*vv(3)
				   
                 theta_w  = atan2theta(wper2_,wper1_) 
                 theta_v  = atan2theta(vw(3),vw(2)) 



			     send(1,num)   = mark
			     send(2,num)   = xyz(1)
			     send(3,num)   = xyz(2)
			     send(4,num)   = xyz(3)
			     send(5,num)   = vtot
			     send(6,num)   = angle *180./pi
			     send(7,num)   = Btot
			     send(8,num)   = bw(1)      
			     send(9,num)   = sqrt(bw(2)**2+bw(3)**2)      
			     send(10,num)  = ew(1)      
			     send(11,num)  = sqrt(ew(2)**2+ew(3)**2)      
			     send(12,num)  = vdotE
			     send(13,num)  = vdotB

                 send(14,num)  = csi_ 
                 send(15,num)  = theta_v
                 send(16,num)  = theta_w

		       endif
	        enddo
         enddo
  enddo

     CALL MPI_ALLREDUCE(num, num0, 1, MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
     allocate(nsum(numberpe),nsum0(numberpe))
	 allocate(save(ind,num0))


     nsum           = 0
     nsum(mype + 1) = num
     CALL MPI_ALLREDUCE(nsum, nsum0,size(nsum),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
!    call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)

     if(mype/=0 .and. num >0) then
        call mpi_send(send(:,1:num), ind*num, mpi_real, &
  	                        0, mype+100, MPI_COMM_WORLD,IERROR)
     endif

     if(mype==0) then
        save(:,1:num) = send(:,1:num)
        do m=1,numberpe-1
	       if(nsum0(m+1) > 0) then
              num   = sum(nsum0(1:m))
		      allocate(recv(ind,nsum0(m+1)))
		      call mpi_recv(recv,size(recv),mpi_real, &
			                       m, m+100, MPI_COMM_WORLD,status,IERROR)
              save(:,num+1:num+nsum0(m+1))  = recv
			  deallocate(recv)
           endif
        enddo
        
		allocate(plot(ind, num0))
	    do m = 1,num0
           mark = int(save(1,m)/10)
		   plot(:,mark) = save(:,m) ! two informations in save(1,L) index*10 + mn
	    enddo

!  ----  shell apply ---------------------------
!        allocate (shell(num0),indi(num0),TMsav(ind,num0))
!		TMsav  = plot
!		shell  = plot(1,:)
!		call shellIs(num0,shell,1,indi)
!		do m=1,num0
!		   plot(:,m) = TMsav(:,indi(m))
!        enddo
!		deallocate(shell,indi,TMsav)
! -------------------------------------------
        do iplot = 1, 1
           call check( nf90_put_var(ncid, varid(7+iplot), plot, start = istart, &
                              count = icount) )
        enddo

        call check( nf90_close(ncid) )
        deallocate(plot)
	 endif 

  deallocate(send,save)
  deallocate(nsum,nsum0)


  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  call allocate_bvector(ev1_,-3)

  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_Orbit_testparticle1
! -------------------------------------------------------------------------------------------
! fv(vpar,vper),fv(vper1,vper2),fv(dvper1,dvper2),,fv(vpar,csi,vper),den(csi,q),  je,jb(q)
subroutine output_fvcsi_testparticles1(rstart)  
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use vector_functions
  use ionRingbeam_Case
  use trace_particles
  use test_particles

  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer    :: bv0_,bfld,bv1_,ev1_
  real,pointer,dimension(:,:,:,:) :: fve,fva,fvv,fvet,fvat,fvvt
  real,pointer,dimension(:,:,:)  :: plot,dene,dena,denv
  real,  pointer,dimension(:) :: ech,ach,vpar,Qarray

  real(p2) :: vminp,vmaxp,p,q,w,x,y,z,vmaxper,rstart
  real(p2) :: dcsi,dvpar,dvper,dvper1,scale
  real(8)  :: xyz(3),xyz0(3),xyz1(3),pqw(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),&
              vb(3),vb1(3),ve1(3),phase,det,vb_eq(3),B_eq,B_
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k,s,nregion,nre,mn,nalongQ1
  real(p2) :: v1,v2,v3,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,ww0,ww1
  real(p2) :: bw(3),vw(3),ew(3),vpar_,vper_ ,vper1_,vper2_,weight,gamma,dweight, &
              dmc1,dmc2,dmc3,dmc4,dmc5,dmc6,dmc7,ee,factor,dfactor,amass,vpermax,&
			  abw,aew,charge,bb,bb0,btot,pll,vb0(3),ee0,vtot,angle,angle0, &
			  Epar_,Eper_(3),vdotEpar,vdotEper,vdotB,csi_,vper(3),flux(13),dQ

  real(p2) :: vde_,vdb_, wbper1(3),wbper2(3),wx2b(3,3),wb2x(3,3),wper1_,wper2_,qflux(5), &
              upar_,uper1_,uper2_,uper_,upar_eq,uper_eq
  integer  :: ic,m,ierror,L,n,nbuff,ic0,iplot,iflag,kinde,ic0_d,num,numt
  integer  :: ind(20),index,kind,mark,ncsi,nsum,nvper,imc1,imc2,imc3,imc4,imc5,imc6,imc7

  character*5 fldname(33)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(60),varid(60), &
             dimids1(4),icount1(4),istart1(4),icount0(1),istart0(1), &
             dimids2(4),icount2(4),istart2(4),&
             dimids3(4),icount3(4),istart3(4)


  real  :: param(80)

  logical file_exist,flagi,flage

  save ic,ic0,ic0_d
  save dimid,varid,ncid,dimids1,icount1,istart1,icount0,istart0,&
       dimids2,icount2,istart2,dimids3,icount3,istart3

  data ic/0/,ic0/0/,ic0_d/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
! ---------------------------------------------------
  nalongQ1    = int(iftestParticle(195))
  if(nalongQ1 ==0) return

  call allocate_bvector(bfld,3)
  call allocate_bvector(bv0_,3)
  call allocate_bvector(bv1_,3)
  call allocate_bvector(ev1_,3)
  do m= 1, mblocks
     bfld(m)%vector     = bv0(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv0_)
  do m= 1, mblocks
     bfld(m)%vector     = bv1(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv1_)

  do m= 1, mblocks
     bfld(m)%vector     = ev1(m)%vector !+ bv1(m)%vector
	 if(nx<=2 .and. nz<=2) bfld(m)%vector(2,:,:,:) = 0.
  enddo

  call bcovcarE(bfld, ev1_)
  call allocate_bvector(bfld,-3)
! --------------------------------------------------------------
  allocate(ech(ne_tps),ach(na_tps))
  do i = 1,ne_tps
     ech(i) = Emin_tps + (i-1.)*de
  enddo
  do i = 1,na_tps
     ach(i) = amin_tps + (i-1.)*dalpha
  enddo

  nv         =  101
  vminp      = -2.5
  vmaxp      =  2.5
  dvpar      = (vmaxp-vminp)/real(nv-1)
  dQ         = ny/(nalongQ1-1.)

  allocate(vpar(nv),Qarray(nalongQ1))
  do i = 1, nv
     vpar(i) = vminp + (i-1)*dvpar
  enddo
  do i = 1, nalongQ1
     Qarray(i) = 0 + (i-1)*dQ 
	 imc1      = min(int(0 + (i-1.)*dQ) +1,ny)
	 dmc1      = imc1 - (0 + (i-1.)*dQ)
	 if(coordinate == -3 ) Qarray(i) = lenth_inq(imc1)*dmc1 + lenth_inq(imc1+1)*(1.-dmc1) &
	                                 - lenth_inq(ny/2+1) 
  enddo

  param      = 0.
  param(1)   = ne_tps
  param(2)   = na_tps
  param(3)   = emin_tps
  param(4)   = emax_tps
  param(5)   = amin_tps
  param(6)   = amax_tps
  param(7)   = vvmin(1)
  param(8)   = vvmax(1)
  param(9)   = dE
  param(10)  = dalpha
  param(11)  = nr_tps
  param(21)  = utime
  param(22)  = uspeed
  param(23)  = ut0
  param(24)  = cspeed
  param(25)  = ulength



  fldname = (/"DEe  ","DAe  ","DLe  ","DDEe ","DDAe ","DDLe ", &
              "VdE1e","VdE2e","VdBe ","Csie ","Dene ", &
			  "DEa  ","DAa  ","DLa  ","DDEa ","DDAa ","DDLa ", &
              "VdE1a","VdE2a","VdBa ","Csia ","Dena ", & 
			  "DEv  ","DAv  ","DLv  ","DDEv ","DDAv ","DDLv ", &
              "VdE1v","VdE2v","VdBv ","Csiv ","Denv "/)


  write(ct,'(I6.6)')(int((istep-1)*dt)/ndiagf_tps) * ndiagf_tps

  filename = 'TracingFtestparticles1_'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if((.not. file_exist) .or. abs(rstart)< 5* dt ) then
      ic0_save(42) = 0 
  endif

  ic0_save(42)= ic0_save(42) +1
  ic          = ic0_save(42)

! ------------------------------------------------------------------------------

  icount0 = (/1/)
  istart0 = (/ic/)
  istart1 = (/1, 1, 1,  ic   /)   !(nr,ne(a),nalongq,t)
  istart2 = (/1, 1, 1 , ic   /)   !(nr,ne(a),nalongq,t)
  istart3 = (/1, 1, 1 , ic   /)   !(nr,nv,nalongq,t)

  if(mype==0 ) then
        icount1 = (/nr_tps, ne_tps,   nalongQ1,   1 /)
        icount2 = (/nr_tps, na_tps,   nalongQ1,   1 /)
        icount3 = (/nr_tps, nv,       nalongQ1,   1 /)

!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'nr', nr_tps, dimid(1)) )
          call check( nf90_def_var(ncid, 'nr', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nalongQ1', nalongQ1, dimid(2)) )
          call check( nf90_def_var(ncid, 'nalongQ1', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'q') )


          call check( nf90_def_dim(ncid, 'nE', ne_tps, dimid(3)) )
          call check( nf90_def_var(ncid, 'nE', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'eV') )


          call check( nf90_def_dim(ncid, 'nvpar', nv, dimid(7)) )
          call check( nf90_def_var(ncid, 'nvpar', NF90_REAL, dimid(7), varid(7)) )
          call check( nf90_put_att(ncid,  varid(7), UNITS, 'm/s') )

          call check( nf90_def_dim(ncid, 'nalpha', na_tps, dimid(4)) )
          call check( nf90_def_var(ncid, 'nalpha', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'deg') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )

! -----------------------------------------------

          dimids1  = (/dimid(1),dimid(3),dimid(2),dimid(5)/)
          dimids2  = (/dimid(1),dimid(4),dimid(2),dimid(5)/)
          dimids3  = (/dimid(1),dimid(7),dimid(2),dimid(5)/)
          do iplot = 1, 11
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids1, varid(10+iplot)) )
           call check( nf90_put_att(ncid, varid(10+iplot), UNITS, 'flux') )
          enddo

          do iplot = 1, 11
           call check( nf90_def_var(ncid, fldname(iplot+11), NF90_REAL, dimids2, varid(21+iplot)) )
           call check( nf90_put_att(ncid, varid(21+iplot), UNITS, 'flux') )
          enddo

          do iplot = 1, 11
           call check( nf90_def_var(ncid, fldname(iplot+22), NF90_REAL, dimids3, varid(32+iplot)) )
           call check( nf90_put_att(ncid, varid(32+iplot), UNITS, 'flux') )
          enddo


          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(3), ech) )
          call check( nf90_put_var(ncid, varid(4), ach) )
          call check( nf90_put_var(ncid, varid(2), Qarray) )
          call check( nf90_put_var(ncid, varid(7), vpar) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart0) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart0) )
          
		  do iplot = 1, 11
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(10+iplot)) )
          enddo
          do iplot = 1, 11
             call check( nf90_inq_varid(ncid, fldname(iplot+11), varid(21+iplot)) )
          enddo
          do iplot = 1, 11
             call check( nf90_inq_varid(ncid, fldname(iplot+22), varid(32+iplot)) )
          enddo

	    endif
  endif
! ---------------------------------------------------------------------------

  deallocate(ech,ach,Qarray,vpar)

  allocate(fve(13,nr_tps,ne_tps,nalongQ1) , fvet(13,nr_tps,ne_tps,nalongQ1))
  allocate(fva(13,nr_tps,na_tps,nalongQ1),  fvat(13,nr_tps,na_tps,nalongQ1))
  allocate(fvv(13,nr_tps,nv,nalongQ1),      fvvt(13,nr_tps,nv,nalongQ1))

  fve   = 0.
  fva   = 0.
  fvv   = 0.
  num   = 0

  do mn = 1,nr_tps
     do m=1,mblocks
           do L = 1,teste(mn)%block(m)%mi
			  amass = teste(mn)%amass
	          charge   = teste(mn)%charge 
		      nre    = mn
			  p      = teste(mn)%block(m)%qv(1,L)
			  q      = teste(mn)%block(m)%qv(2,L)
			  w      = teste(mn)%block(m)%qv(3,L)
			  weight = teste(mn)%block(m)%qv(7,L)

			  wz1      = w +1.0
              kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	          kp       = min(block(m)%nz1,kk+1)
			  wz1      = wz1 - kk
	          wz0      = 1.0 - wz1

			  wx1      = p +1.0
              ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	          ip       = min(block(m)%nx1,ii+1)
			  wx1      = wx1 - ii
	          wx0      = 1.0 - wx1

			  wy1      = q +1.0
              jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	          jp       = min(block(m)%ny1,jj+1)
			  wy1      = wy1 - jj
	          wy0      = 1.0 - wy1

              vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
              vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
              ve1   = ev1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             ev1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 ev1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 ev1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 ev1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             ev1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 ev1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 ev1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

             if(dipole) then
		        xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !
                x    = xyz(1)
                y    = xyz(2)
                z    = xyz(3)
			    vb   = dipole3DB(x,y,z)
				bb   = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
		        xyz0 = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !
                x    = xyz0(1)
                y    = xyz0(2)
                z    = xyz0(3)
			    vb0  = dipole3DB(x,y,z)
				bb0  = sqrt(vb0(1)**2+vb0(2)**2+vb0(3)**2)
			 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
				 btot    = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
				 uu      = teste(mn)%block(m)%qv(4:6,L)/amass
! --------- check D_LL ----------------				 
	             xyz1    = a_cross_b(uu,vb)/bb/charge*amass
	             x       = xyz(1)+xyz1(1)
	             y       = xyz(2)+xyz1(2)
	             z       = xyz(3)+xyz1(3)
	             pqw     = xyz2pqw(0,x,y,z)

	             pLL     = (pqw(1)-pmin)/deltax

! ---------------------------------
            if(dipole) then
		        xyz  = pqw2xyz(0,pmin+pLL*deltax,qmin+q*deltay,wmin+w*deltaz) !
                x    = xyz(1)
                y    = xyz(2)
                z    = xyz(3)
			    vb   = dipole3DB(x,y,z)
				bb   = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
		        xyz0 = pqw2xyz(0,pmin+pLL*deltax,qequator,wmin+w*deltaz) !
                x    = xyz0(1)
                y    = xyz0(2)
                z    = xyz0(3)
			    vb0  = dipole3DB(x,y,z)
				bb0  = sqrt(vb0(1)**2+vb0(2)**2+vb0(3)**2)
		    endif
		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
				 btot    = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)

! ---------------------------------



				 vbper1  = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

                 wbper2   = a_cross_b(vb,vb1)
                 wbper1   = a_cross_b(wbper2,vb)
                 wbper1   = wbper1/max(sqrt(wbper1(1)**2+wbper1(2)**2+wbper1(3)**2),1.e-8_p2) 
                 wbper2   = wbper2/max(sqrt(wbper2(1)**2+wbper2(2)**2+wbper2(3)**2),1.e-8_p2) 
! --------------------------------------------
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
		         ee       = (gamma-1.)*cspeed**2*amass * uT0  
		         ee       = log10(ee)
		         ee0      = teste(mn)%block(m)%qv(8,L)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma
				 vtot     = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
				 vpar_    = vv(1)*vb(1) + vv(2)*vb(2) + vv(3)*vb(3) 
				 vper     = vv - vpar_*vb
				 vper_    = sqrt(max(vtot**2-vpar_**2,0._p2)) * sqrt(BB0/max(BB,1.e-10))
			     vpar_    = vpar_/max(abs(vpar_),1.e-10_p2)*sqrt(abs(vtot**2-vper_**2))
				 angle    = vpar_/max(vtot,1.e-12)
				 if(ifrad) angle = acos(angle)*180./pi

!				 angle    = acos( max(min((vpar)/max(vtot,1.e-6),1.),-1.)) *180./pi   ![0., 90]
! check equatorial pitch angle ------------
                 

                 angle0   = teste(mn)%block(m)%qv(9,L)   ![-1,1]

                 Epar_    = vb(1)*ve1(1) + vb(2)*vE1(2) + vb(3)*vE1(3)
				 Eper_    = vE1 - Epar_*vb

				 vdotEpar = Epar_*vpar_
				 vdotEper = Eper_(1)*vper(1) + Eper_(2)*vper(2) + Eper_(3)*vper(3) 
				 vdotB    = vb1(1)*vv(1) + vb1(2)*vv(2) + vb1(3)*vv(3) 

				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)


				 ew(1)    = x2b(1,1)*ve1(1) + x2b(1,2)*ve1(2)+x2b(1,3)*ve1(3)
				 ew(2)    = x2b(2,1)*ve1(1) + x2b(2,2)*ve1(2)+x2b(2,3)*ve1(3)
				 ew(3)    = x2b(3,1)*ve1(1) + x2b(3,2)*ve1(2)+x2b(3,3)*ve1(3)


                 csi_    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 csi_    = mod(csi_ + pi2*2,  pi2)
!   ----------------------------------------------------------------------------------------------------------
                 imc1     = min(max(int(max(ee0-emin_tps,0.)/de +1.),1),ne_tps-1)
		         dmc1     = imc1-(ee0 - emin_tps)/de
                 imc2     = min(max(int(max(angle0-amin_tps,0.)/dalpha +1.),1),na_tps-1)
		         dmc2     = imc2-(angle0 - amin_tps)/dalpha
                 imc4     = min(int((q-qmin)/dQ + 1.),nalongQ1-1)
                 dmc4     = imc4 - (q-qmin)/dQ

				 flux(1)  = 10.**ee - 10.**ee0   !4 DDe
				 flux(2)  = abs(angle) - abs(angle0)       !4 DDa
                 if(ifrad) flux(2) = acos(abs(cos(angle*pi/180.)))*180./pi - &
				                     acos(abs(cos(angle0*pi/180.)))*180./pi
				 flux(3)  = (pLL - teste(mn)%block(m)%qv(10,L))*deltax/ulength
				 flux(4)  = flux(1)**2
				 flux(5)  = flux(2)**2
				 flux(6)  = flux(3)**2
				 flux(7)  = vdotEpar
				 flux(8)  = vdotEper
				 flux(9)  = vdotB
				 flux(10) = csi_
				 flux(11) = 1.


                 fve(1:11,mn,imc1,imc4)     = fve(1:11,mn,imc1,imc4)    + dmc1*dmc4*flux(1:11)	          
                 fve(1:11,mn,imc1,imc4+1)   = fve(1:11,mn,imc1,imc4+1)  + dmc1*(1.-dmc4)*flux(1:11)		          
                 fve(1:11,mn,imc1+1,imc4)   = fve(1:11,mn,imc1+1,imc4)  + (1.-dmc1)*dmc4*flux(1:11)		          
                 fve(1:11,mn,imc1+1,imc4+1) = fve(1:11,mn,imc1+1,imc4+1)+ (1.-dmc1)*(1.-dmc4)*flux(1:11)	          

                 fva(1:11,mn,imc2,imc4)     = fva(1:11,mn,imc2,imc4)    + dmc2*dmc4*flux(1:11)	          
                 fva(1:11,mn,imc2,imc4+1)   = fva(1:11,mn,imc2,imc4+1)  + dmc2*(1.-dmc4)*flux(1:11)		          
                 fva(1:11,mn,imc2+1,imc4)   = fva(1:11,mn,imc2+1,imc4)  + (1.-dmc2)*dmc4*flux(1:11)		          
                 fva(1:11,mn,imc2+1,imc4+1) = fva(1:11,mn,imc2+1,imc4+1)+ (1.-dmc2)*(1.-dmc4)*flux(1:11)	          

                 if(vpar_ >= vminp .and. vpar_< vmaxp) then
                    imc3     = min(int((vpar_-vminp)/dvpar + 1.),nv-1)
                    dmc3     = imc3 - (vpar_-vminp)/dvpar

                    fvv(1:11,mn,imc3,imc4)     = fvv(1:11,mn,imc3,imc4)    + dmc3*dmc4*flux(1:11)	          
                    fvv(1:11,mn,imc3,imc4+1)   = fvv(1:11,mn,imc3,imc4+1)  + dmc3*(1.-dmc4)*flux(1:11)		          
                    fvv(1:11,mn,imc3+1,imc4)   = fvv(1:11,mn,imc3+1,imc4)  + (1.-dmc3)*dmc4*flux(1:11)		          
                    fvv(1:11,mn,imc3+1,imc4+1) = fvv(1:11,mn,imc3+1,imc4+1)+ (1.-dmc3)*(1.-dmc4)*flux(1:11)
				 endif	          
                 num   = num + 1


221           continue
		  enddo
     enddo
  enddo

  CALL MPI_ALLREDUCE(fve,fvet,size(fve),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fva,fvat,size(fva),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  CALL MPI_ALLREDUCE(fvv,fvvt,size(fvv),MPI_real,MPI_SUM,mpi_comm_world,IERROR)
  fvat(:,:,na_tps-1,:) = (fvat(:,:,1,:) + fvat(:,:,na_tps-1,:))/2. 
  fvat(:,:,1,:)        =  fvat(:,:,na_tps-1,:) 
  fvat(:,:,na_tps,:)   =  fvat(:,:,2,:)

  allocate(dene(nr_tps,ne_tps,nalongQ1))
  allocate(dena(nr_tps,na_tps,nalongQ1))
  allocate(denv(nr_tps,nv,nalongQ1))
  dene    = max(fvet(11,:,:,:),1.e-12)   
  dena    = max(fvat(11,:,:,:),1.e-12) 
  denv    = max(fvvt(11,:,:,:),1.e-12)


  if(mype==0) then
     allocate(plot(nr_tps,ne_tps,nalongQ1))
     do iplot = 1, 11
        if(iplot == 1)  plot = fvet(1,:,:,:)/dene
        if(iplot == 2)  plot = fvet(2,:,:,:)/dene
        if(iplot == 3)  plot = fvet(3,:,:,:)/dene
        if(iplot == 4)  plot = fvet(4,:,:,:)/dene - (fvet(1,:,:,:)/dene)**2  
        if(iplot == 5)  plot = fvet(5,:,:,:)/dene - (fvet(2,:,:,:)/dene)**2
        if(iplot == 6)  plot = fvet(6,:,:,:)/dene - (fvet(3,:,:,:)/dene)**2
        if(iplot == 7)  plot = fvet(7,:,:,:)/dene
        if(iplot == 8)  plot = fvet(8,:,:,:)/dene
        if(iplot == 9)  plot = fvet(9,:,:,:)/dene
        if(iplot == 10) plot = fvet(10,:,:,:)/dene
        if(iplot == 11) plot = fvet(11,:,:,:)/iftestParticle(6)

        call check( nf90_put_var(ncid, varid(10+iplot), plot, start = istart1, &
                              count = icount1) )
     enddo
     deallocate(plot)
     allocate(plot(nr_tps,na_tps,nalongQ1))
     do iplot = 1, 11
        if(iplot == 1)  plot = fvat(1,:,:,:)/dena
        if(iplot == 2)  plot = fvat(2,:,:,:)/dena
        if(iplot == 3)  plot = fvat(3,:,:,:)/dena
        if(iplot == 4)  plot = fvat(4,:,:,:)/dena - (fvat(1,:,:,:)/dena)**2  
        if(iplot == 5)  plot = fvat(5,:,:,:)/dena - (fvat(2,:,:,:)/dena)**2
        if(iplot == 6)  plot = fvat(6,:,:,:)/dena - (fvat(3,:,:,:)/dena)**2
        if(iplot == 7)  plot = fvat(7,:,:,:)/dena
        if(iplot == 8)  plot = fvat(8,:,:,:)/dena
        if(iplot == 9)  plot = fvat(9,:,:,:)/dena
        if(iplot == 10) plot = fvat(10,:,:,:)/dena
        if(iplot == 11) plot = fvat(11,:,:,:)/iftestParticle(6)

        call check( nf90_put_var(ncid, varid(21+iplot), plot, start = istart2, &
                              count = icount2) )
     enddo
     deallocate(plot)
     allocate(plot(nr_tps,nv,nalongQ1))
     do iplot = 1, 11
        if(iplot == 1)  plot = fvvt(1,:,:,:)/denv
        if(iplot == 2)  plot = fvvt(2,:,:,:)/denv
        if(iplot == 3)  plot = fvvt(3,:,:,:)/denv
        if(iplot == 4)  plot = fvvt(4,:,:,:)/denv - (fvvt(1,:,:,:)/denv)**2  
        if(iplot == 5)  plot = fvvt(5,:,:,:)/denv - (fvvt(2,:,:,:)/denv)**2
        if(iplot == 6)  plot = fvvt(6,:,:,:)/denv - (fvvt(3,:,:,:)/denv)**2
        if(iplot == 7)  plot = fvvt(7,:,:,:)/denv
        if(iplot == 8)  plot = fvvt(8,:,:,:)/denv
        if(iplot == 9)  plot = fvvt(9,:,:,:)/denv
        if(iplot == 10) plot = fvvt(10,:,:,:)/denv
        if(iplot == 11) plot = fvvt(11,:,:,:)/iftestParticle(6)

        call check( nf90_put_var(ncid, varid(32+iplot), plot, start = istart3, &
                              count = icount3) )
     enddo
     deallocate(plot)
  endif



  if(mype==0) call check( nf90_close(ncid) )


  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  call allocate_bvector(ev1_,-3)
  deallocate(dene,dena,denv,fve,fvet,fva,fvat,fvv,fvvt)
  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_fvcsi_testparticles1

! -----------------------------------------------------------------------------------
! fv(nx,vpar, nalong q) of E*V, B*V, E 
subroutine output_4fvparL_ncfile(rstart)  
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use vector_functions
  use ionRingbeam_Case
  use trace_particles
  use test_particles

  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer    :: bv0_,bfld,bv1_,ev1_
  type(particle_type)   :: particle
  real,pointer,dimension(:,:,:,:,:) :: fv1,fv1t
  real,pointer,dimension(:,:,:,:)  :: plot,den
  real,  pointer,dimension(:) :: vpar,Qarray

  real(p2) :: vminp,vmaxp,p,q,w,x,y,z,vmaxper,rstart
  real(p2) :: dcsi,dvpar,dvper,dvper1,scale
  real(8)  :: xyz(3),xyz0(3),xyz1(3),pqw(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),&
              vb(3),vb1(3),ve1(3),phase,det,vb_eq(3),B_eq,B_,mu
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k,s,nregion,nre,mn,kinde,nalongQ1
  real(p2) :: v1,v2,v3,dQ,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,ww0,ww1
  real(p2) :: bw(3),vw(3),ew(3),vpar_,vper_ ,vper1_,vper2_,weight,gamma,dweight, &
              dmc1,dmc2,dmc3,dmc4,dmc5,dmc6,dmc7,ee,factor,dfactor,amass,vpermax,&
			  abw,aew,charge,bb,bb0,btot,pll,vb0(3),ee0,vtot,angle,angle0, &
			  Epar_,Eper_(3),vdotEpar,vdotEper,vdotB,csi_,qflux(13),pval(10),dpval(10)

  real(p2) :: vde1_,vde2_,vdb_, wbper1(3),wbper2(3),wx2b(3,3),wb2x(3,3),wper1_,wper2_, &
              upar_,uper1_,uper2_,uper_,upar_eq,uper_eq,frac
  integer  :: ic,m,ierror,L,n,nbuff,ic0,iplot,iflag,ic0_d,num,numt
  integer  :: ind(20),index,kind,mark,ncsi,nsum,nvper,imc1,imc2,imc3,imc4,imc5,imc6,imc7

  character*5 fldname(6)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(40),varid(40), &
             dimids1(5),icount1(5),istart1(5),icount0(1),istart0(1)
  real  :: param(80)
  logical file_exist,flagi,flage
  save ic,ic0,ic0_d
  save dimid,varid,ncid,dimids1,icount1,istart1,icount0,istart0

  data ic/0/,ic0/0/,ic0_d/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
! ---------------------------------------------------
  call allocate_bvector(bfld,3)
  call allocate_bvector(bv0_,3)
  call allocate_bvector(bv1_,3)
  call allocate_bvector(ev1_,3)
  do m= 1, mblocks
     bfld(m)%vector     = bv0(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv0_)
  do m= 1, mblocks
     bfld(m)%vector     = bv1(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv1_)
  call bcovcarE(ev1, ev1_)
  call allocate_bvector(bfld,-3)
! --------------------------------------------------------------
  kinde       = int(pvdomain1(3))
  if(kinde == 0) return

  nre         = int(pvdomain1(91))
  if(nre  == 0) return

  nalongQ1    = int(pvdomain1(92))
  nv          = 101
  vmaxp       =  3.5
  vminp       = -3.5
  dvpar       = (vmaxp-vminp)/real(nv-1)
  dQ          = ny/(nalongQ1-1.)
  L           = 93
  do n=1,nre
     pval(n)  = pvdomain1(L)*nx
     dpval(n) = pvdomain1(L+1)
	 L        = L + 2
  enddo

  do m=1,kinde
     param(10+m)           =  int(pvdomain(3+m))	
     do n = 1, kinds
	    if(	int(param(10+m)) == n   ) ind(n)    = m  ! ions
	    if(	int(param(10+m)) == n+10) ind(n+10) = m  ! eles
	 enddo	 
  enddo

  
  allocate(vpar(nv),Qarray(nalongQ1))
  do i = 1, nv
     vpar(i) = vminp + (i-1)*dvpar
  enddo
  do i = 1, nalongQ1
     Qarray(i) = 0 + (i-1)*dQ 
	 imc1      = min(int(0 + (i-1.)*dQ) +1,ny)
	 dmc1      = imc1 - (0 + (i-1.)*dQ)
	 if(coordinate == -3 ) Qarray(i) = lenth_inq(imc1)*dmc1 + lenth_inq(imc1+1)*(1.-dmc1) &
	                                 - lenth_inq(ny/2+1) 
  enddo



  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed
  param(8)   = va
  param(9)   = alpha
  param(10)  = kinds
  param(11:kinds+10)      = mions(1:kinds)

  fldname = (/"vde1 ","vde2 ","vdB  ","Csi  ","ee   ","mu   " /)


  write(ct,'(I6.6)')(int((istep-1)*dt)/max(int(pvdomain1(99)),1) ) * max(int(pvdomain1(99)),1)

  filename = 'FvparQ'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if((.not. file_exist) .or. abs(rstart)< 5* dt ) then
      ic0_save(51) = 0 
  endif

  ic0_save(51)= ic0_save(51) +1
  ic          = ic0_save(51)

! ------------------------------------------------------------------------------

  icount0 = (/1/)
  istart0 = (/ic/)
  istart1 = (/1, 1, 1, 1, ic   /)   !(nkind,nre,vpar,nalongq,t)

  if(mype==0 ) then
        icount1 = (/kinde,nre, nv, nalongQ1, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'kind', kinde, dimid(1)) )
          call check( nf90_def_var(ncid, 'kind', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nregion', nre, dimid(2)) )
          call check( nf90_def_var(ncid, 'nregion', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'q') )


          call check( nf90_def_dim(ncid, 'nvpar', nv, dimid(3)) )
          call check( nf90_def_var(ncid, 'nvpar', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'eV') )

          call check( nf90_def_dim(ncid, 'nalongQ1', nalongQ1, dimid(4)) )
          call check( nf90_def_var(ncid, 'nalongQ1', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'deg') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )

! -----------------------------------------------

          dimids1  = (/dimid(1),dimid(2),dimid(3),dimid(4),dimid(5)/)
          do iplot = 1, 6
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids1, varid(10+iplot)) )
           call check( nf90_put_att(ncid, varid(10+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(3), vpar) )
          call check( nf90_put_var(ncid, varid(4), Qarray) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart0) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart0) )
          
		  do iplot = 1, 6
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(10+iplot)) )
          enddo
	    endif
  endif
! ---------------------------------------------------------------------------

  allocate(fv1(7,kinde,nre,nv,nalongQ1) , fv1t(7,kinde,nre,nv,nalongQ1))
  fv1   = 0.

  do mn=1,kinds
     particle = eles(mn)
	 amass    = particle%amass
	 frac     = particle%frac
	 scale    = particle%vth
	 if(ind(mn+10) ==0) goto 220
     do m=1,mblocks
           do L = 1,particle%block(m)%mi

                 p      = particle%block(m)%qv(1,L)
                 q      = particle%block(m)%qv(2,L)
                 w      = particle%block(m)%qv(3,L)
			     weight   = particle%block(m)%qv(13,L) * particle%block(m)%qv(14,L)  * frac

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1

                 vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 ve1   = ev1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             ev1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 ev1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 ev1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 ev1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             ev1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 ev1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 ev1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
					 B_   = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
!		             xyz  = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !get magnetic field at equator
!                     x    = xyz(1)
!                     y    = xyz(2)
!                     z    = xyz(3)
!			         vb_eq= dipole3DB(x,y,z)
!					 B_eq = sqrt(vb_eq(1)**2+vb_eq(2)**2+vb_eq(3)**2)
				 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb       = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2) 
                 vbper1   = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

				 uu       = particle%block(m)%qv(4:6,L)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 ee       = amass*cspeed**2*(gamma -1.)*ut0

				 vv       = uu /gamma     ! scale factor here
				 vv       = vv /scale 

				 vw(1)    = x2b(1,1)*uu(1) + x2b(1,2)*uu(2)+x2b(1,3)*uu(3)
				 vw(2)    = x2b(2,1)*uu(1) + x2b(2,2)*uu(2)+x2b(2,3)*uu(3)
				 vw(3)    = x2b(3,1)*uu(1) + x2b(3,2)*uu(2)+x2b(3,3)*uu(3)
				 mu       = (vw(2)**2+vw(3)**2)/2./B_

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)
                 angle    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 angle    = mod(angle + pi2*2,  pi2)    ! southward ????  --- csi ----


				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 ew(1)    = x2b(1,1)*ve1(1) + x2b(1,2)*ve1(2)+x2b(1,3)*ve1(3)
				 ew(2)    = x2b(2,1)*ve1(1) + x2b(2,2)*ve1(2)+x2b(2,3)*ve1(3)
				 ew(3)    = x2b(3,1)*ve1(1) + x2b(3,2)*ve1(2)+x2b(3,3)*ve1(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

! ----------------------------------------------------------------------------------
                 vde1_    = vw(1) * ew(1)
                 vde2_    = vw(2) * ew(2) + vw(3) * ew(3)
                 vdb_     = a_dot_b(vv, vb1) /max(sqrt(vb1(1)**2+vb1(2)**2+vb1(3)**2),1.e-8_p2)

				 qflux(1:7)    = (/vde1_,vde2_,vdb_, angle, ee, mu, 1._p2/)

                 if(vpar_ < vminp .or. vpar_ >=vmaxp ) goto 221 
                 do n = 1, nre
				    if(abs(p-pval(n)) <= dpval(n)) then
                       imc1  = min(max(int(max(vpar_ - vminp,0.)/dvpar + 1.),1),nv-1)
		               dmc1  = imc1-(vpar_ - vminp)/dvpar
                       imc2  = min(int(q/dQ + 1.),nalongQ1-1)
		               dmc2  = imc2 - q/dQ


                       fv1(1:7,mn,n,imc1,imc2)     = fv1(1:7,mn,n,imc1,imc2)    + dmc1*dmc2*qflux(1:7)	          
                       fv1(1:7,mn,n,imc1,imc2+1)   = fv1(1:7,mn,n,imc1,imc2+1)  + dmc1*(1.-dmc2)*qflux(1:7)		          
                       fv1(1:7,mn,n,imc1+1,imc2)   = fv1(1:7,mn,n,imc1+1,imc2)  + (1.-dmc1)*dmc2*qflux(1:7)		          
                       fv1(1:7,mn,n,imc1+1,imc2+1) = fv1(1:7,mn,n,imc1+1,imc2+1)+ (1.-dmc1)*(1.-dmc2)*qflux(1:7)	          


					endif
                 enddo

221           continue
		  enddo
     enddo
220  continue
  enddo

  CALL MPI_ALLREDUCE(fv1,fv1t,size(fv1),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  allocate(den(kinde,nre,nv,nalongQ1))
  
  den   = max(fv1t(7,:,:,:,:),1.e-12)  !/iftestParticle(6)


  if(mype==0) then
     allocate(plot(kinde,nre,nv,nalongQ1))
     do iplot = 1, 6
        plot  = fv1t(iplot,:,:,:,:)/den
        call check( nf90_put_var(ncid, varid(10+iplot), plot, start = istart1, &
                              count = icount1) )
     enddo
     deallocate(plot)
  endif



  if(mype==0) call check( nf90_close(ncid) )


  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  call allocate_bvector(ev1_,-3)
  deallocate(den,fv1,fv1t)
  deallocate(vpar,Qarray)
  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_4fvparL_ncfile
! -------------------------------------------------------
! -----------------------------------------------------------------------------------
! fv(kinds,vpar, nalong q) of E*V, B*V, E  ====== for ions =====
! (/kinde,nre, nv, nalongQ1, 1 /)
subroutine output_4fvparL_ncfile_case58(rstart)  
  use netcdf
  use global_parameters
  use grid_data
  use diagnos
  use vector_functions
  use ionRingbeam_Case
  use trace_particles
  use test_particles

  implicit none
  include 'mpif.h'
  type(bvector_type), dimension(:), pointer    :: bv0_,bfld,bv1_,ev1_
  type(particle_type)   :: particle
  real,pointer,dimension(:,:,:,:,:) :: fv1,fv1t
  real,pointer,dimension(:,:,:,:)  :: plot,den
  real,  pointer,dimension(:) :: vpar,Qarray

  real(p2) :: vminp,vmaxp,p,q,w,x,y,z,vmaxper,rstart
  real(p2) :: dcsi,dvpar,dvper,dvper1,scale
  real(8)  :: xyz(3),xyz0(3),xyz1(3),pqw(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),&
              vb(3),vb1(3),ve1(3),phase,det,vb_eq(3),B_eq,B_,mu
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k,s,nregion,nre,mn,kinde,nalongQ1
  real(p2) :: v1,v2,v3,dQ,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,ww0,ww1
  real(p2) :: bw(3),vw(3),ew(3),vpar_,vper_ ,vper1_,vper2_,weight,gamma,dweight, &
              dmc1,dmc2,dmc3,dmc4,dmc5,dmc6,dmc7,ee,factor,dfactor,amass,vpermax,&
			  abw,aew,charge,bb,bb0,btot,pll,vb0(3),ee0,vtot,angle,angle0, &
			  Epar_,Eper_(3),vdotEpar,vdotEper,vdotB,csi_,qflux(13),pval(10),dpval(10)

  real(p2) :: vde1_,vde2_,vdb_, wbper1(3),wbper2(3),wx2b(3,3),wb2x(3,3),wper1_,wper2_, &
              upar_,uper1_,uper2_,uper_,upar_eq,uper_eq,frac
  integer  :: ic,m,ierror,L,n,nbuff,ic0,iplot,iflag,ic0_d,num,numt
  integer  :: ind(20),index,kind,mark,ncsi,nsum,nvper,imc1,imc2,imc3,imc4,imc5,imc6,imc7

  character*5 fldname(7)
  character*6 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(40),varid(40), &
             dimids1(5),icount1(5),istart1(5),icount0(1),istart0(1)
  real  :: param(80)
  logical file_exist,flagi,flage
  save ic,ic0,ic0_d
  save dimid,varid,ncid,dimids1,icount1,istart1,icount0,istart0

  data ic/0/,ic0/0/,ic0_d/0/

  interface 
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB
      subroutine bcovcarE(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bcovcarE

      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
  end interface
! output the ion flux distribution as function of flux(x,energy,alpha,time), where alpha is the pitch ange
! ---------------------------------------------------
  call allocate_bvector(bfld,3)
  call allocate_bvector(bv0_,3)
  call allocate_bvector(bv1_,3)
  call allocate_bvector(ev1_,3)
  do m= 1, mblocks
     bfld(m)%vector     = bv0(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv0_)
  do m= 1, mblocks
     bfld(m)%vector     = bv1(m)%vector !+ bv1(m)%vector
  enddo
  call bconcarB(bfld,bv1_)
  call bcovcarE(ev1, ev1_)
  call allocate_bvector(bfld,-3)
! --------------------------------------------------------------
  kinde       = int(pvdomain1(3))
  if(kinde == 0) return

  nre         = int(pvdomain1(91))
  if(nre  == 0) return

  nalongQ1    = int(pvdomain1(92))
  nv          = 101
  vmaxp       =  5
  vminp       = -5
  dvpar       = (vmaxp-vminp)/real(nv-1)
  dQ          = nx/(nalongQ1-1.)

  L           = 93
  do n=1,nre
     pval(n)  = pvdomain1(L)*nx
     dpval(n) = pvdomain1(L+1)
	 L        = L + 2
  enddo

  do m=1,kinde
     param(10+m)           =  int(pvdomain(3+m))	
     do n = 1, kinds
	    if(	int(param(10+m)) == n   ) then
		    ind(n)    = m  ! ions
			flagi      = .true.
			flage      = .false.
        endif
	    if(	int(param(10+m)) == n+10) then
		    ind(n+10) = m  ! eles
			flage      = .true.
			flagi      = .false.
        endif
	 enddo	 
  enddo

  
  allocate(vpar(nv),Qarray(nalongQ1))
  do i = 1, nv
     vpar(i) = vminp + (i-1)*dvpar
  enddo
  do i = 1, nalongQ1
     Qarray(i) = 0 + (i-1)*dQ 
	 imc1      = min(int(0 + (i-1.)*dQ) +1,nx)
	 dmc1      = imc1 - (0 + (i-1.)*dQ)
	 if(coordinate == -3 ) Qarray(i) = lenth_inq(imc1)*dmc1 + lenth_inq(imc1+1)*(1.-dmc1) &
	                                 - lenth_inq(ny/2+1) 
  enddo



  param      = 0.
  param(1)   = utime
  param(2)   = uspeed
  param(3)   = ut0
  param(4)   = va
  param(5)   = ulength
  param(6)   = uRe
  param(7)   = cspeed
  param(8)   = va
  param(9)   = alpha
  param(10)  = kinds
  param(11:kinds+10)      = mions(1:kinds)

  fldname = (/"vde1 ","vde2 ","vdB  ","Csi  ","ee   ","mu   ","fv   " /)


  write(ct,'(I6.6)')(int((istep-1)*dt)/max(int(pvdomain1(99)),1) ) * max(int(pvdomain1(99)),1)

  filename = 'FvparQ'//ct//'.nc'


  inquire(file=filename,exist=file_exist)
  if((.not. file_exist) .or. abs(rstart)< 5* dt ) then
      ic0_save(51) = 0 
  endif

  ic0_save(51)= ic0_save(51) +1
  ic          = ic0_save(51)

! ------------------------------------------------------------------------------

  icount0 = (/1/)
  istart0 = (/ic/)
  istart1 = (/1, 1, 1, 1, ic   /)   !(nkind,nre,vpar,nalongq,t)

  if(mype==0 ) then
        icount1 = (/kinde,nre, nv, nalongQ1, 1 /)
!       Create the file. 
        if(ic == 1 ) then
          call check( nf90_create(FILENAME, nf90_clobber, ncid) )
! define parameters
          call check( nf90_def_dim(ncid, 'kind', kinde, dimid(1)) )
          call check( nf90_def_var(ncid, 'kind', NF90_REAL, dimid(1), varid(1)) )
          call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

          call check( nf90_def_dim(ncid, 'nregion', nre, dimid(2)) )
          call check( nf90_def_var(ncid, 'nregion', NF90_REAL, dimid(2), varid(2)) )
          call check( nf90_put_att(ncid,  varid(2), UNITS, 'q') )


          call check( nf90_def_dim(ncid, 'nvpar', nv, dimid(3)) )
          call check( nf90_def_var(ncid, 'nvpar', NF90_REAL, dimid(3), varid(3)) )
          call check( nf90_put_att(ncid,  varid(3), UNITS, 'eV') )

          call check( nf90_def_dim(ncid, 'nalongQ1', nalongQ1, dimid(4)) )
          call check( nf90_def_var(ncid, 'nalongQ1', NF90_REAL, dimid(4), varid(4)) )
          call check( nf90_put_att(ncid,  varid(4), UNITS, 'deg') )

          call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(5)) )
          call check( nf90_def_var(ncid, 'time', NF90_REAL, dimid(5), varid(5)) )

          call check( nf90_def_dim(ncid, 'param', 80, dimid(6)) )
          call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(6), varid(6)) )
          call check( nf90_put_att(ncid,  varid(6), UNITS, 'non') )

! -----------------------------------------------

          dimids1  = (/dimid(1),dimid(2),dimid(3),dimid(4),dimid(5)/)
          do iplot = 1, 7
           call check( nf90_def_var(ncid, fldname(iplot), NF90_REAL, dimids1, varid(10+iplot)) )
           call check( nf90_put_att(ncid, varid(10+iplot), UNITS, 'flux') )
          enddo

          call check( nf90_enddef(ncid) )

          call check( nf90_put_var(ncid, varid(6), param) )
          call check( nf90_put_var(ncid, varid(3), vpar) )
          call check( nf90_put_var(ncid, varid(4), Qarray) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart0) )

        
		else
          call check( nf90_open(FILENAME, nf90_write, ncid) )
          call check( nf90_inq_varid(ncid, 'time', varid(5)) )
          call check( nf90_put_var(ncid, varid(5), stime*aion, start = istart0) )
          
		  do iplot = 1, 7
             call check( nf90_inq_varid(ncid, fldname(iplot), varid(10+iplot)) )
          enddo
	    endif
  endif
! ---------------------------------------------------------------------------

  allocate(fv1(7,kinde,nre,nv,nalongQ1) , fv1t(7,kinde,nre,nv,nalongQ1))
  fv1   = 0.

  do mn=1,kinds
     if(flagi)    particle = eles(mn)  
     if(flage) particle = ions(mn)  
	 amass    = particle%amass
	 frac     = particle%frac
	 scale    = particle%vth
	 if((ind(mn+10) ==0 .and. flage) .or. (ind(mn) ==0 .and. flagi)) goto 220
     do m=1,mblocks
           do L = 1,particle%block(m)%mi

                 p      = particle%block(m)%qv(1,L)
                 q      = particle%block(m)%qv(2,L)
                 w      = particle%block(m)%qv(3,L)
			     weight   = particle%block(m)%qv(13,L) * particle%block(m)%qv(14,L)  * frac

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1

                 vb    = bv0_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv0_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv0_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv0_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv0_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv0_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv0_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv0_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 vb1   = bv1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bv1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bv1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bv1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bv1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bv1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bv1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bv1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  
                 ve1   = ev1_(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             ev1_(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 ev1_(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 ev1_(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 ev1_(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             ev1_(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 ev1_(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 ev1_(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
					 B_   = sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
!		             xyz  = pqw2xyz(0,pmin+p*deltax,qequator,wmin+w*deltaz) !get magnetic field at equator
!                     x    = xyz(1)
!                     y    = xyz(2)
!                     z    = xyz(3)
!			         vb_eq= dipole3DB(x,y,z)
!					 B_eq = sqrt(vb_eq(1)**2+vb_eq(2)**2+vb_eq(3)**2)
				 endif

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb       = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2) 
                 vbper1   = a_cross_b(vb,yaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,xaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 
                 x2b(1,:) = vb
                 x2b(2,:) = vbper1
                 x2b(3,:) = vbper2

				 uu       = particle%block(m)%qv(4:6,L)/amass
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 ee       = amass*cspeed**2*(gamma -1.)*ut0

				 vv       = uu /gamma     ! scale factor here
				 vv       = vv /scale 

				 vw(1)    = x2b(1,1)*uu(1) + x2b(1,2)*uu(2)+x2b(1,3)*uu(3)
				 vw(2)    = x2b(2,1)*uu(1) + x2b(2,2)*uu(2)+x2b(2,3)*uu(3)
				 vw(3)    = x2b(3,1)*uu(1) + x2b(3,2)*uu(2)+x2b(3,3)*uu(3)
				 mu       = (vw(2)**2+vw(3)**2)/2./B_

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)
                 angle    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 
                 angle    = mod(angle + pi2*2,  pi2)    ! southward ????  --- csi ----


				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)
				 ew(1)    = x2b(1,1)*ve1(1) + x2b(1,2)*ve1(2)+x2b(1,3)*ve1(3)
				 ew(2)    = x2b(2,1)*ve1(1) + x2b(2,2)*ve1(2)+x2b(2,3)*ve1(3)
				 ew(3)    = x2b(3,1)*ve1(1) + x2b(3,2)*ve1(2)+x2b(3,3)*ve1(3)
				 vpar_    = vw(1)
				 vper_    = sqrt(vw(2)**2+vw(3)**2)

! ----------------------------------------------------------------------------------
                 vde1_    = vw(1) * ew(1)
                 vde2_    = vw(2) * ew(2) + vw(3) * ew(3)
                 vdb_     = a_dot_b(vv, vb1) /max(sqrt(vb1(1)**2+vb1(2)**2+vb1(3)**2),1.e-8_p2)

				 qflux(1:7)    = (/vde1_,vde2_,vdb_, angle, ee, mu, 1._p2/) * weight

                 if(vpar_ < vminp .or. vpar_ >=vmaxp ) goto 221 
                 do n = 1, nre
!				    if(abs(p-pval(n)) <= dpval(n)) then
                       imc1  = min(max(int(max(vpar_ - vminp,0.)/dvpar + 1.),1),nv-1)
		               dmc1  = imc1-(vpar_ - vminp)/dvpar
                       imc2  = min(int(p/dQ + 1.),nalongQ1-1)
		               dmc2  = imc2 - p/dQ


                       fv1(1:7,mn,n,imc1,imc2)     = fv1(1:7,mn,n,imc1,imc2)    + dmc1*dmc2*qflux(1:7)	          
                       fv1(1:7,mn,n,imc1,imc2+1)   = fv1(1:7,mn,n,imc1,imc2+1)  + dmc1*(1.-dmc2)*qflux(1:7)		          
                       fv1(1:7,mn,n,imc1+1,imc2)   = fv1(1:7,mn,n,imc1+1,imc2)  + (1.-dmc1)*dmc2*qflux(1:7)		          
                       fv1(1:7,mn,n,imc1+1,imc2+1) = fv1(1:7,mn,n,imc1+1,imc2+1)+ (1.-dmc1)*(1.-dmc2)*qflux(1:7)	          


!					endif
                 enddo

221           continue
		  enddo
     enddo
220  continue
  enddo

  CALL MPI_ALLREDUCE(fv1,fv1t,size(fv1),MPI_real,MPI_SUM,mpi_comm_world,IERROR)

  allocate(den(kinde,nre,nv,nalongQ1))
  
  den   = max(fv1t(7,:,:,:,:),1.e-12)  !/iftestParticle(6)


  if(mype==0) then
     allocate(plot(kinde,nre,nv,nalongQ1))
     do iplot = 1, 6
        plot  = fv1t(iplot,:,:,:,:)/den
        call check( nf90_put_var(ncid, varid(10+iplot), plot, start = istart1, &
                              count = icount1) )
     enddo
        iplot = 7
        plot  = den
        call check( nf90_put_var(ncid, varid(10+iplot), plot, start = istart1, &
                              count = icount1) )

     deallocate(plot)
  endif



  if(mype==0) call check( nf90_close(ncid) )


  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  call allocate_bvector(ev1_,-3)
  deallocate(den,fv1,fv1t)
  deallocate(vpar,Qarray)
  
  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine output_4fvparL_ncfile_case58



!=================================================================================================
subroutine diagnos_particle_nc
  use netcdf
  use global_parameters
  use grid_data
  use vector_functions
  use diagnos
  implicit none
  include 'mpif.h'
  integer :: nions
  real,pointer,dimension(:,:) :: send,recv
  integer,pointer,dimension(:) :: num_recv,request
  type(bvector_type), dimension(:), pointer    :: Beff
  real(p2),pointer,dimension(:,:) :: qvp_writi,qvp_write,dummy
  integer :: num_diag,num_particle_diag,num_particle_diag_pe,num,num0,numt, &
             num_skip,m,L,icount,ierror
  integer :: nregion,n,kind,loop,Li_tot,Le_tot,Lt,i,j,k,ind,mkind,mn
  real(p2) :: r,theta,zeta,vx,vy,vz,x,y,z,vper,amass,charge,frac,p,q,w
  real(p2) :: vb(3),vv(3),vbper1(3),vbper2(3),vpar,vper1,vper2,&
               W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,angle
  integer  :: ii,jj,kk,ip,jp,kp,ier,nbuff,nbuffi,nbuffe

  integer,pointer,dimension(:)   :: L_start,L_ok
  real,pointer,dimension(:,:) :: plot
  real,pointer,dimension(:,:) :: bb 
  real     :: bs(3),bsa(3)
  real  :: param(180)
  character*6 fldname(74)
  character*4 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(4),varid(30),dimids(2)
  logical :: yes
  integer :: status(mpi_status_size)
  real(8) :: xyz(3)


  interface
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine bconcarB(a,b)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: a,b
      end subroutine bconcarB

     function pqw2xyz(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: pqw2xyz(3)
		integer :: ml
	 end function pqw2xyz
     function xyz2pqw(ml,a,b,c)
        use constants
	    real*8 :: a,b,c
	    real*8 :: xyz2pqw(3)
		integer :: ml
	 end function xyz2pqw
      function dipole3DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole3DB(3)
	  end function dipole3DB

  end interface
!------------------------
  ind      = 11
  num_diag = max(int(pdomain(10)),4000)
  nregion  = int(pdomain(1)+0.001)
  if(n_pregion <1) return
  allocate(bb(n_pregion,3))

  call allocate_bvector(beff,3)
  do m=1,mblocks
    Beff(m)%vector = bv0(m)%vector + bv1(m)%vector
  enddo
  call bconcarB(beff,Bcar)
  call allocate_bvector(beff,-3)

    do m=1,n_pregion
       i = min(max(int(domain_pregion(m,1)+domain_pregion(m,2))/2,1),nx)
       j = min(max(int(domain_pregion(m,3)+domain_pregion(m,4))/2,1),nx)
       k = min(max(int(domain_pregion(m,5)+domain_pregion(m,6))/2,1),nx)
       bs = 0.
	   do n=1,mblocks
          if(i>=block(n)%i0 .and. i<=block(n)%i1 .and. &
		     j>=block(n)%j0 .and. j<=block(n)%j1 .and. &
			 k>=block(n)%k0 .and. k<=block(n)%k1) bs = bcar(n)%vector(:,i,j,k)
	   enddo

       CALL MPI_ALLREDUCE(bs,bsa,3,mpi_real,MPI_SUM,myComm1,IERROR)
       bb(m,:) = bsa
    enddo

  if(mype==0) then 
     nbuffi = 50000
     nbuffe = 50000
  	 allocate(QVP_writi(ind,nbuffi)) !nbuff
  	 allocate(QVP_write(ind,nbuffe)) !nbuff
  endif


! do ions first
  yes     = .false.
  do i = 1,int(pdomain(3))
     if(pdomain(3+i) < 10) yes = .true.
  enddo

!  ---------- do ions -----------------------------------
!每个粒子11个信息
!mn：粒子种类，整数
!p,q,w：位置坐标，单位是格点数
!vx,vy,vz：三个速度：模拟中速度的单位
!vpar：平行速度，平行于磁场方向的速度
!vper1、vper2：垂直方向速度，汪老师知道垂直方向怎么选的
!frac：粒子的权重，汪老师知道 ions(mn)%block%...取默认值1



  if(yes) then
     mkind = size(ions)

     nbuff = max(maxval(ions(1)%block%mi),2000)*2
     allocate(send(11, nbuff),num_recv(numberpe),request(numberpe))

     num =0
	 mkind = size(ions)

	 do mn = 1, mkind
        amass  = ions(mn)%amass
        charge = ions(mn)%charge
        frac   = ions(mn)%frac

     do m=1,mblocks
        do L = 1,ions(mn)%block(m)%mi
           p     = ions(mn)%block(m)%qv(1,L)
           q     = ions(mn)%block(m)%qv(2,L)
           w     = ions(mn)%block(m)%qv(3,L)

           vx     = ions(mn)%block(m)%qv(4,L)/amass
           vy     = ions(mn)%block(m)%qv(5,L)/amass
           vz     = ions(mn)%block(m)%qv(6,L)/amass

           do n = 1,nregion
			  if(p>= domain_pregion(n,1) .and. p<= domain_pregion(n,2) .and. & 
			     q>= domain_pregion(n,3) .and. q<= domain_pregion(n,4) .and. & 
			     w>= domain_pregion(n,5) .and. w<= domain_pregion(n,6) )then
                 num = num +1

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1
                 vb    = bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1 
				 if(coordinate==0) vb =(/br0,bt0,bz0/)		  

		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vbper1  = a_cross_b(vb,xaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,zaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 


				 vv      = ions(mn)%block(m)%qv(4:6,L)/amass
				 vpar    = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)
				 vper1   = vbper1(1)*vv(1) + vbper1(2)*vv(2) + vbper1(3)*vv(3)
				 vper2   = vbper2(1)*vv(1) + vbper2(2)*vv(2) + vbper2(3)*vv(3)


				 if(num >=  nbuff) then
                    allocate (dummy(ind,num-1))
					dummy   = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(ind,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
                 endif
				 send(1,num) = mn
				 send(2,num) = p
				 send(3,num) = q
				 send(4,num) = w
				 send(5,num) = vx
				 send(6,num) = vy
				 send(7,num) = vz
				 send(8,num) = vpar
				 send(9,num) = vper1
				 send(10,num) = vper2
				 send(11,num) = ions(mn)%block(m)%qv(13,L) * ions(mn)%block(m)%qv(14,L) * frac
				 goto 110
		       endif  
			enddo  !end of n --region
110			continue
	    enddo ! end of L
     enddo ! end of mblock
	 enddo 
!    ------------------------- reduce number of particles here ---------------
     call MPI_ALLREDUCE(num,numt,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IER)
	 num0     = num
	 num_skip = max(numt/num_diag,1)
	 num      = 0
	 do m=1,num0,num_skip
        num = num+1
		send(:,num) = send(:,m)
	 enddo
!    --------------------------------------------------------
     call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)

     if(mype > 0 .and. num > 0) then
	       call mpi_isend(send(:,1:num),ind*num, mpi_real,0, &
		                           mype, MPI_COMM_WORLD, request(mype),ier)

     endif
     if(mype==0) then
	    Li_tot = 0
        do i=1,num
		   Li_tot = Li_tot +1
		   if(Li_tot >= nbuffi) then
              allocate(dummy(ind,Li_tot-1))
			  dummy = QVP_writi(:,1:Li_tot-1)
			  deallocate(QVP_writi)
			  nbuffi = nbuffi*1.5
			  allocate(QVP_writi(ind,nbuffi))
			  QVP_writi(:,1:Li_tot-1) = dummy
			  deallocate(dummy)
		   endif
           QVP_writi(:,Li_tot) = send(:,i)
        enddo
     endif

     do L=1,numberpe-1
	    if(num_recv(L+1) > 0) then
	       if(mype==0) then
		      allocate(recv(ind,num_recv(L+1)))
              call mpi_recv(recv,ind*num_recv(L+1), mpi_real,L, &
		                           L, MPI_COMM_WORLD, status, ier)
		  
		      do i=1,num_recv(L+1)
		         Li_tot = Li_tot +1
		         if(Li_tot >= nbuffi) then
                    allocate(dummy(ind,Li_tot-1))
			        dummy = QVP_writi(:,1:Li_tot-1)
			        deallocate(QVP_writi)
			        nbuffi = nbuffi*1.5
			        allocate(QVP_writi(ind,nbuffi))
			        QVP_writi(:,1:Li_tot-1) = dummy
			        deallocate(dummy)
		         endif
                 QVP_writi(:,Li_tot) = recv(:,i)
              enddo
			  deallocate(recv)
           elseif(mype==L)then
		      call MPI_Wait(request(L), status, ier)
		   endif
		endif
	 enddo 
	 deallocate(send,request,num_recv)
  endif

! now do electrons
! num_skip = max(1,me/NUM_PARTICLE_DIAG_PE)
! if(qsload == 1) num_skip = 1


  yes     = .false.
  do i = 1,int(pdomain(3))
     if(pdomain(3+i) > 10) yes = .true.
  enddo
! ------------------ do electrons -----------------------------------
  if(yes) then
	 mkind = size(eles)
     nbuff = max(maxval(eles(1)%block%mi),2000)*2
     allocate(send(ind, nbuff),num_recv(numberpe),request(numberpe))
     num =0

	 do mn=1,mkind
        amass  = eles(mn)%amass
        charge = eles(mn)%charge
        frac   = eles(mn)%frac

     do m=1,mblocks
        do L = 1,eles(mn)%block(m)%mi
           p     = eles(mn)%block(m)%qv(1,L)
           q     = eles(mn)%block(m)%qv(2,L)
           w     = eles(mn)%block(m)%qv(3,L)
           vx    = eles(mn)%block(m)%qv(4,L)/amass
           vy    = eles(mn)%block(m)%qv(5,L)/amass
           vz    = eles(mn)%block(m)%qv(6,L)/amass
           do n = 1,nregion
			  if(p>= domain_pregion(n,1) .and. p<= domain_pregion(n,2) .and. & 
			     q>= domain_pregion(n,3) .and. q<= domain_pregion(n,4) .and. & 
			     w>= domain_pregion(n,5) .and. w<= domain_pregion(n,6) )then
                 num = num +1

			     wz1      = w +1.0
                 kk       = max(block(m)%nz0,min(block(m)%nz1-1,int(wz1)))
	             kp       = min(block(m)%nz1,kk+1)
			     wz1      = wz1 - kk
	             wz0      = 1.0 - wz1

			     wx1      = p +1.0
                 ii       = max(block(m)%nx0,min(block(m)%nx1-1,int(wx1)))
	             ip       = min(block(m)%nx1,ii+1)
			     wx1      = wx1 - ii
	             wx0      = 1.0 - wx1

			     wy1      = q +1.0
                 jj       = max(block(m)%ny0,min(block(m)%ny1-1,int(wy1)))
	             jp       = min(block(m)%ny1,jj+1)
			     wy1      = wy1 - jj
	             wy0      = 1.0 - wy1
                 vb    = bcar(m)%vector(:,ii,jj,kk) * wx0*wy0*wz0  + &
			             bcar(m)%vector(:,ip,jj,kk) * wx1*wy0*wz0  + &
						 bcar(m)%vector(:,ii,jp,kk) * wx0*wy1*wz0  + &
						 bcar(m)%vector(:,ip,jp,kk) * wx1*wy1*wz0  + &
						 bcar(m)%vector(:,ii,jj,kp) * wx0*wy0*wz1  + &
			             bcar(m)%vector(:,ip,jj,kp) * wx1*wy0*wz1  + &
						 bcar(m)%vector(:,ii,jp,kp) * wx0*wy1*wz1  + &
						 bcar(m)%vector(:,ip,jp,kp) * wx1*wy1*wz1  

                 if(dipole) then
		             xyz  = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz) !get magnetic field at equator
                     x    = xyz(1)
                     y    = xyz(2)
                     z    = xyz(3)
			         vb   = dipole3DB(x,y,z)
				 endif


		         if(sqrt(vb(1)**2+vb(2)**2+vb(3)**2)<1.e-5) vb = zaix
                 vb      = vb/sqrt(vb(1)**2+vb(2)**2+vb(3)**2)
                 vbper1  = a_cross_b(vb,xaix)
		         if(sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) ==0.)vbper1  = a_cross_b(vb,zaix)
                 vbper2   = a_cross_b(vb,vbper1)
                 vbper1   = vbper1/sqrt(vbper1(1)**2+vbper1(2)**2+vbper1(3)**2) 
                 vbper2   = vbper2/sqrt(vbper2(1)**2+vbper2(2)**2+vbper2(3)**2) 


				 vv      = eles(mn)%block(m)%qv(4:6,L)/amass
				 vpar    = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)
				 vper1   = vbper1(1)*vv(1) + vbper1(2)*vv(2) + vbper1(3)*vv(3)
				 vper2   = vbper2(1)*vv(1) + vbper2(2)*vv(2) + vbper2(3)*vv(3)

				 if(num >=  nbuff) then
                    allocate (dummy(ind,num-1))
					dummy   = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(ind,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
                 endif
				 send(1,num) = block(m)%ele(L)%kind
				 send(2,num) = p
				 send(3,num) = q
				 send(4,num) = w
				 send(5,num) = eles(mn)%block(m)%qv(4,L)/amass
				 send(6,num) = eles(mn)%block(m)%qv(5,L)/amass
				 send(7,num) = eles(mn)%block(m)%qv(6,L)/amass
				 send(8,num) = vpar
				 send(9,num) = vper1
				 send(10,num) = vper2
				 send(11,num) = eles(mn)%block(m)%qv(13,L) * eles(mn)%block(m)%qv(14,L) * frac
				 goto 111
		       endif  
			enddo  !end of n --region
111			continue
	    enddo ! end of L
     enddo ! end of mblock
	 enddo


!    ------------------------- reduce number of particles here ---------------
     call MPI_ALLREDUCE(num,numt,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IER)
	 num0     = num
	 num_skip = max(numt/num_diag,1)
	 num      = 0
	 do m=1,num0,num_skip
        num = num+1
		send(:,num) = send(:,m)
	 enddo
!    --------------------------------------------------



     call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)
     if(mype > 0 .and. num > 0) then
	       call mpi_isend(send(:,1:num),ind*num, mpi_real,0, &
		                           mype, MPI_COMM_WORLD, request(mype),ier)

     endif

     if(mype==0) then
	    Le_tot = 0
        do i=1,num
		   Le_tot = Le_tot +1
		   if(Le_tot >= nbuffe) then
              allocate(dummy(ind,Le_tot-1))
			  dummy = QVP_write(:,1:Le_tot-1)
			  deallocate(QVP_write)
			  nbuffe = nbuffe*1.5
			  allocate(QVP_write(ind,nbuffe))
			  QVP_write(:,1:Le_tot-1) = dummy
			  deallocate(dummy)
		   endif
           QVP_write(:,Le_tot) = send(:,i)
        enddo
     endif

     do L=1,numberpe-1
	    if(num_recv(L+1) > 0) then
	       if(mype==0) then
		      allocate(recv(ind,num_recv(L+1)))
              call mpi_recv(recv,ind*num_recv(L+1), mpi_real,L, &
		                           L, MPI_COMM_WORLD, status, ier)
		  
		      do i=1,num_recv(L+1)
		         Le_tot = Le_tot +1
		         if(Le_tot >= nbuffe) then
                    allocate(dummy(ind,Le_tot-1))
			        dummy = QVP_write(:,1:Le_tot-1)
			        deallocate(QVP_write)
			        nbuffe = nbuffe*1.5
			        allocate(QVP_write(ind,nbuffe))
			        QVP_write(:,1:Le_tot-1) = dummy
			        deallocate(dummy)
		         endif
                 QVP_write(:,Le_tot) = recv(:,i)
              enddo
			  deallocate(recv)
           elseif(mype==L)then
		      call MPI_Wait(request(L), status, ier)
		   endif
		endif
	 enddo 
	 deallocate(send,request,num_recv)
  endif

! ---------------------------------------------------------------
  param      = 0.
  param(1)   = stime   ! has been scaled to omega_i
  param(2)   = Li_tot  ! has been scaled to omega_i
  param(3)   = Le_tot  ! has been scaled to omega_i

  param(4)   = ind
  param(5)   = kinds
  param(6)   = betae
  param(7)   = nx
  param(8)   = ny
  param(9)   = nz

  param(10:kinds+9)   = fraci(1:kinds)
  param(20:kinds+19)  = frace(1:kinds)
  param(30)   = nregion
  param(31:30+nregion)   = domain_pregion(1:nregion,1)
  param(41:40+nregion)   = domain_pregion(1:nregion,2)
  param(51:50+nregion)   = domain_pregion(1:nregion,3)
  param(61:60+nregion)   = domain_pregion(1:nregion,4)
  param(71:70+nregion)   = domain_pregion(1:nregion,5)
  param(81:80+nregion)   = domain_pregion(1:nregion,6)
  param(91:90+nregion)   = bb(1:nregion,1)
  param(101:100+nregion) = bb(1:nregion,2)
  param(111:110+nregion) = bb(1:nregion,3)


  param(121)  = ulength*deltax
  param(122)  = ulength*deltay
  param(123)  = ulength*deltaz
  param(124)  = uspeed
  param(125)  = Va

! ------------------------------
  idiagP = idiagP + 1
  write(ct,'(I4.4)')IDIAGP 
  if(Li_tot+Le_tot == 0 ) goto 831


  if(mype==0) then
        filename = 'particle'//ct//'.nc'
        call check( nf90_create(FILENAME, nf90_clobber, ncid) )

        call check( nf90_def_dim(ncid, 'param', 180, dimid(1)) )
        call check( nf90_def_var(ncid, 'param', NF90_REAL, dimid(1), varid(1)) )
        call check( nf90_put_att(ncid,  varid(1), UNITS, 'non') )

        call check( nf90_def_dim(ncid, 'var', ind, dimid(2)) )
        call check( nf90_def_var(ncid, 'var', NF90_REAL, dimid(2), varid(2)) )
        call check( nf90_put_att(ncid,  varid(2), UNITS, 'non') )

        call check( nf90_def_dim(ncid, 'nion', Li_tot, dimid(3)) )
        call check( nf90_def_var(ncid, 'nion', NF90_REAL, dimid(3), varid(3)) )
        call check( nf90_put_att(ncid,  varid(3), UNITS, 'non') )

        call check( nf90_def_dim(ncid, 'nele', Le_tot, dimid(4)) )
        call check( nf90_def_var(ncid, 'nele', NF90_REAL, dimid(4), varid(4)) )
        call check( nf90_put_att(ncid,  varid(4), UNITS, 'non') )


        dimids=(/dimid(2),dimid(3)/)
        call check( nf90_def_var(ncid, 'ion', NF90_REAL, dimids, varid(5)) )
        call check( nf90_put_att(ncid, varid(5), UNITS, 'non') )
        dimids=(/dimid(2),dimid(4)/)
        call check( nf90_def_var(ncid, 'ele', NF90_REAL, dimids, varid(6)) )
        call check( nf90_put_att(ncid, varid(6), UNITS, 'non') )


        call check( nf90_enddef(ncid) )

        call check( nf90_put_var(ncid, varid(1), param) )
		allocate(plot(ind,Li_tot))
		plot = QVP_writi(:,1:Li_tot)
        call check( nf90_put_var(ncid, varid(5), plot)  )

        deallocate(plot)
		allocate(plot(ind,Le_tot))
		plot = QVP_write(:,1:Le_tot)

        call check( nf90_put_var(ncid, varid(6), plot) )
		deallocate(plot)

        call check( nf90_close(ncid) )
  endif
831  continue
   
  if(mype ==0 ) deallocate(qvp_writi,qvp_write)
  deallocate(bb)

  return

  contains
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
end subroutine diagnos_particle_nc


