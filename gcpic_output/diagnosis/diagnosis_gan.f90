subroutine diagnos_field(box9,iflag)
  use global_parameters
  use grid_data
  use diagnos
  use global_case
  implicit none
  include 'mpif.h'
  type(diag_type)   :: box9

  real,    pointer,dimension(:,:,:) :: plot
  real(p2) :: param(80)
  integer :: ic,i,iplot,ifield,m,ifchecki,ifchecke,n,ifcheck(30),iflag

  type(bscalar_type), dimension(:), pointer    :: temp,beta_i,beta_e
  type(bvector_type), dimension(:), pointer    :: bfld,efld,ui,ue,delJ,dumm,av1,delB,delE

  CHARACTER*5 CT
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
    subroutine gather_vector(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
       real(p2),pointer,dimension(:,:,:,:) :: b
    end subroutine gather_vector
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

    subroutine bbcellecenter(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bbcellecenter

	  subroutine bounds(a,char)
        use constants
		use grid_data
        type(bvector_type), dimension(:), pointer    :: a
	    character*1 char
	  end subroutine bounds
    subroutine bconvert_rtheta(a)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
    end subroutine bconvert_rtheta

      subroutine Adconvert(a,b,c)
		use grid_data
	    type(bscalar_type), dimension(:), pointer    :: a
        real,    pointer,dimension(:,:,:)            :: b
	    type(diag_type)                              :: c
      end subroutine Adconvert

    subroutine getUT(m,a)
       use constants
	   use grid_data
	   type(particle_type), dimension(:), pointer    :: a
    end subroutine getUT


  end interface
!------------------------
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

  if( ghybrid) write(ct,'(I5.5)')int(stime + 2.*dt)


! do the first time of coefficient
!  if(ic == 1) then
!     allocate(diagweight(6,nxd*nyd*nzd),LL(nxd,nyd,nzd))
!	 allocate(xdiag(nxd),ydiag(nyd),zdiag(nzd))
!	 call stconvert
!  endif
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
  if(mod(simtype(2),100) == 51) param(20) = 1./lamda_de
  param(21) = cva
  param(22) = lamda_i/rhos
  param(23) = lamda_e/rhos
  param(24) = uRe/rhos
  param(25) = lamda_de/rhos
  param(29) = va
  param(30) = alpha

  if(mod(simtype(2),100) >= 90) param(21:28) = bnvt_sw

! some key parameters
  param(41) = omega_lh/omega_i
  param(75-box9%nfield+1:75)= box9%DiagQuantity(1:box9%nfield) 
  param(77) = simtype(2)
  param(78) = simtype(1)
  param(79) = coordinate
  param(80) = box9%nfield   ! number of field quantity are recorded


  IF(mype==0) THEN
!     if(box9%idiagbox ==1)then
	    if(iflag == 0) open(11, FILE = 'field'//ct//'.dat')      !, FORM="unformatted")
	    if(iflag == 1) open(11, FILE = 'fieldfine'//ct//'.dat')      !, FORM="unformatted")
	    if(iflag == 2) open(11, FILE = 'fieldbox'//ct//'.dat')      !, FORM="unformatted")
        write(11,*)3,80
        write(11,*)box9%NXD,box9%NYD,box9%NZD
        write(11,*)param
        write(11,1102)box9%xdiag,box9%ydiag,box9%zdiag
!     else
!	    open(11, FILE = 'field'//ct//'.dat')      !, FORM="unformatted")
!           write(11,*)3,80
!           write(11,*)cuv%NXD,cuv%NYD,cuv%NZD
!           write(11,*)param
!           write(11,1102)cuv%xdiag,cuv%ydiag,cuv%zdiag

!	 endif
  ENDIF
!--------------------------------------------------------
  allocate(plot(box9%nxd,box9%nyd,box9%nzd))
! if(idiagbox ==0) allocate(plot(cuv%nxd,cuv%nyd,cuv%nzd))

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
  if(box9%idiagbox ==1 .or. box9%idiagbox ==3) then
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
  if(box9%idiagbox ==1) then
     do m=1,mblocks
        dumm(m)%vector = ev1(m)%vector 
     enddo
	 call bcovcarE(dumm,efld)
	 if(mod(simtype(2),100)==51) call bounds(efld,'e')
  else
     do m=1,mblocks
	    do n=1,3
           dumm(m)%vector(n,:,:,:) =  ev1(m)%vector(n,:,:,:)/block(m)%node%h(3,n,:,:,:)   
        enddo
     enddo
     call becellecenter(dumm,efld)
	 if(mod(simtype(2),100)==51) call bounds(efld,'e')
  endif
! (3) ui-field -------
  if(box9%idiagbox ==1) then
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
     if(kinds_fi > 0) then
     do m=1,mblocks
	    ue(m)%vector(1,:,:,:) = block(m)%ionf(1,:,:,:)%v(1)
	    ue(m)%vector(2,:,:,:) = block(m)%ionf(1,:,:,:)%v(2)
	    ue(m)%vector(3,:,:,:) = block(m)%ionf(1,:,:,:)%v(3)
	 enddo
	 endif
  else
     if(box9%idiagbox ==1) then
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
	 if((ifield >= 40 .and. ifield <= 58) .or. (ifield >= 87 .and. ifield<=92)  ) then
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
		   if(mod(simtype(1),10) >=6) then   ! hybrid simulation j= curlB
  	          do m=1,mblocks
		         dumm(m)%vector    = bv1(m)%vector  + bv0(m)%vector 
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

	 if(ifield  == 7 .or. ifield  == 8 .or. ifield  == 9 ) then  ! A to plot magnetic field lines
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

              if(Ldevice) call bconvert_rtheta(delB)

		   else
  	          call field3(delB)
           endif
        endif
	 endif

	 if(ifield  == 27 .or. ifield  == 28 .or. ifield  == 29 ) then  ! delta_E
	    ifcheck(6) = ifcheck(6) +1
		if(ifcheck(6) ==1) then
           if(abs(coordinate)  ==1 .or. abs(coordinate)  ==2 .or. abs(coordinate)  ==3 .or. &
		      abs(coordinate)  ==8 .or. abs(coordinate)  ==5 .or. Ldevice) then
		      do m=1,mblocks
	             do n =1,3
                    dumm(m)%vector(n,:,:,:)  = ev1(m)%vector(n,:,:,:) / block(m)%node%h(3,n,:,:,:)
		         enddo
			  enddo
	          call becellecenter(dumm,delE)
              if(Ldevice) call bconvert_rtheta(delE)
		   else

           endif
        endif
	 endif


	 if(ifield  == 71) then  ! beta_i
        do m=1,mblocks
		   beta_i(m)%scalar = 0.
		   do n=1,kinds
		      beta_i(m)%scalar = beta_i(m)%scalar + block(m)%ini(n,:,:,:)%ni0 * &
			                     (block(m)%ini(n,:,:,:)%Ti0(1)+block(m)%ini(n,:,:,:)%Ti0(2)*2)/3.
		   enddo
		   temp(m)%scalar = bfld(m)%vector(1,:,:,:)**2+bfld(m)%vector(2,:,:,:)**2+bfld(m)%vector(3,:,:,:)**2
        enddo
        do m=1,mblocks
		   beta_i(m)%scalar = beta_i(m)%scalar/max(temp(m)%scalar,1.e-4_p2)*alpha
        enddo
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
		      temp(m)%scalar = ne(m)%scalar+nef(m)%scalar
		      if(simtype(1) == 6 .or. simtype(1) == 7 .or. simtype(1) == 66 .or. simtype(1) == 67) &
              temp(m)%scalar = nif(m)%scalar
		   enddo
        case(6)                               !phi
		   do m=1,mblocks
		      temp(m)%scalar = phi(m)%scalar
		   enddo
        case(7,8,9)                           !a1
		   do m=1,mblocks
		      temp(m)%scalar = av1(m)%vector(ifield-6,:,:,:)
		   enddo
        case(10)                              !Ue||--- Je||


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
        case(21,22,23)                        !delta J
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
!        case(27,28,29)                        !dj X B0
!          allocate(delJ(3,nz0:nz1,mgrid0:mgrid1))
!		   call get_j(bv1,delJ)
!		   call cross(delJ,Bv_bar,vec1)
!           temp = vec1(ifield-26,:,:)
!	       call dconvert(temp,plot)
!		   deallocate(delJ)	 
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
           if(mype==0) write(11,1102)PLOT	 
        case(39)                        !divPpe/ne
!           temp = dUedt(3,:,:)
! ----------------------------------------------plasma informations ------------------
        case(41,42,43)        !ui(1-3,:)
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%vi0(ifield-40)
		   enddo
        case(44)              !Ni
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ni0
		   enddo
        case(45,46)           !ui(1-3,:),Ti
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ti0(ifield-44)
		   enddo
        case(47,48,49)        !ui(1-3,:) kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%vi0(ifield-46)
		   enddo
        case(50)              !ui(1-3,:),Ni, kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ni0
		   enddo
        case(51,52)           !ui(1-3,:),Ti, kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ti0(ifield-50)
		   enddo

        case(53,54,55)        !ui(1-3,:) kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%vi0(ifield-52)
		   enddo
        case(56)              !ui(1-3,:),Ni, kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%ni0
		   enddo
        case(57,58)           !ui(1-3,:),Ti, kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%ti0(ifield-56)
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
        case(81,82)              ! for electron Three kinds 81-86
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%edotv(ifield-80)
		   enddo

        case(83,84)              
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%edotv(ifield-82)
		   enddo

        case(85,86)              ! for electron Three kinds 81-86
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%edotv(ifield-84)
		   enddo

        case(87,88)              ! for electron Three kinds 81-86
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%edotv(ifield-86)
		   enddo

        case(89,90)              
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%edotv(ifield-88)
		   enddo

        case(91,92)              ! for electron Three kinds 81-86
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%edotv(ifield-90)
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
	       
!     call dconvert(temp,plot,box9)
!	 if(idiagbox ==0) call dconvert(temp,plot,cuv)
     if(mype==0) write(11,1102)PLOT	 
  enddo
  if(mype==0) close(11)
! 1102  format(5f14.5)
  1102  format(5e16.5)
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
  return
end subroutine diagnos_field

!============================================================================================

subroutine dconvert(c,dall,cc)
  use global_parameters
  use grid_data
  implicit none
  include 'mpif.h'
  type(bscalar_type), dimension(:), pointer    :: c
  type(diag_type)   :: cc
  real,    pointer,dimension(:,:,:) :: d,dall 
  integer :: ii,jj,kk,ip,jp,kp,L,i,j,k,ierror,m
  real(p2) :: wx0,wx1,wy0,wy1,wz0,wz1
!-----------------
  if(kmype /= 0) return  
  allocate(d(cc%nxd,cc%nyd,cc%nzd))

  d   = 0.;  
	  L=0;  
      do I = 1, cc%nxd
	     do J = 1, cc%nyd
	        do K = 1, cc%nzd
			   if(cc%LL(i,j,k) >0 ) Then
			     L      = L+1;
				 m      = cc%LL(i,j,k)
                 ii     = int(cc%diagweight(1,L))
                 jj     = int(cc%diagweight(3,L))
                 kk     = int(cc%diagweight(5,L))
				 ip     = min(ii+1,block(m)%nx1)
				 jp     = min(jj+1,block(m)%ny1)
				 kp     = min(kk+1,block(m)%nz1)
                 wx1    = cc%diagweight(2,L)
                 wy1    = cc%diagweight(4,L)
                 wz1    = cc%diagweight(6,L)
                 wx0    = 1.0-wx1
                 wy0    = 1.0-wy1
                 wz0    = 1.0-wz1



                 d(i,j,k)   = c(m)%scalar(ii,jj,kk) * wx0*wy0*wz0 + &
                              c(m)%scalar(ip,jj,kk) * wx1*wy0*wz0 + &
                              c(m)%scalar(ii,jp,kk) * wx0*wy1*wz0 + &
                              c(m)%scalar(ip,jp,kk) * wx1*wy1*wz0 + &
							  c(m)%scalar(ii,jj,kp) * wx0*wy0*wz1 + &
                              c(m)%scalar(ip,jj,kp) * wx1*wy0*wz1 + &
                              c(m)%scalar(ii,jp,kp) * wx0*wy1*wz1 + &
                              c(m)%scalar(ip,jp,kp) * wx1*wy1*wz1
                endif
             enddo
          enddo
       enddo
! if(p2==8)CALL MPI_ALLREDUCE(d,dall,size(d),mpi_double_precision,MPI_SUM,myComm1,IERROR)
  CALL MPI_ALLREDUCE(d,dall,size(d),mpi_real,MPI_SUM,myComm1,IERROR)
  DEALLOCATE(d)
 
  return
end subroutine dconvert

!============================================================================================

subroutine Adconvert(c,dall,cc)
  use global_parameters
  use grid_data
  implicit none
  include 'mpif.h'
  type(bscalar_type), dimension(:), pointer    :: c
  type(diag_type)   :: cc
  real,    pointer,dimension(:,:,:) :: d,dall 
  integer :: ii,jj,kk,ip,jp,kp,L,i,j,k,ierror,m
  real(p2) :: wx0,wx1,wy0,wy1,wz0,wz1
!-----------------
  if(kmype /= 0) return  
  allocate(d(cc%nxd,cc%nyd,cc%nzd))

  d   = 0.;  
  do m=1,mblocks1
     do i=block(m)%i0,block(m)%i1p
        do j=block(m)%j0,block(m)%j1p
	       do k=block(m)%k0,block(m)%k1p
              if(cc%adiag(m)%vector(7,i,j,k)> 0) then
			  ii  = cc%adiag(m)%vector(1,i,j,k) 
			  jj  = cc%adiag(m)%vector(2,i,j,k) 
			  kk  = cc%adiag(m)%vector(3,i,j,k) 
			  ip  = min(cc%nxd,ii+1)
			  jp  = min(cc%nyd,jj+1)
			  kp  = min(cc%nzd,kk+1)
			  wx1 = cc%adiag(m)%vector(4,i,j,k) 
			  wy1 = cc%adiag(m)%vector(5,i,j,k) 
			  wz1 = cc%adiag(m)%vector(6,i,j,k) 
              wx0    = 1.0-wx1
              wy0    = 1.0-wy1
              wz0    = 1.0-wz1
              d(ii,jj,kk)   = d(ii,jj,kk) + c(m)%scalar(i,j,k) * wx0*wy0*wz0
              d(ii,jj,kp)   = d(ii,jj,kp) + c(m)%scalar(i,j,k) * wx0*wy0*wz1
              d(ii,jp,kk)   = d(ii,jp,kk) + c(m)%scalar(i,j,k) * wx0*wy1*wz0
              d(ii,jp,kp)   = d(ii,jp,kp) + c(m)%scalar(i,j,k) * wx0*wy1*wz1
              d(ip,jj,kk)   = d(ip,jj,kk) + c(m)%scalar(i,j,k) * wx1*wy0*wz0
              d(ip,jj,kp)   = d(ip,jj,kp) + c(m)%scalar(i,j,k) * wx1*wy0*wz1
              d(ip,jp,kk)   = d(ip,jp,kk) + c(m)%scalar(i,j,k) * wx1*wy1*wz0
              d(ip,jp,kp)   = d(ip,jp,kp) + c(m)%scalar(i,j,k) * wx1*wy1*wz1
			  endif
           enddo
        enddo
     enddo
  enddo

  CALL MPI_ALLREDUCE(d,dall,size(d),mpi_real,MPI_SUM,myComm1,IERROR)
  DEALLOCATE(d)

	  L=0;  
      do I = 1, cc%nxd
	     do J = 1, cc%nyd
	        do K = 1, cc%nzd
			     L      = L+1;
                 dall(i,j,k)   = dall(i,j,k) /max(cc%diagweight(1,L),1.e-6_p2)
           enddo
        enddo
     enddo
  


 
  return
end subroutine Adconvert



!=================================================================================================
! convert vector quantity from cylindricall coordinate to cartesian
subroutine cconvert(c,d)
  use global_parameters
  use diagnos
  implicit none
  include 'mpif.h'
  real(p2),pointer,dimension(:,:,:) :: c,d 
  integer :: m,k
!-----------------  
  d = c
  do m=mgrid0,mgrid1
     do k = nz0,nz1
        d(1,k,m) = c(1,k,m) * diagCartesian(1,k,m) + &
		           c(2,k,m) * diagCartesian(2,k,m) + &
		           c(3,k,m) * diagCartesian(3,k,m)
        d(2,k,m) = c(1,k,m) * diagCartesian(4,k,m) + &
		           c(2,k,m) * diagCartesian(5,k,m) + &
		           c(3,k,m) * diagCartesian(6,k,m)
        d(3,k,m) = c(1,k,m) * diagCartesian(7,k,m) + &
		           c(2,k,m) * diagCartesian(8,k,m) + &
		           c(3,k,m) * diagCartesian(9,k,m)
     enddo
  enddo
  return
end subroutine cconvert

!=================================================================================================

subroutine checkdiag()
  use global_parameters
  use grid_data
  use diagnos
  implicit none
  include 'mpif.h'
  integer :: flag,ierror,ndims(3)
  real(p2)    :: ddomains(6) 
  real    :: timecpu0(0:20)

  flag = 0
  call rereadParameters(flag)
  if(flag >=20) then

!     if(allocated(diagweight))    deallocate(diagweight)
!     if(allocated(LL))            deallocate(LL)
!     if(allocated(xdiag))         deallocate(xdiag)
!     if(allocated(ydiag))         deallocate(ydiag)
!     if(allocated(zdiag))         deallocate(zdiag)
!	 if(allocated(diagcartesian)) deallocate(diagcartesian)

!     allocate(diagweight(6,nxd*nyd*nzd),LL(nxd,nyd,nzd))
!	 allocate(xdiag(nxd),ydiag(nyd),zdiag(nzd))
!    allocate(diagcartesian(9,nz0:nz1,mgrid))

!	 call stconvert
     ndims    = (/nxd,nyd,nzd/)
     ddomains = ddomain  !(/dxmin,dxmax,dymin,dymax,dzmin,dzmax/)
     call setupDiag(box,ndims,ddomains)

  endif


  if(mype==0 .and. stdout /= 6) then
     close(unit=stdout)
     open(stdout,file='stdout.out',status='unknown',position='append')
  endif
  call mpi_allreduce(timecpu,timecpu0,21,mpi_real,mpi_sum,mpi_comm_world, ierror)


  if(mype==0) then
     write(stdout,'(5a14)') 'initial','ions','electrons','field','other'  !,'operator'
     write(stdout,'(5f14.1)') timecpu0(1)/60.,timecpu0(2)/60.,timecpu0(3)/60.,timecpu0(4)/60.,timecpu0(5)/60.
  endif

  return
end subroutine checkdiag

!=================================================================================================
subroutine diagnos_particle
  use global_parameters
  use grid_data
  use vector_functions
  use diagnos
  implicit none
  include 'mpif.h'
  integer :: nions
  real,pointer,dimension(:,:) :: send,recv,dummy,save
  integer,pointer,dimension(:) :: num_recv,request

  real(p2),pointer,dimension(:,:) :: temp 
  integer :: i,num,m,L,ierror,Lt
  integer :: nregion,n,kind,nbuff,ier,iplot,loop,loops
  real(p2) :: r,theta,zeta,vx,vy,vz,x,y,z,vper,amass,p,q,w,param(20)
  real(p2) :: vb(3),vv(3),vpar,&
               W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,angle
  integer  :: ii,jj,kk,ip,jp,kp

  integer :: ncid,dimid(5),varid(30),dimids(2),istart(2),icount(2)
  integer :: status(mpi_status_size)
  real(8) :: xyz(3)
  character*5 fldname(40)
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename


  CHARACTER*5 CT
  interface
    subroutine dconvert(a,b)
       use constants
       real(p2),pointer,dimension(:,:,:) :: a,b
    end subroutine dconvert
    subroutine cconvert(a,b)
       use constants
       real(p2),pointer,dimension(:,:,:) :: a,b
    end subroutine cconvert

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
!------------------------
  if(n_pregion==0) return

  idiagP = idiagP + 1
  write(ct,'(I5.5)')IDIAGP 
  if( ghybrid .and. mod(ndiagP*dt,1._p2) < 0.1  ) write(ct,'(I5.5)')int(stime + 2.*dt)

  
  fldname = (/"reg01","reg02","reg03","reg04","reg05","reg06","reg07","reg08","reg09","reg10", &
              "reg11","reg12","reg13","reg14","reg15","reg16","reg17","reg18","reg19","reg20", &
              "reg21","reg22","reg23","reg24","reg25","reg26","reg27","reg28","reg29","reg30", &
              "reg31","reg32","reg33","reg34","reg35","reg36","reg37","reg38","reg39","reg40"/)


  param(1) = aion
  param(2) = aelectron
  param(3) = ulength
  param(4) = utime
  param(5) = uspeed
  param(6) = utemperature
  param(7) = ubfield
  param(8) = uefield
  param(9) = cspeed

! ---------------------------------------------------------------
  loop     = 2     ! including ion & electron
  if(ghybrid .or. hybrid) loop = 1
  filename = 'particle'//ct//'.dat'

  Loop = 1
  do i = 1,int(pdomain(3))
     if(pdomain(3+i)>10 ) loop = 2  
  enddo

  if( mype == 0) then
      OPEN(33, FILE = 'particle'//ct//'.dat')  
      Write(33,*)Stime, n_pregion,6,9,kinds
      Write(33,1106)domain_pregion(1:n_pregion,:)
	  write(33,*) 20
	  write(33,*) param
	  write(33,*)loop
  endif

  do Loops = 1,Loop

     if( mype == 0) then
	    if(loops==1) write(33,*)fraci(1:kinds),qions(1:kinds),mions(1:kinds)
	    if(loops==2) write(33,*)frace(1:kinds),qions(1:kinds)*0-1.,mions(1:kinds)*0+1.
     endif

  do n=1,n_pregion
     nbuff = max(maxval(block%mi),2000)*2
     allocate(send(9, nbuff),num_recv(numberpe),request(numberpe))
     num = 0

     if(Loops ==1) then
     do m=1,mblocks
           do L = 1,block(m)%mi
              kind  = block(m)%ion(L)%kind
	          amass = mions(kind)*aion
              p     = block(m)%ion(L)%p(1)
              q     = block(m)%ion(L)%p(2)
              w     = block(m)%ion(L)%p(3)
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

			  if(x>= domain_pregion(n,1) .and. x< domain_pregion(n,2) .and. & 
			     y>= domain_pregion(n,3) .and. y< domain_pregion(n,4) .and. & 
			     z>= domain_pregion(n,5) .and. z< domain_pregion(n,6) )then
				 num = num +1
				 if(num >= nbuff) then
                    allocate(dummy(9,num-1))
					dummy = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(9,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
				 endif

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

				 send(1,num) = block(m)%ion(L)%kind
				 send(2,num) = x
				 send(3,num) = y
				 send(4,num) = z
				 send(5,num) = block(m)%ion(L)%v(1)/amass
				 send(6,num) = block(m)%ion(L)%v(2)/amass
				 send(7,num) = block(m)%ion(L)%v(3)/amass
				 send(8,num) = vpar
				 send(9,num) = block(m)%ion(L)%w *fraci(kind)
			  endif	  			   
		   enddo
     enddo

     else if(Loops==2) then
     do m=1,mblocks
           do L = 1,block(m)%me
              kind  = block(m)%ele(L)%kind
              p     = block(m)%ele(L)%p(1)
              q     = block(m)%ele(L)%p(2)
              w     = block(m)%ele(L)%p(3)
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

			  if(x>= domain_pregion(n,1) .and. x< domain_pregion(n,2) .and. & 
			     y>= domain_pregion(n,3) .and. y< domain_pregion(n,4) .and. & 
			     z>= domain_pregion(n,5) .and. z< domain_pregion(n,6) )then
				 num = num +1
				 if(num >= nbuff) then
                    allocate(dummy(9,num-1))
					dummy = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(9,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
				 endif

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
				 vv      = block(m)%ele(L)%v/aelectron
				 vpar    = vb(1)*vv(1) + vb(2)*vv(2) + vb(3)*vv(3)

				 send(1,num) = block(m)%ele(L)%kind
				 send(2,num) = x
				 send(3,num) = y
				 send(4,num) = z
				 send(5,num) = block(m)%ele(L)%v(1)/aelectron
				 send(6,num) = block(m)%ele(L)%v(2)/aelectron
				 send(7,num) = block(m)%ele(L)%v(3)/aelectron
				 send(8,num) = vpar
				 send(9,num) = block(m)%ele(L)%w *frace(kind)
			  endif	  			   
		   enddo
     enddo
	 endif


     call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)
     if(mype > 0 .and. num > 0) then
	       call mpi_isend(send(:,1:num),9*num, mpi_real,0, &
		                           mype, MPI_COMM_WORLD, request(mype),ier)

     endif
     if(mype==0) then
	    allocate(save(9,nbuff)) !nbuff
	    iplot = 0
        do i=1,num
		   iplot = iplot +1
		   istart = (/1,iplot/)
		   if(iplot >= nbuff) then
              allocate(dummy(9,iplot-1))
			  dummy = save(:,1:iplot-1)
			  deallocate(save)
			  nbuff = nbuff*1.5
			  allocate(save(9,nbuff))
			  save(:,1:iplot-1) = dummy
			  deallocate(dummy)
		   endif
           save(:,iplot) = send(:,i)
        enddo
     endif

     do L=1,numberpe-1
	    if(num_recv(L+1) > 0) then
	       if(mype==0) then
		      allocate(recv(9,num_recv(L+1)))
              call mpi_recv(recv,9*num_recv(L+1), mpi_real,L, &
		                           L, MPI_COMM_WORLD, status, ier)
		  
		      do i=1,num_recv(L+1)
		         iplot = iplot +1
		         istart = (/1,iplot/)
		         if(iplot >= nbuff) then
                    allocate(dummy(9,iplot-1))
			        dummy = save(:,1:iplot-1)
			        deallocate(save)
			        nbuff = nbuff*1.5
			        allocate(save(9,nbuff))
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
	    write(33,*)9,iplot
		write(33,1106)save(:,1:iplot)
		deallocate(save)
	 endif
     deallocate(send,request,num_recv)
  enddo
  enddo
  if(mype==0) close(33)

1106  format(5e12.4)
   
  return
end subroutine diagnos_particle

! =================================================================================================
! output fe, dfe, J=int(vf dv) dJ=int(v df dv) at each time
! v in unit of vthe  and int(f dv) = 1

subroutine DiagEdistributionLHWcurentDrive
  use global_parameters
  use vector_functions
  implicit none
  include 'mpif.h'
  integer :: nions
  real(p2),pointer,dimension(:) :: vpar, fe,dfe,vfe,vdfe, fet,dfet,vfet,vdfet
  real(p2),pointer,dimension(:,:) :: QVP_DIAG
  real(p2),pointer,dimension(:,:) :: QVP_DIAG_ALL
  integer :: nv,i,num,numt,ic,m,ierror
  real(p2) :: vminp,vmaxp,dv,je,dje,jet,djet,qq,w1,w2,vfe_,vdfe_,vx,vy,vz,we,v
  save ic
  data ic/0/
!------------------------
  ic     = ic +1
! design a vpar [-10,10] and dv = 0.1
  vmaxp   =  10.
  vminp   = -10.
  nv      = 201
  dv      = (vmaxp-vminp)/real(nv-1)
  allocate(vpar(nv),fe(nv),dfe(nv),vfe(nv),vdfe(nv),fet(nv),dfet(nv),vfet(nv),vdfet(nv))
  do i = 1, nv
     vpar(i) = vminp + (i-1)*dv
  enddo
  fe    = 0.
  dfe   = 0.
  vfe   = 0.
  vdfe  = 0.
  je    = 0.
  dje   = 0.
  num   = 0

  DO m=1,me
        num   = num +1
		if(vx < vminp)  vx = vminp
		if(vx >=vmaxp)  vx = vmaxp-0.01
		qq    = (vx -vminp)/dv+1.

		i     = int(qq)
		w2    = qq - i
		w1    = 1.- w2
  ENDDO

  CALL MPI_ALLREDUCE(num,numt,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
  if(p2==8) then
  call mpi_allreduce(fe,fet,nv,mpi_double_precision,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(dfe,dfet,nv,mpi_double_precision,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(vfe,vfet,nv,mpi_double_precision,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(vdfe,vdfet,nv,mpi_double_precision,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(je,jet,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(dje,djet,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD, ierror)
  else
  call mpi_allreduce(fe,fet,nv,mpi_real,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(dfe,dfet,nv,mpi_real,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(vfe,vfet,nv,mpi_real,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(vdfe,vdfet,nv,mpi_real,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(je,jet,1,mpi_real,mpi_sum,MPI_COMM_WORLD, ierror)
  call mpi_allreduce(dje,djet,1,mpi_real,mpi_sum,MPI_COMM_WORLD, ierror)
  endif

  fet   = fet/real(numt)
  dfet  = dfet/real(numt)
  vfet  = vfet/real(numt)
  vdfet = vdfet/real(numt)

! int_(vfe dv)
  vfe_    = 0.
  vdfe_   = 0.
  do i = 1,nv-1
     v     = (vpar(i)+vpar(i+1))/2.
	 vfe_  = vfe_ + v*(fe(i)+fe(i+1))/2. * dv
	 vdfe_ = vdfe_ + v*(dfe(i)+dfe(i+1))/2. * dv
  enddo


  if(ic == 1) then
     IF(mype==0) then
        OPEN(11, FILE = 'fe.dat')      !, FORM="unformatted")
        write(11,*)nv
        write(11,1106)vpar
     ENDIF
  else
      if(mype==0) OPEN(11, FILE = 'fe.dat',position='append')      !, FORM="unformatted")
  endif

  IF(mype==0) then
     write(11,*)stime
     write(11,1106)fet,dfet
     write(11,1106)jet/real(numt),djet/real(numt),vfe_,vdfe_
     close(11)
  ENDIF

1106  format(5e12.4)
   
  DEALLOCATE(vpar, fe,dfe,vfe,vdfe, fet,dfet,vfet,vdfet)
  return
end subroutine DiagEdistributionLHWcurentDrive

!=================================================================================================
subroutine particle_trace
  use global_parameters
  use grid_data
  implicit none
  include 'mpif.h'
  integer :: nions,ic,nskip,l,m,mn
  real(p2),pointer,dimension(:,:) :: temp 
  real(p2) :: vx,vy,vz,vb(3),vb0(3),vpar,vper,mu,mu0
  real*8   ::  p,q,w,     xyz(3)

  CHARACTER*4 CT
  save ic
  data ic/0/
  interface
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

      function dipole2DB(a,b,c)
        use constants
	    real(p2):: a,b,c,dipole2DB(3)
	  end function dipole2DB

  end interface
!------------------------
  if(kpe > 1  ) return
  if(mype /= 0) return
  ic=ic+1
  write(ct,'(I4.4)')ic
! do the first time of coefficient
  nions = 50
  nskip = max(1,block(1)%mi/nions)
  nions = block(1)%mi/nskip

  if(nions <= 1) then
     nions = 50
     nskip = max(1,sum(block(:)%me/nions))
     nions = sum(block(:)%me)/nskip
  endif

  if(qsload == 1) nskip =1

!  nions = 10
!  nskip = 1

  if(ic == 1) then
     IF(mype==0) then
        OPEN(11, FILE = 'FeFi_ions.dat')      
        OPEN(12, FILE = 'FeFi_eles.dat')      
		if(coordinate == 0) then
		  write(11,*)8,nions,pmin,pmax,wmin,wmax
		  write(12,*)8,nions,pmin,pmax,wmin,wmax
		else
		  write(11,*)8,nions,-a1,a1,-a1,a1
		  write(12,*)8,nions,-a1,a1,-a1,a1
		endif
     ENDIF

     do m=1,mblocks
        L = 0
        do mn=1,block(m)%me
           if(mod(mn,nskip)==0 ) then
	          L=L+1
              block(m)%ele(L)%p  = block(m)%ele(mn)%p
              block(m)%ele(L)%p0 = block(m)%ele(mn)%p0
              block(m)%ele(L)%v  = block(m)%ele(mn)%v
              block(m)%ele(L)%v0 = block(m)%ele(mn)%v0
              block(m)%ele(L)%mark(1) = L + nions*(m-1)
           endif
        enddo
		block(m)%me = L
     enddo




  else
      if(mype==0) OPEN(11, FILE = 'FeFi_ions.dat',position='append')      !, FORM="unformatted")
      if(mype==0) OPEN(12, FILE = 'FeFi_eles.dat',position='append')      !, FORM="unformatted")
  endif
  allocate(temp(8,nions))
  temp  = 0.
  L     = 0
  do m=1,mblocks
     do mn=1,block(m)%mi
        if(mod(mn,nskip)==0 .and. L < nions) then
	    L=L+1
		p     = block(m)%ion(mn)%p(1)
        q     = block(m)%ion(mn)%p(2)
	    w     = block(m)%ion(mn)%p(3)
		vx    = block(m)%ion(mn)%v(1)/aion
		vy    = block(m)%ion(mn)%v(2)/aion
		vz    = block(m)%ion(mn)%v(3)/aion
		xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)


		temp(1,L) = xyz(1)
		temp(2,L) = xyz(2)
		temp(3,L) = xyz(3)
		temp(4,L) = vx
		temp(5,L) = vy
		temp(6,L) = vz

      endif
    enddo
  enddo
      write(11,*) sTime
      write(11,1102)temp
	  close(11)

! do electrons -------------------------
  temp  = 0.
  L     = 0
  do m=1,mblocks
     do mn =1,block(m)%me
!       if(mod(mn,nskip)==0 .and. L < nions) then
	    L=L+1
		p     = block(m)%ele(mn)%p(1)
        q     = block(m)%ele(mn)%p(2)
	    w     = block(m)%ele(mn)%p(3)
		vx    = block(m)%ele(mn)%v(1)/aelectron
		vy    = block(m)%ele(mn)%v(2)/aelectron
		vz    = block(m)%ele(mn)%v(3)/aelectron
		xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)

		temp(1,L) = xyz(1)
		temp(2,L) = xyz(2)
		temp(3,L) = xyz(3)
		temp(4,L) = vx
		temp(5,L) = vy
		temp(6,L) = vz
        temp(7,L) = block(m)%ele(mn)%mark(1)
        temp(8,L) = p
!       endif
     enddo
  enddo
      write(12,*) sTime
      write(12,1102)temp
	  close(12)

  deallocate(temp)
   

1102  format(5e12.4)
  return
end subroutine particle_trace
! ----------------------------------------------------------------------------------
subroutine diagnos_field1
  use global_parameters
  use grid_data
  use diagnos
  implicit none
  include 'mpif.h'
  real,    pointer,dimension(:,:,:) :: plot
  real(p2),pointer,dimension(:,:,:) :: test1
  real(p2),pointer,dimension(:,:,:,:) :: vec1,vec2,test3
  real(p2) :: param(80)
  integer :: ic,i,iplot,ifield,m,ifchecki,ifchecke,n,ifcheck(30)

  type(bscalar_type), dimension(:), pointer    :: temp,beta_i,beta_e
  type(bvector_type), dimension(:), pointer    :: bfld,efld,ui,ue,delJ,dumm,av1,delB,delE

  CHARACTER*5 CT
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
    subroutine gather_vector(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a
       real(p2),pointer,dimension(:,:,:,:) :: b
    end subroutine gather_vector
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

    subroutine bbcellecenter(a,b)
       use constants
	   use grid_data
	   type(bvector_type), dimension(:), pointer    :: a,b
    end subroutine bbcellecenter

  end interface
!------------------------
  ic = ic +1
  idiagf1 = idiagf1 + 1
  write(ct,'(I5.5)')IDIAGf1 
! do the first time of coefficient
!  if(ic == 1) then
!     allocate(diagweight(6,nxd*nyd*nzd),LL(nxd,nyd,nzd))
!	 allocate(xdiag(nxd),ydiag(nyd),zdiag(nzd))
!	 call stconvert
!  endif
  param     = 0.
  param(1)  = stime  ! has been scaled to omega_i
  param(2)  = aelectron
  param(3)  = aion
  param(4)  = br0
  param(5)  = bt0
  param(6)  = bz0
  param(7)  = betae
  param(8)  = dxd
  param(9)  = dyd
  param(10) = dzd
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
  param(23) = ulength*rhos
  param(28) = ua0
  param(29) = ua0
  param(30) = uj0
  param(31) = ua0

! some key parameters
  param(41) = omega_lh/omega_i
  param(75-nfield+1:75)= DiagQuantity(1:nfield) 
  param(77) = simtype(2)
  param(78) = simtype(1)
  param(79) = coordinate
  param(80) = nfield   ! number of field quantity are recorded


  IF(mype==0) THEN
     if(idiagbox1 ==1)then
	    open(11, FILE = 'wave'//ct//'.dat')      !, FORM="unformatted")
           write(11,*)3,80
           write(11,*)box1%NXD,box1%NYD,box1%NZD
           write(11,*)param
           write(11,1102)box1%xdiag,box1%ydiag,box1%zdiag
     else
	    open(11, FILE = 'wave'//ct//'.dat')      !, FORM="unformatted")
           write(11,*)3,80
           write(11,*)cuv1%NXD,cuv1%NYD,cuv1%NZD
           write(11,*)param
           write(11,1102)cuv1%xdiag,cuv1%ydiag,cuv1%zdiag

	 endif
  ENDIF
!--------------------------------------------------------
  if(idiagbox1 ==1) allocate(plot(box1%nxd,box1%nyd,box1%nzd))
  if(idiagbox1 ==0) allocate(plot(cuv1%nxd,cuv1%nyd,cuv1%nzd))

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
  if(idiagbox ==1) then
     do m=1,mblocks
        dumm(m)%vector = bv1(m)%vector  + bv0(m)%vector 
     enddo
	 call bconcarB(dumm,bfld)
  else
     do m=1,mblocks
	    do n=1,3
           dumm(m)%vector(n,:,:,:) = (bv1(m)%vector(n,:,:,:)  + bv0(m)%vector(n,:,:,:)) * &
                                      block(m)%node%h(2,n,:,:,:)
        enddo
     enddo
	 call bbcellecenter(dumm,bfld)
  endif
! (2) e-field -------
  if(idiagbox ==1) then
     do m=1,mblocks
        dumm(m)%vector = ev1(m)%vector 
     enddo
	 call bcovcarE(dumm,efld)
  else
     do m=1,mblocks
	    do n=1,3
           dumm(m)%vector(n,:,:,:) =  ev1(m)%vector(n,:,:,:)   / &
                                      block(m)%node%h(1,n,:,:,:)
        enddo
     enddo
     call becellecenter(dumm,efld)
  endif
! (3) ui-field -------
  if(idiagbox ==1) then
     do m=1,mblocks
        dumm(m)%vector = ji(m)%vector 
     enddo
	 call bcovcarE(dumm,ui)
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
        if(.not. df)ui(m)%vector(n,:,:,:) = ui(m)%vector(n,:,:,:)/max(ni(m)%scalar,1.e-4_p2)
        if(df)ui(m)%vector(n,:,:,:) = ui(m)%vector(n,:,:,:)/max(ni0(m)%scalar,1.e-4_p2)
	 enddo
  enddo
! (4) ue-field --------------
  if(simtype(1) == 6 .or. simtype(1) == 7 .or. simtype(1) == 66 .or. simtype(1) == 67) then
	 ue(m)%vector(1,:,:,:) = block(m)%ionf(1,:,:,:)%v(1)
	 ue(m)%vector(2,:,:,:) = block(m)%ionf(1,:,:,:)%v(2)
	 ue(m)%vector(3,:,:,:) = block(m)%ionf(1,:,:,:)%v(3)
  else
     if(idiagbox ==1) then
        do m=1,mblocks
           dumm(m)%vector = je(m)%vector 
        enddo
	    call bcovcarE(dumm,ue)
     else
        do m=1,mblocks
	       do n=1,3
              dumm(m)%vector(n,:,:,:) = je(m)%vector(n,:,:,:)*block(m)%node%h(3,n,:,:,:)
           enddo
        enddo
	    call becellecenter(dumm,ue)
     endif

     do m=1,mblocks
	    do n=1,3
           if(.not. df)ue(m)%vector(n,:,:,:) = ue(m)%vector(n,:,:,:)/max(ne(m)%scalar,1.e-4_p2)
           if(df)ue(m)%vector(n,:,:,:) = ue(m)%vector(n,:,:,:)/max(ne0(m)%scalar,1.e-4_p2)
	    enddo
     enddo
  endif



  ifchecki = 0
  ifchecke = 0
  ifcheck  = 0
  do iplot = 1,nfield1
     ifield = DiagQuantity1(iplot)
	 if(ifield >= 40 .and. ifield <= 58 ) then
	    ifchecki =   ifchecki + 1
	    if(ifchecki ==1) call getUT(0,ions)
     endif
	 if(ifield >= 59 .and. ifield<=70 ) then
	    ifchecke =   ifchecke + 1
	    if(ifchecke ==1  ) call getUT(0,eles)
     endif


	 if(ifield  == 21 .or. ifield  == 22 .or. ifield  == 23 ) then  ! delta_j
	    ifcheck(3) = ifcheck(3) +1
		if(ifcheck(3) ==1) then
  	       do m=1,mblocks
		      dumm(m)%vector    = bv1(m)%vector  + bv0(m)%vector 
		   enddo
		   call get_j(dumm,delJ)
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
	 if(ifield  == 24 .or. ifield  == 25 .or. ifield  == 26 ) then  ! delta_j
	    ifcheck(5) = ifcheck(5) +1
		if(ifcheck(5) ==1) then
           if(dipole .or. coordinate  ==8 ) then
		      do m=1,mblocks
	             do n =1,3
                 delB(m)%vector(n,:,:,:) = bv1(m)%vector(n,:,:,:) * block(m)%node%h(1,n,:,:,:)
		         enddo
			  enddo
		   else
  	          call field3(delB)
           endif
        endif
	 endif

	 if(ifield  == 27 .or. ifield  == 28 .or. ifield  == 29 ) then  ! delta_E
	    ifcheck(6) = ifcheck(6) +1
		if(ifcheck(6) ==1) then
           if( dipole .or. coordinate  ==8 ) then
		      do m=1,mblocks
	             do n =1,3
                 delE(m)%vector(n,:,:,:) = ev1(m)%vector(n,:,:,:) / block(m)%node%h(2,n,:,:,:)
		         enddo
			  enddo
		   else
           endif
        endif
	 endif
	 if(ifield  == 71) then  ! beta_i
        do m=1,mblocks
		   beta_i(m)%scalar = 0.
		   do n=1,kinds
		      beta_i(m)%scalar = beta_i(m)%scalar + block(m)%ini(n,:,:,:)%ni0 * &
			                     (block(m)%ini(n,:,:,:)%Ti0(1)+block(m)%ini(n,:,:,:)%Ti0(2)*2)/3.
		   enddo
		   temp(m)%scalar = bfld(m)%vector(1,:,:,:)**2+bfld(m)%vector(2,:,:,:)**2+bfld(m)%vector(3,:,:,:)**2
        enddo
        do m=1,mblocks
		   beta_i(m)%scalar = beta_i(m)%scalar/max(temp(m)%scalar,1.e-4_p2)*alpha
        enddo
	 endif  

  enddo

  do iplot = 1, nfield1
     ifield = DiagQuantity1(iplot)
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
		   enddo
        case(18,19,20)                        !ue
		   do m=1,mblocks
		      temp(m)%scalar = ue(m)%vector(ifield-17,:,:,:)
		   enddo
        case(21,22,23)                        !delta J
		   do m=1,mblocks
		      temp(m)%scalar = delJ(m)%vector(ifield-20,:,:,:)
		   enddo
        case(24,25,26)                        !dBvec,b||and Bper
		   do m=1,mblocks
		      temp(m)%scalar = delB(m)%vector(ifield-23,:,:,:)
		   enddo
        case(27,28,29)                        !delE
		   do m=1,mblocks
		      temp(m)%scalar = delE(m)%vector(ifield-26,:,:,:)
		   enddo
!        case(27,28,29)                        !dj X B0
!          allocate(delJ(3,nz0:nz1,mgrid0:mgrid1))
!		   call get_j(bv1,delJ)
!		   call cross(delJ,Bv_bar,vec1)
!           temp = vec1(ifield-26,:,:)
!	       call dconvert(temp,plot)
!		   deallocate(delJ)	 
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
           if(mype==0) write(11,1102)PLOT	 
        case(39)                        !divPpe/ne
!           temp = dUedt(3,:,:)
! ----------------------------------------------plasma informations ------------------
        case(41,42,43)        !ui(1-3,:)
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%vi0(ifield-40)
		   enddo
        case(44)              !Ni
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ni0
		   enddo
        case(45,46)           !ui(1-3,:),Ti
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ti0(ifield-44)
		   enddo
        case(47,48,49)        !ui(1-3,:) kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%vi0(ifield-46)
		   enddo
        case(50)              !ui(1-3,:),Ni, kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ni0
		   enddo
        case(51,52)           !ui(1-3,:),Ti, kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ti0(ifield-50)
		   enddo

        case(53,54,55)        !ui(1-3,:) kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%vi0(ifield-52)
		   enddo
        case(56)              !ui(1-3,:),Ni, kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%ni0
		   enddo
        case(57,58)           !ui(1-3,:),Ti, kind=3
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(3,:,:,:)%ti0(ifield-56)
		   enddo
! -------------------------------------------------------------------------
        case(59,60,61)        !ue(1-3,:), kind=1
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ve0(ifield-58)
		   enddo
        case(62)              !Ne
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%ne0
		   enddo
        case(63,64)           !Te
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(1,:,:,:)%te0(ifield-62)
		   enddo
        case(65,66,67)        !ue(1-3,:),kind=2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ve0(ifield-64)
		   enddo
        case(68)              !Ne --- 2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%ne0
		   enddo
        case(69,70)           !Te ----2
		   do m=1,mblocks
		      temp(m)%scalar = block(m)%ini(2,:,:,:)%te0(ifield-68)
		   enddo
        case(71)              !beta_i
		   do m=1,mblocks
		      temp(m)%scalar = beta_i(m)%scalar
		   enddo
        case(72)              !beta_e
		   do m=1,mblocks
		      temp(m)%scalar = beta_e(m)%scalar
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
	       
     if(idiagbox1 ==1) call dconvert(temp,plot,box1)
	 if(idiagbox1 ==0) call dconvert(temp,plot,cuv1)
     if(mype==0) write(11,1102)PLOT	 

  enddo
  if(mype==0) close(11)
! 1102  format(5f14.5)
  1102  format(5e16.5)
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
  return
end subroutine diagnos_field1

! --------------------------------------------------------------------

! f(vapr,vper),f(E,alpha)
subroutine Diag_Pdistributions
  use global_parameters
  use grid_data
  implicit none
  include 'mpif.h'
  integer :: nions
  real(p2),pointer,dimension(:) :: xvpar, xvper,xee,xalpha,nnn
  real(p2),pointer,dimension(:,:,:) :: fvxy,fEalpha,afvxy,afEalpha

  integer :: nv,i,j,ic,m,n,ierror,nregion,kk,ij,mkind,mn,L
  real(p2) :: vminp,vmaxp,dv,dv2,dE,dje,jet,djet,qq,w1,w2,w3,w4, &
              w11,w12,w21,w22,vfe_,vdfe_,vx,we,v,vv(3),avi_sbbar(3),qx,qy,avi_btbar,ee,&
			  dvpar,dvper,dalpha,emin,emax,weight,theta
  real(p2) :: x,y,z,wz1,wz0,wp0,wt00,wt10,wp1,wt01,wt11,charge,amass,vpar,vtot,vper,num(3),anum(3)

  
  CHARACTER*5 CT
  save ic
  data ic/0/
!------------------------

  ic     = ic +1
  write(ct,'(I5.5)')IDIAGf 

! f(vx,vy),f(vx,vz),f(vy,vz), f(v||),f(vper),f(E) in three region 

  vmaxp   =  10.
  vminp   = -10.
  nv      =  201

  dvpar   =  (vmaxp-vminp)/real(nv-1)
  dvper   =  (vmaxp-0.)  /real(nv-1)

  Emin    = -4.  !alog (E)
  Emax    =  2.
  dE      =  (Emax-Emin)/real(nv-1)
  dalpha  =  180./real(nv-1)

  allocate(xvpar(nv),xvper(nv),xee(nv),xalpha(nv),nnn(kinds))
  allocate(fvxy(kinds,nv,nv),fEalpha(kinds,nv,nv))
  allocate(afvxy(kinds,nv,nv),afEalpha(kinds,nv,nv))

  do i = 1, nv
     xvpar(i)  =  vminp  + (i-1)*dvpar
     xvper(i)  = (i-1.) * dvper
     xee(i)    =  Emin  + (i-1.) * dE
	 xalpha(i) = (i-1.) * dalpha
  enddo


  fvxy      = 0.
  fEalpha   = 0.

  num       = 0.
  nnn       = 0.

  DO m=1,mblocks 
     do L =1, block(m)%mi
           mkind    = block(m)%ion(L)%kind
           we       = block(m)%ion(L)%w * fraci(mkind)
	       charge   = qion*qions(mkind)
	       amass    = aion*mions(mkind)

	       vv(1)    = block(m)%ion(L)%v(1)/amass/va
	       vv(2)    = block(m)%ion(L)%v(2)/amass/va
	       vv(3)    = block(m)%ion(L)%v(3)/amass/va

	       vpar     = (vv(1)*br0+vv(2)*bt0+vv(3)*bz0)/sqrt(br0**2+bt0**2+bz0**2)

		   vtot     = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
		   theta    = acos(vpar/vtot) *180./pi
		   vper     = sqrt(max(vtot**2-vpar**2,0._p2))
		   ee       = log10(0.5*amass*vtot**2 / amass)
!		   if(ee<emin) ee = emin


! --- f(vpar,vper) ----------------------------------
           if(abs(vpar)<= vmaxp .and. vper <=vmaxp) then   
		      nnn(mkind)= nnn(mkind) + we
		      qx    = (vpar -vminp) /dvpar + 1.
		      qy    = (vper -0.  ) /dvper + 1.
		      i     = min(int(qx),nv-1)
		      j     = min(int(qy),nv-1)
		      w2    = qx - i
		      w1    = 1.- w2
			  w4    = qy - j    
			  w3    = 1. -w4    
              w11   = w1*w3
              w12   = w1*w4
              w21   = w2*w3
              w22   = w2*w4

			  fvxy(mkind,i,j)      = fvxy(mkind,i,j)       + w11 * we /vper
			  fvxy(mkind,i,j+1)    = fvxy(mkind,i,j+1)     + w12 * we /vper
			  fvxy(mkind,i+1,j)    = fvxy(mkind,i+1,j)     + w21 * we /vper
			  fvxy(mkind,i+1,j+1)  = fvxy(mkind,i+1,j+1)   + w22 * we /vper
           endif
! -- f(vy,vz) -----------------------------------------   
           if(ee >= emin .and. ee < emax) then    
		      qx    = (theta -0.)/dalpha  + 1.
		      qy    = (ee -Emin)/de + 1.
		      i     = min(int(qx),nv-1)
		      j     = min(int(qy),nv-1)
		      w2    = qx - i
		      w1    = 1.- w2
			  w4    = qy - j    
			  w3    = 1. -w4    
              w11   = w1*w3
              w12   = w1*w4
              w21   = w2*w3
              w22   = w2*w4

			  fEalpha(mkind,i,j)      = fEalpha(mkind,i,j)       + w11 * we
			  fEalpha(mkind,i,j+1)    = fEalpha(mkind,i,j+1)     + w12 * we
			  fEalpha(mkind,i+1,j)    = fEalpha(mkind,i+1,j)     + w21 * we
			  fEalpha(mkind,i+1,j+1)  = fEalpha(mkind,i+1,j+1)   + w22 * we
            endif
      ENDDO
  enddo


  if(p2==8) then
    call mpi_allreduce(fvxy,afvxy,nv*nv*kinds,mpi_double_precision,mpi_sum,mpi_comm_world, ierror)
    call mpi_allreduce(fEalpha,afEalpha,nv*nv*kinds,mpi_double_precision,mpi_sum,mpi_comm_world, ierror)
  else
    call mpi_allreduce(fvxy,afvxy,nv*nv*kinds,mpi_real,mpi_sum,mpi_comm_world, ierror)
    call mpi_allreduce(fEalpha,afEalpha,nv*nv*kinds,mpi_real,mpi_sum,mpi_comm_world, ierror)
  endif


  if(mype==0) then
	 open(11, FILE = 'fion'//ct//'.dat')      !, FORM="unformatted")
        write(11,*)stime
        write(11,*)kinds,nv,nv
        write(11,1106)xvpar,xvper,xalpha,xee
        write(11,1106)afvxy
        write(11,1106)afEalpha
     close(11)
  endif

1106  format(5e12.4)
  deallocate(fvxy,fEalpha,afvxy,afEalpha)
  deallocate(xvpar,xvper,xee,xalpha,nnn)

  return
end subroutine Diag_Pdistributions
! ------------------------------------------------------------------------------------
! f(vapr,vper),f(E,alpha)
subroutine Diag_edistributions
  use global_parameters
  use grid_data
  implicit none
  include 'mpif.h'
  integer :: nions
  real(p2),pointer,dimension(:) :: xvpar, xvper,xee,xalpha,nnn
  real(p2),pointer,dimension(:,:,:) :: fvxy,fEalpha,afvxy,afEalpha

  integer :: nv,i,j,ic,m,n,ierror,nregion,kk,ij,mkind,mn,L
  real(p2) :: vminp,vmaxp,dv,dv2,dE,dje,jet,djet,qq,w1,w2,w3,w4, &
              w11,w12,w21,w22,vfe_,vdfe_,vx,we,v,vv(3),avi_sbbar(3),qx,qy,avi_btbar,ee,&
			  dvpar,dvper,dalpha,emin,emax,weight,theta
  real(p2) :: x,y,z,wz1,wz0,wp0,wt00,wt10,wp1,wt01,wt11,charge,amass,vpar,vtot,vper,num(3),anum(3)

  
  CHARACTER*5 CT
  save ic
  data ic/0/
!------------------------

  ic     = ic +1
  write(ct,'(I5.5)')IDIAGf 

! f(vx,vy),f(vx,vz),f(vy,vz), f(v||),f(vper),f(E) in three region 

  vmaxp   =  10.
  vminp   = -10.
  nv      =  201

  dvpar   =  (vmaxp-vminp)/real(nv-1)
  dvper   =  (vmaxp-0.)  /real(nv-1)

  Emin    = -4.  !alog (E)
  Emax    =  2.
  dE      =  (Emax-Emin)/real(nv-1)
  dalpha  =  180./real(nv-1)

  allocate(xvpar(nv),xvper(nv),xee(nv),xalpha(nv),nnn(kinds))
  allocate(fvxy(kinds,nv,nv),fEalpha(kinds,nv,nv))
  allocate(afvxy(kinds,nv,nv),afEalpha(kinds,nv,nv))

  do i = 1, nv
     xvpar(i)  =  vminp  + (i-1)*dvpar
     xvper(i)  = (i-1.) * dvper
     xee(i)    =  Emin  + (i-1.) * dE
	 xalpha(i) = (i-1.) * dalpha
  enddo


  fvxy      = 0.
  fEalpha   = 0.

  num       = 0.
  nnn        = 0.

  DO m=1,mblocks 
     do L =1, block(m)%me
           mkind    = block(m)%ele(L)%kind
           we       = block(m)%ele(L)%w * frace(mkind)
	       charge   = qelectron
	       amass    = aelectron*meles(mkind)

	       vv(1)    = block(m)%ele(L)%v(1)/amass !/va
	       vv(2)    = block(m)%ele(L)%v(2)/amass !/va
	       vv(3)    = block(m)%ele(L)%v(3)/amass !/va

	       vpar     = (vv(1)*br0+vv(2)*bt0+vv(3)*bz0)/max(sqrt(br0**2+bt0**2+bz0**2),1.e-4_p2)
		   vtot     = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
		   if(abs(vpar)>= vtot) vtot = vtot * 1.001
		   theta    = acos(vpar/vtot) *180./pi
		   vper     = sqrt(max(vtot**2-vpar**2,0._p2))
		   ee       = log10(0.5*amass*vtot**2 / amass)
!		   if(ee<emin) ee = emin

! --- f(vpar,vper) ----------------------------------
           if(abs(vpar)<= vmaxp .and. vper <=vmaxp .and. vpar >= vminp) then   
		      nnn(mkind)= nnn(mkind) + we
		      qx    = max((vpar -vminp) /dvpar,0._p2) + 1.
		      qy    = max((vper -0.  ) /dvper,1._p2) + 1.
		      i     = min(int(qx),nv-1)
		      j     = min(int(qy),nv-1)
		      w2    = qx - i
		      w1    = 1.- w2
			  w4    = qy - j    
			  w3    = 1. -w4    
              w11   = w1*w3
              w12   = w1*w4
              w21   = w2*w3
              w22   = w2*w4

			  fvxy(mkind,i,j)      = fvxy(mkind,i,j)       + w11 * we /vper
			  fvxy(mkind,i,j+1)    = fvxy(mkind,i,j+1)     + w12 * we /vper
			  fvxy(mkind,i+1,j)    = fvxy(mkind,i+1,j)     + w21 * we /vper
			  fvxy(mkind,i+1,j+1)  = fvxy(mkind,i+1,j+1)   + w22 * we /vper
           endif
! -- f(vy,vz) -----------------------------------------   
           if(ee >= emin .and. ee < emax) then    
		      qx    = (theta -0.)/dalpha  + 1.
		      qy    = (ee -Emin)/de + 1.
		      i     = min(int(qx),nv-1)
		      j     = min(int(qy),nv-1)
		      w2    = qx - i
		      w1    = 1.- w2
			  w4    = qy - j    
			  w3    = 1. -w4    
              w11   = w1*w3
              w12   = w1*w4
              w21   = w2*w3
              w22   = w2*w4

			  fEalpha(mkind,i,j)      = fEalpha(mkind,i,j)       + w11 * we
			  fEalpha(mkind,i,j+1)    = fEalpha(mkind,i,j+1)     + w12 * we
			  fEalpha(mkind,i+1,j)    = fEalpha(mkind,i+1,j)     + w21 * we
			  fEalpha(mkind,i+1,j+1)  = fEalpha(mkind,i+1,j+1)   + w22 * we
            endif
      ENDDO
  enddo


  if(p2==8) then
    call mpi_allreduce(fvxy,afvxy,nv*nv*kinds,mpi_double_precision,mpi_sum,mpi_comm_world, ierror)
    call mpi_allreduce(fEalpha,afEalpha,nv*nv*kinds,mpi_double_precision,mpi_sum,mpi_comm_world, ierror)
  else
    call mpi_allreduce(fvxy,afvxy,nv*nv*kinds,mpi_real,mpi_sum,mpi_comm_world, ierror)
    call mpi_allreduce(fEalpha,afEalpha,nv*nv*kinds,mpi_real,mpi_sum,mpi_comm_world, ierror)
  endif


  if(mype==0) then
     if(ic == 1 ) open(11, FILE = 'fele.dat')
     if(ic >  1 ) open(11, FILE = 'fele.dat',status='unknown',position='append')
	       !, FORM="unformatted")
        write(11,*)stime
        write(11,*)kinds,nv,nv
        write(11,1106)xvpar,xvper,xalpha,xee
        write(11,1106)afvxy
        write(11,1106)afEalpha
     close(11)
  endif

1106  format(5e12.4)
  deallocate(fvxy,fEalpha,afvxy,afEalpha)
  deallocate(xvpar,xvper,xee,xalpha,nnn)

  return
end subroutine Diag_edistributions


!=================================================================================================
subroutine diagnos_particle_gyro
  use global_parameters
  use grid_data
  use vector_functions
  use diagnos
  implicit none
  include 'mpif.h'
  integer :: nions
  real,pointer,dimension(:,:) :: send,recv,dummy,save
  integer,pointer,dimension(:) :: num_recv,request
  type(bvector_type), dimension(:), pointer    :: bv0_,bfld,bv1_

  real(p2),pointer,dimension(:,:) :: temp 
  integer :: i,num,m,L,ierror,Lt
  integer :: nregion,n,kind,nbuff,ier,iplot,loop,loops
  real(p2) :: r,theta,zeta,vx,vy,vz,x,y,z,vper,amass,p,q,w
  integer :: ncid,dimid(5),varid(30),dimids(2),istart(2),icount(2)
  integer :: status(mpi_status_size)
  real(8) :: xyz(3)
  character*5 fldname(40)
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),vb(3),vb1(3),phase,det
  integer  :: ii,jj,kk,ip,jp,kp
  real(p2) :: v1,v2,v3,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,angle
  real(p2) :: bw(3),vw(3),uu(3),gamma 



  CHARACTER*5 CT
  interface
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine dconvert(a,b)
        use constants
        real(p2),pointer,dimension(:,:,:) :: a,b
      end subroutine dconvert
      subroutine cconvert(a,b)
        use constants
        real(p2),pointer,dimension(:,:,:) :: a,b
      end subroutine cconvert
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
  end interface
!------------------------

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
! -------------------------------------------

  idiagP = idiagP + 1
  write(ct,'(I5.5)')IDIAGP 
  if( mod(ndiagP*dt,1._p2) < 0.1  ) write(ct,'(I5.5)')int(stime/omega_i + 2.*dt)


  n_pregion = int(pdomain(1)+0.001)
  allocate(domain_pregion(n_pregion,6))
  
  fldname = (/"reg01","reg02","reg03","reg04","reg05","reg06","reg07","reg08","reg09","reg10", &
              "reg11","reg12","reg13","reg14","reg15","reg16","reg17","reg18","reg19","reg20", &
              "reg21","reg22","reg23","reg24","reg25","reg26","reg27","reg28","reg29","reg30", &
              "reg31","reg32","reg33","reg34","reg35","reg36","reg37","reg38","reg39","reg40"/)


  L  = 1

  if(ghybrid) then
    do m=1,n_pregion
	   L = L +1
	   domain_pregion(m,1) = pdomain(L)
	   L = L +1
	   domain_pregion(m,2) = pdomain(L)
	   L = L +1
	   domain_pregion(m,3) = pdomain(L)
	   L = L +1
	   domain_pregion(m,4) = pdomain(L)
	   L = L +1
	   domain_pregion(m,5) = pdomain(L)
	   L = L +1
	   domain_pregion(m,6) = pdomain(L)
	enddo

  else

    L  =  2
    do m=1,n_pregion
	   domain_pregion(m,1) = real(int(pdomain(L)*nx * pmove)) - &
	                         max(real(int(pdomain(L+1))),1.)
	   domain_pregion(m,2) = real(int(pdomain(L)*nx * pmove)) + &
	                         max(real(int(pdomain(L+1))),1.)!max(real(int(pdomain(L)*nx)),1.)
	   L = L +2
	   domain_pregion(m,3) = real(int(pdomain(L)*ny * qmove))- &
	                         max(real(int(pdomain(L+1))),1.)
	   domain_pregion(m,4) = real(int(pdomain(L)*ny * qmove))+ &
	                         max(real(int(pdomain(L+1))),1.)
	   L = L +2
	   domain_pregion(m,5) = real(int(pdomain(L)*nz * wmove))- &
	                         max(real(int(pdomain(L+1))),1.)
	   domain_pregion(m,6) = max(real(int(pdomain(L)*nz)),1.)+ &
	                         max(real(int(pdomain(L+1))),1.)
       L = L +2
	enddo
  endif
! ---------------------------------------------------------------
  loop     = 2     ! including ion & electron
  if(ghybrid) loop = 1
  filename = 'particle'//ct//'.dat'
  if( mype == 0) then
      OPEN(33, FILE = 'particle'//ct//'.dat')  
      Write(33,*)Stime, n_pregion,6,9,kinds
      Write(33,1106)domain_pregion
	  write(33,*)loop
  endif

  do Loops = 1,Loop
     if( mype == 0) then
	    if(loops==1) write(33,*)fraci(1:kinds),qions(1:kinds),mions(1:kinds)
	    if(loops==2) write(33,*)frace(1:kinds),qions(1:kinds)*0-1.,mions(1:kinds)*0+1.
     endif

  do n=1,n_pregion
     nbuff = max(maxval(block%mi),2000)*2
     allocate(send(9, nbuff),num_recv(numberpe),request(numberpe))
     num = 0

     if(Loops ==1) then
     do m=1,mblocks
           do L = 1,block(m)%mi
              kind  = block(m)%ion(L)%kind
	          amass = mions(kind)*aion
              p     = block(m)%ion(L)%p(1)
              q     = block(m)%ion(L)%p(2)
              w     = block(m)%ion(L)%p(3)
              xyz   = pqw2xyz(0,pmin+p*deltax,qmin+q*deltay,wmin+w*deltaz)
		      x     = xyz(1)
		      y     = xyz(2)
		      z     = xyz(3)

			  if(x>= domain_pregion(n,1) .and. x< domain_pregion(n,2) .and. & 
			     y>= domain_pregion(n,3) .and. y< domain_pregion(n,4) .and. & 
			     z>= domain_pregion(n,5) .and. z< domain_pregion(n,6) )then
				 num = num +1
				 if(num >= nbuff) then
                    allocate(dummy(9,num-1))
					dummy = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(9,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
				 endif
				 send(1,num) = block(m)%ion(L)%kind
				 send(2,num) = x
				 send(3,num) = y
				 send(4,num) = z
				 send(5,num) = block(m)%ion(L)%v(1)/amass
				 send(6,num) = block(m)%ion(L)%v(2)/amass
				 send(7,num) = block(m)%ion(L)%v(3)/amass
				 send(8,num) = block(m)%ion(L)%w
			  endif	  			   
		   enddo
     enddo
     else if(Loops==2) then
     do m=1,mblocks
           do L = 1,block(m)%me
              kind  = block(m)%ele(L)%kind
              p     = block(m)%ele(L)%p(1)
              q     = block(m)%ele(L)%p(2)
              w     = block(m)%ele(L)%p(3)
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


			  if(x>= domain_pregion(n,1) .and. x< domain_pregion(n,2) .and. & 
			     y>= domain_pregion(n,3) .and. y< domain_pregion(n,4) .and. & 
			     z>= domain_pregion(n,5) .and. z< domain_pregion(n,6)  .and. &
				 kind == 1)then
				 num = num +1
				 if(num >= nbuff) then
                    allocate(dummy(9,num-1))
					dummy = send(:,1:num-1)
					deallocate(send)
					nbuff = nbuff*1.5
					allocate(send(9,nbuff))
					send(:,1:num-1) = dummy
					deallocate(dummy)
				 endif
				 send(1,num) = block(m)%ele(L)%kind
				 send(2,num) = x
				 send(3,num) = y
				 send(4,num) = z

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

				 uu       = block(m)%ele(L)%v(:)/aelectron
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma


				 vw(1)    = x2b(1,1)*vv(1) + x2b(1,2)*vv(2)+x2b(1,3)*vv(3)
				 vw(2)    = x2b(2,1)*vv(1) + x2b(2,2)*vv(2)+x2b(2,3)*vv(3)
				 vw(3)    = x2b(3,1)*vv(1) + x2b(3,2)*vv(2)+x2b(3,3)*vv(3)

				 bw(1)    = x2b(1,1)*vb1(1) + x2b(1,2)*vb1(2)+x2b(1,3)*vb1(3)
				 bw(2)    = x2b(2,1)*vb1(1) + x2b(2,2)*vb1(2)+x2b(2,3)*vb1(3)
				 bw(3)    = x2b(3,1)*vb1(1) + x2b(3,2)*vb1(2)+x2b(3,3)*vb1(3)

                 angle    = atan2theta(vw(3),vw(2)) - atan2theta(bw(3),bw(2)) 

				 send(5,num) = vw(1)
				 send(6,num) = vw(2)
				 send(7,num) = vw(3)
				 send(9,num) = angle

				 send(8,num) = block(m)%ele(L)%w
			  endif	  			   
		   enddo
     enddo
	 endif


     call mpi_allgather(num,1,mpi_integer,num_recv,1,mpi_integer,MPI_COMM_WORLD,ier)
     if(mype > 0 .and. num > 0) then
	       call mpi_isend(send(:,1:num),9*num, mpi_real,0, &
		                           mype, MPI_COMM_WORLD, request(mype),ier)

     endif
     if(mype==0) then
	    allocate(save(9,nbuff)) !nbuff
	    iplot = 0
        do i=1,num
		   iplot = iplot +1
		   istart = (/1,iplot/)
		   if(iplot >= nbuff) then
              allocate(dummy(9,iplot-1))
			  dummy = save(:,1:iplot-1)
			  deallocate(save)
			  nbuff = nbuff*1.5
			  allocate(save(9,nbuff))
			  save(:,1:iplot-1) = dummy
			  deallocate(dummy)
		   endif
           save(:,iplot) = send(:,i)
        enddo
     endif

     do L=1,numberpe-1
	    if(num_recv(L+1) > 0) then
	       if(mype==0) then
		      allocate(recv(9,num_recv(L+1)))
              call mpi_recv(recv,9*num_recv(L+1), mpi_real,L, &
		                           L, MPI_COMM_WORLD, status, ier)
		  
		      do i=1,num_recv(L+1)
		         iplot = iplot +1
		         istart = (/1,iplot/)
		         if(iplot >= nbuff) then
                    allocate(dummy(9,iplot-1))
			        dummy = save(:,1:iplot-1)
			        deallocate(save)
			        nbuff = nbuff*1.5
			        allocate(save(9,nbuff))
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
	    write(33,*)iplot
		write(33,1106)save(:,1:iplot)
		deallocate(save)
	 endif
     deallocate(send,request,num_recv)
  enddo
  enddo
  if(mype==0) close(33)

1106  format(5e12.4)
   
  deallocate(domain_pregion)
  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  return
end subroutine diagnos_particle_gyro
! ---------------------------------------------------------------------------
subroutine DiagEdistribution_gyro
  use global_parameters
  use grid_data
  use vector_functions
  use diagnos

  implicit none
  include 'mpif.h'
  integer :: nions
  type(bvector_type), dimension(:), pointer    :: bv0_,bfld,bv1_

  real(p2),pointer,dimension(:,:,:,:) :: fv,fvt
  real(p2),pointer,dimension(:) :: vpar,vper,csi
  real,pointer,dimension(:,:) :: num,numt

  integer  :: nv,ic,m,ierror,L,kind,n,nbuff
  real(p2) :: vminp,vmaxp,p,q,w,x,y,z
  real(p2) :: dcsi,dvpar,dvper
  real(8)  :: xyz(3)
  real(p2) :: vbper1(3),vbper2(3),x2b(3,3),b2x(3,3),m33(3,3),vv(3),uu(3),vb(3),vb1(3),phase,det
  integer  :: ii,jj,kk,ip,jp,kp,i,j,k
  real(p2) :: v1,v2,v3,&
              W000,W001,W010,W011,W100,W101,W110,W111,WX0,WX1,WY0,WY1,WZ0,WZ1,angle
  real(p2) :: bw(3),vw(3),vpar_,vper_ ,weight,gamma
  real     :: vtot_,ee_
  integer  :: n_E, n_alpha
  real(p2) :: delta_e,delta_alpha,E_min,E_max
  real(p2),pointer,dimension(:) :: array_e,array_alpha,volum
  real(p2),pointer,dimension(:,:,:) :: ff,fft


  save ic
  data ic/0/
  interface
      subroutine allocate_bvector(c,m)
        use constants
		use grid_data
	    type(bvector_type), dimension(:), pointer    :: c
		integer :: m
      end subroutine allocate_bvector
      subroutine dconvert(a,b)
        use constants
        real(p2),pointer,dimension(:,:,:) :: a,b
      end subroutine dconvert
      subroutine cconvert(a,b)
        use constants
        real(p2),pointer,dimension(:,:,:) :: a,b
      end subroutine cconvert
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
  end interface
!------------------------

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
  n_pregion = int(pdomain(1)+0.001)
  allocate(domain_pregion(n_pregion,6))

  L  = 1
  if(ghybrid) then
    do m=1,n_pregion
	   L = L +1
	   domain_pregion(m,1) = pdomain(L)
	   L = L +1
	   domain_pregion(m,2) = pdomain(L)
	   L = L +1
	   domain_pregion(m,3) = pdomain(L)
	   L = L +1
	   domain_pregion(m,4) = pdomain(L)
	   L = L +1
	   domain_pregion(m,5) = pdomain(L)
	   L = L +1
	   domain_pregion(m,6) = pdomain(L)
	enddo

  else
    do m=1,n_pregion
	   L = L +1
	   domain_pregion(m,1) = real(int(pdomain(L)*nx * pmove))
	   L = L +1
	   domain_pregion(m,2) = max(real(int(pdomain(L)*nx)),1.)
	   L = L +1
	   domain_pregion(m,3) = real(int(pdomain(L)*ny * qmove))
	   L = L +1
	   domain_pregion(m,4) = max(real(int(pdomain(L)*ny)),1.)
	   L = L +1
	   domain_pregion(m,5) = real(int(pdomain(L)*nz * wmove))
	   L = L +1
	   domain_pregion(m,6) = max(real(int(pdomain(L)*nz)),1.)
	enddo
  endif


  ic     = ic +1
! design a vpar [-10,10] and dv = 0.1
! f(vpar,theta,vper,region)
! vpar [-5,5]; theta[0,pi2]; vper[0,10]; region [3]
  nv      = 101
  n_e     = 101
  n_alpha = 30

  vmaxp   =  6.
  vminp   = -6.
  
  E_min   = alog10(0.0001)
  E_max   = alog10(1.e3/2.)


  dvpar   = (vmaxp-vminp)/real(nv-1)
  dcsi    = (pi*2)/real(nv-1)
  dvper   = 1.
  delta_e      = (e_max-e_min)/(n_e-1.)
  delta_alpha  = pi/2./(n_alpha-1.) 

  allocate(fv(nv,nv,6,n_pregion),fvt(nv,nv,6,n_pregion),vpar(nv),csi(nv),vper(6))
  allocate(num(6,n_pregion),numt(6,n_pregion))

  allocate(ff(n_e,n_alpha,n_pregion),fft(n_e,n_alpha,n_pregion))
  allocate(array_e(n_e),array_alpha(n_alpha),volum(n_pregion))


  do i = 1, nv
     vpar(i) = vminp + (i-1)*dvpar
	 csi(i)  = (i-1) * dcsi
  enddo
  do i=1,n_alpha
     array_alpha(i) = delta_alpha*(i-1.)
  enddo
  do i=1,n_e
     array_e(i)     = e_min + delta_e*(i-1.)
  enddo
   
  do n=1,n_pregion
     volum(n)  = (domain_pregion(n,4) - domain_pregion(n,3))*volmax0
  enddo

  fv    = 0.
  fvt   = 0.
  num   = 0
  numt  = 0
  ff    = 0.
  fft   = 0.

  do n=1,n_pregion
     do m=1,mblocks
           do L = 1,block(m)%me
              kind  = block(m)%ele(L)%kind
              p     = block(m)%ele(L)%p(1)
              q     = block(m)%ele(L)%p(2)
              w     = block(m)%ele(L)%p(3)
			  weight= block(m)%ele(L)%w
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
			  if(x>= domain_pregion(n,1) .and. x< domain_pregion(n,2) .and. & 
			     y>= domain_pregion(n,3) .and. y< domain_pregion(n,4) .and. & 
			     z>= domain_pregion(n,5) .and. z< domain_pregion(n,6)) then
!				 num = num +1
!				 send(1,num) = block(m)%ele(L)%kind
!				 send(2,num) = x
!				 send(3,num) = y
!				 send(4,num) = z

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
				 uu       = block(m)%ele(L)%v(:)/aelectron
		         gamma    =  sqrt(1+(uu(1)**2+uu(2)**2+uu(3)**2)/cspeed**2)
	             if(ifrelastic == 0) gamma   = 1.
				 vv       = uu/gamma


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

                 if(kind == 1) then 
                    i        = max(min(int((vpar_ - vminp)/dvpar + 1.),  nv-1),1)
				    wx0      = (vpar_ - vminp)/dvpar + 1. - i
				    j        = min(int(angle/dcsi + 1.),nv-1)
				    wy0      = angle/dcsi + 1. - j
				    k        = min(int(vper_/dvper +1.),6)

				    fv(i,j,k,n)     = fv(i,j,k,n)     + wx0 * wy0 * weight *frace(kind)
				    fv(i+1,j,k,n)   = fv(i+1,j,k,n)   + wx1 * wy0 * weight *frace(kind)
				    fv(i,j+1,k,n)   = fv(i,j+1,k,n)   + wx0 * wy1 * weight *frace(kind)
				    fv(i+1,j+1,k,n) = fv(i+1,j+1,k,n) + wx1 * wy1 * weight *frace(kind)
				    num(k,n)        = num(k,n) +  1. *frace(kind)
                 endif
                 vtot_ = sqrt(vpar_**2+vper_**2)
				 ee_   = (0.5*vtot_**2  )
		         if(ifrelastic > 0) ee_   = aelectron*cspeed**2*(gamma-1.)  !
				 ee_   = alog10(ee_)
				 if(ee_ < e_min) goto 222
			     angle = acos( max(min(abs(vpar_)/max(vtot_,1.e-6),1.),-1.))   ![0., pi/2]

                 i        = max(min(int((ee_ - e_min)/delta_e + 1.),  n_e-1),1)
				 wx0      = (ee_ - e_min)/delta_e + 1. - i
				 j        = min(int(angle/delta_alpha + 1.),n_alpha-1)
				 wy0      = angle/delta_alpha + 1. - j


				 ff(i,j,n)     = ff(i,j,n)     + wx0 * wy0 * frace(kind)
				 ff(i+1,j,n)   = ff(i+1,j,n)   + wx1 * wy0 * frace(kind)
				 ff(i,j+1,n)   = ff(i,j+1,n)   + wx0 * wy1 * frace(kind)
				 ff(i+1,j+1,n) = ff(i+1,j+1,n) + wx1 * wy1 * frace(kind)
222              continue
!                endif

			  endif	  			   
		   enddo
     enddo
  enddo



  CALL MPI_ALLREDUCE(num,numt,6*n_pregion,mpi_real,MPI_SUM,MPI_COMM_WORLD,IERROR)
  if(p2==8) then
     call mpi_allreduce(fv,fvt,nv*nv*6*n_pregion,mpi_double_precision,mpi_sum,MPI_COMM_WORLD, ierror)
     call mpi_allreduce(ff,fft,n_e*n_alpha*n_pregion,mpi_double_precision,mpi_sum,MPI_COMM_WORLD, ierror)
  else
     call mpi_allreduce(fv,fvt,nv*nv*6*n_pregion,mpi_real,mpi_sum,MPI_COMM_WORLD, ierror)
     call mpi_allreduce(ff,fft,n_e*n_alpha*n_pregion,mpi_real,mpi_sum,MPI_COMM_WORLD, ierror)
  endif

  do n=1,n_pregion
     do k=1,6
	    fvt(:,:,k,n) = fvt(:,:,k,n)/numt(k,n)
     enddo
	 do i=1,n_alpha
	    fft(:,i,:) = fft(:,i,:)/max(sin(array_alpha(i)),0.1)
     enddo
  enddo

  if(ic == 1) then
     IF(mype==0) then
        OPEN(11, FILE = 'fve.dat')      !, FORM="unformatted")
        write(11,*)nv,nv,6,n_pregion
        write(11,1106)vpar,csi

        OPEN(13, FILE = 'fluxe.dat')      !, FORM="unformatted")
        write(13,*)n_e,n_alpha,n_pregion
        write(13,1106)array_e,array_alpha,volum

     ENDIF
  else
      if(mype==0) OPEN(11, FILE = 'fve.dat',position='append')      !, FORM="unformatted")
      if(mype==0) OPEN(13, FILE = 'fluxe.dat',position='append')      !, FORM="unformatted")
  endif

  IF(mype==0) then
     write(11,*)stime
     write(11,1106)fvt
     close(11)
     write(13,*)stime
     write(13,1106)fft
     close(11)
  ENDIF

1106  format(5e12.4)
   
  DEALLOCATE(vpar, fv,fvt,num,numt,csi,vper)
  call allocate_bvector(bv0_,-3)
  call allocate_bvector(bv1_,-3)
  deallocate(domain_pregion)


  deallocate(ff,fft)
  deallocate(array_e,array_alpha,volum)


  return
end subroutine DiagEdistribution_gyro

!=================================================================================================
subroutine output_preciptation()  ! test particles
  use global_parameters
  use grid_data
  use global_case
! use diagnos
  implicit none
  include 'mpif.h'
  real,  pointer,dimension(:,:) :: send
  real(p2),  pointer,dimension(:,:) :: recv
  integer,  pointer,dimension(:) :: num_this,num_tot
  integer :: ic,ic0,i,iplot,m,n,iflag
  character*3 fldname(2)
  character*3 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  ::  param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,angle,vtot
  real(p2)  :: xyz(3),p,q,w,amass,wx1,wy1,wz1,wx0,wy0,wz0
  real(p2)  :: Emin,Emax,de,dalpha,uu(3),gamma,vvmin,vvmax,dv,amin,amax
  integer   :: num
  integer :: j,k,L,imc1,imc2,index,IERROR,ier
  logical file_exist
  integer :: status(mpi_status_size)

  save ic,ic0
  data ic/0/,ic0/0/

  ic       = ic +1
  param    = 0
  param(1) = utime
  param(2) = uspeed
  param(3) = ut0
  param(4) = cspeed

  if(mype==0) then 
     if(ic==1 .and. irun == 0) then
        open(11, FILE = 'preciptationi.dat')      
		write(11,1102) param    
    else
        open(11, FILE = 'preciptationi.dat',position='append')      
	endif
  endif


  CALL MPI_ALLREDUCE(mi_preciptation,num,1,MPI_integer,MPI_SUM,mpi_comm_world,IERROR)
  
  if (num ==0) then
	 if(mype==0) write(11,*) 8,num    

  else
     allocate(num_this(numberpe),num_tot(numberpe),send(8,num))
	 num_this = 0
	 num_this(mype+1) = mi_preciptation
     CALL MPI_ALLREDUCE(num_this,num_tot,size(num_this),MPI_integer,MPI_SUM,mpi_comm_world,IERROR)

     if(mype/=0) then
	     if(mi_preciptation>=1) then
           if(p2==8) call mpi_send(preciptationi(:,1:mi_preciptation), &
		                          8*mi_preciptation, mpi_double_precision, &
		                          0, 2, mpi_comm_world, ier)
           if(p2==4) call mpi_send(preciptationi(:,1:mi_preciptation), &
		                          8*mi_preciptation, mpi_real, &
		                          0, 2, mpi_comm_world, ier)
          endif
     endif
     
	 if(mype==0) then
	    mi = 0
	    send(:,mi+1:mi+mi_preciptation) = preciptationi(:,1:mi_preciptation)
		mi = mi + num_tot(1)
        do m=1,numberpe-1 
		   if(num_tot(m+1) >0) then
           allocate(recv(8,num_tot(m+1)))
           if(p2==8) call mpi_recv(recv,size(recv),mpi_double_precision, &
			              m,2,mpi_comm_world,status,ier)

           if(p2==4) call mpi_recv(recv,size(recv),mpi_real, &
			              m,2,mpi_comm_world,status,ier)
           
	       send(:,mi+1:mi+num_tot(m+1)) = recv(:,:)
		   mi = mi + num_tot(m+1)
		   deallocate(recv)
		   endif
        enddo
	 endif

	 if(mype==0) write(11,*) 8,num    
	 if(mype==0) write(11,1102) send    

     deallocate(num_this,num_tot,send)
	 mi_preciptation = 0
  endif

  if(mype==0) close(11)

  if(mod(simtype(1),10) >= 6) return

! for electrons -------------------------------------------------
  
  if(mype==0) then 
     if(ic==1 .and. irun == 0) then
        open(11, FILE = 'preciptatione.dat')      
		write(11,1102) param    
    else
        open(11, FILE = 'preciptatione.dat',position='append')      
	endif
  endif


  CALL MPI_ALLREDUCE(me_preciptation,num,1,MPI_integer,MPI_SUM,mpi_comm_world,IERROR)
  
  if (num ==0) then
	 if(mype==0) write(11,*) 8,num    

  else
     allocate(num_this(numberpe),num_tot(numberpe),send(8,num))
	 num_this = 0
	 num_this(mype+1) = me_preciptation
     CALL MPI_ALLREDUCE(num_this,num_tot,size(num_this),MPI_integer,MPI_SUM,mpi_comm_world,IERROR)

     if(mype/=0) then
	     if(me_preciptation>=1) then
           if(p2==8) call mpi_send(preciptatione(:,1:me_preciptation), &
		                          8*me_preciptation, mpi_double_precision, &
		                          0, 2, mpi_comm_world, ier)
           if(p2==4) call mpi_send(preciptatione(:,1:me_preciptation), &
		                          8*me_preciptation, mpi_real, &
		                          0, 2, mpi_comm_world, ier)
          endif
     endif
     
	 if(mype==0) then
	    mi = 0
	    send(:,mi+1:mi+me_preciptation) = preciptatione(:,1:me_preciptation)
		mi = mi + num_tot(1)
        do m=1,numberpe-1 
		   if(num_tot(m+1) >0) then
           allocate(recv(8,num_tot(m+1)))
           if(p2==8) call mpi_recv(recv,size(recv),mpi_double_precision, &
			              m,2,mpi_comm_world,status,ier)

           if(p2==4) call mpi_recv(recv,size(recv),mpi_real, &
			              m,2,mpi_comm_world,status,ier)
           
	       send(:,mi+1:mi+num_tot(m+1)) = recv(:,:)
		   mi = mi + num_tot(m+1)
		   deallocate(recv)
		   endif
        enddo
	 endif

	 if(mype==0) write(11,*) 8,num    
	 if(mype==0) write(11,1102) send    

     deallocate(num_this,num_tot,send)
	 me_preciptation = 0
  endif

  if(mype==0) close(11)


  1102  format(5e16.5)

  return
end subroutine output_preciptation

!=================================================================================================
subroutine output_preciptationA()  ! test particles
  use global_parameters
  use grid_data
  use global_case
! use diagnos
  implicit none
  include 'mpif.h'
  real,  pointer,dimension(:,:) :: send
  real(p2),  pointer,dimension(:,:) :: recv
  integer,  pointer,dimension(:) :: num_this,num_tot
  integer :: ic(8),ic0,i,iplot,m,n,iflag,ifile
  character*3 fldname(2)
  character*3 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  ::  param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,angle,vtot
  real(p2)  :: xyz(3),p,q,w,amass,wx1,wy1,wz1,wx0,wy0,wz0
  real(p2)  :: Emin,Emax,de,dalpha,uu(3),gamma,vvmin,vvmax,dv,amin,amax
  integer   :: num
  integer :: j,k,L,imc1,imc2,index,IERROR,ier
  logical file_exist
  integer :: status(mpi_status_size)

  save ic
  data ic/0,0,0,0,0,0,0,0/

  ic        = ic +1
  param     = 0
  param(1)  = utime
  param(2)  = uspeed
  param(3)  = ut0
  param(4)  = cspeed
  param(11) = micell
  param(12) = mecell
  param(21:21+kinds-1) = mions(1:kinds)
  param(31:31+kinds-1) = meles(1:kinds)

  do ifile =1,8
  if (ifpreciptation(ifile) > 0 ) then
     if(mype==0) then 
        if(ifile ==1) filename ='preciptation_i.dat'
        if(ifile ==2) filename ='preciptation_e.dat'
        if(ifile ==3) filename ='preciptation_O.dat'
        if(ifile ==4) filename ='preciptation_He.dat'
        if(ifile ==5) filename ='preciptation_H.dat'
        if(ifile ==6) filename ='preciptation_teste1.dat'
        if(ifile ==7) filename ='preciptation_teste2.dat'
        if(ifile ==8) filename ='preciptation_teste3.dat'
     inquire(file=filename,exist=file_exist)	    
     if((ic(ifile)==1 .and. irun == 0) .or. (.not. file_exist)) then
        open(11, FILE = filename) 
		if(ifile ==1) param(20) = aion     
		if(ifile ==2) param(20) = aelectron    
		write(11,1102) param    
    else
        open(11, FILE = filename,position='append')      
	endif
  endif

  CALL MPI_ALLREDUCE(loss(ifile)%mi,num,1,MPI_integer,MPI_SUM,mpi_comm_world,IERROR)
  
  if (num ==0) then
	 if(mype==0) write(11,*) 8,num    !
  else
     allocate(num_this(numberpe),num_tot(numberpe),send(8,num))
	 num_this = 0
	 num_this(mype+1) = loss(ifile)%mi
     CALL MPI_ALLREDUCE(num_this,num_tot,size(num_this),MPI_integer,MPI_SUM,mpi_comm_world,IERROR)

     if(mype/=0) then
	     if(loss(ifile)%mi>=1) then
           if(p2==8) call mpi_send(loss(ifile)%preciptation(:,1:loss(ifile)%mi), &
		                          8*loss(ifile)%mi, mpi_double_precision, &
		                          0, 2, mpi_comm_world, ier)
           if(p2==4) call mpi_send(loss(ifile)%preciptation(:,1:loss(ifile)%mi), &
		                          8*loss(ifile)%mi, mpi_real, &
		                          0, 2, mpi_comm_world, ier)
          endif
     endif
      
	 if(mype==0) then
	    mi = 0
	    send(:,mi+1:mi+loss(ifile)%mi) = loss(ifile)%preciptation(:,1:loss(ifile)%mi)
		mi = mi + num_tot(1)
        do m=1,numberpe-1 
		   if(num_tot(m+1) >0) then
           allocate(recv(8,num_tot(m+1)))
           if(p2==8) call mpi_recv(recv,size(recv),mpi_double_precision, &
			              m,2,mpi_comm_world,status,ier)

           if(p2==4) call mpi_recv(recv,size(recv),mpi_real, &
			              m,2,mpi_comm_world,status,ier)
           
	       send(:,mi+1:mi+num_tot(m+1)) = recv(:,:)
		   mi = mi + num_tot(m+1)
		   deallocate(recv)
		   endif
        enddo
	 endif
	 if(mype==0) write(11,*) 8,num    
	 if(mype==0) write(11,1102) send    

     deallocate(num_this,num_tot,send)
	 loss(ifile)%mi = 0
  endif

  if(mype==0) close(11)

  endif
  enddo

  1102  format(5e16.5)

  return
end subroutine output_preciptationA

!=================================================================================================
subroutine diagfine_fld_position()  ! test particles
  use global_parameters
  use grid_data
  use global_case
  use diagnos
! use diagnos
  implicit none
  include 'mpif.h'
  real,  pointer,dimension(:,:) :: ebfld
  integer,  pointer,dimension(:,:) :: ijk,isd
  integer,  pointer,dimension(:)   :: id_cpu,id_cpu0
  integer :: ic,ic0,i,iplot,m,n,iflag
  character*3 fldname(2)
  character*2 ct
  character (len = *), parameter :: UNITS = "units" 
  character*100 filename
  integer :: ncid,dimid(9),dimids2(2),icount2(2),varid(30), &
             dimids(3),icount(3),istart(3),icount1(1),istart1(1)
  real  ::  param(80)
  real(p2)::minx(3),maxx(3),x,y,z,r,factor,dmc1,dmc2,&
            bx,by,bz,vx,vy,vz,vpar,vper,ee,angle,vtot
  real(p2)  :: xyz(3),p,q,w,amass,wx1,wy1,wz1,wx0,wy0,wz0
  real(p2)  :: Emin,Emax,de,dalpha,uu(3),gamma,vvmin,vvmax,dv,amin,amax
  real      :: fsave(4,6,500),fsave_t(500)
  integer   :: num,nposition,icc
  integer :: j,k,L,imc1,imc2,index,IERROR,ier
  logical file_exist
  integer :: status(mpi_status_size)

  save ic,ic0,id_cpu,icc,fsave,fsave_t
  data ic/0/,ic0/0/,icc/0/

  if(kmype /= 0) return  

  if(istep < positionData(1)) return

  nposition = positionData(2)
  allocate(ebfld(nposition,6),ijk(nposition,3))
  L = 2
  do m=1,nposition
	 L = L +1
     ijk(m,1) = positionData(L)
	 L = L +1
     ijk(m,2) = positionData(L)
	 L = L +1
     ijk(m,3) = positionData(L)
  enddo

  ic       = ic +1
  icc      = icc +1
  param    = 0
  param(1) = utime
  param(2) = uspeed
  param(3) = ut0
  param(4) = cspeed
  param(5) = uEfield
  param(6) = uBfield
  param(7) = aion/aelectron  !mass ratio
  param(10)= nposition
  fsave_t(icc) = stime

  if(ic == 1) allocate(id_cpu(nposition),id_cpu0(nposition))  
  ebfld = 0.
  
  if(ic == 1) then
     id_cpu0  = 0
     do m=1,mblocks
        do i= block(m)%i0,block(m)%i1
           do j= block(m)%j0,block(m)%j1
              do k= block(m)%k0,block(m)%k1
		         do n=1,nposition
                    if(i==ijk(n,1) .and. j==ijk(n,2) .and. k==ijk(n,3) ) then
                       ebfld(n,1:3) = bv1(m)%vector(:,i,j,k) * block(m)%node%h(4,:,i,j,k)
                       ebfld(n,4:6) = ev1(m)%vector(:,i,j,k) / block(m)%node%h(3,:,i,j,k)
					   id_cpu0(n) = kmype1+1
					   fsave(n,:,icc) = ebfld(n,:) 
				    endif
			     enddo
              enddo
           enddo
        enddo
     enddo
  else
     do n= 1, nposition
        if(kmype1 == id_cpu(n)-1) then
        do m=1,mblocks
           do i= block(m)%i0,block(m)%i1
              do j= block(m)%j0,block(m)%j1
                 do k= block(m)%k0,block(m)%k1
                    if(i==ijk(n,1) .and. j==ijk(n,2) .and. k==ijk(n,3) ) then
                       ebfld(n,1:3) = bv1(m)%vector(:,i,j,k) * block(m)%node%h(4,:,i,j,k)
                       ebfld(n,4:6) = ev1(m)%vector(:,i,j,k) / block(m)%node%h(3,:,i,j,k)
					   fsave(n,:,icc) = ebfld(n,:) 
                    endif
                 enddo
              enddo
            enddo
		enddo
		endif
	 enddo
  endif


  if(ic==1)call mpi_allreduce(id_cpu0,id_cpu,nposition,mpi_integer,mpi_sum,myComm1,ier)
  if(ic==1) deallocate(id_cpu0) 

  if(mod(istep,500) == 0) then
  do n=1,nposition
     write(ct,'(I2.2)')n
     if(kmype1 == id_cpu(n)-1) then 
        inquire(file='fldt_multiple_position'//ct//'.dat',exist=file_exist)	    
        if((ic==1 .and. irun == 0) .or. .not. file_exist) then


		   param(21:23)  =  positionData(3+ (n-1)*3:3+ n*3-1 )
           open(11, FILE = 'fldt_multiple_position'//ct//'.dat')      
		   write(11,*) 80,1,6   
		   write(11,1102) param    
        else
           open(11, FILE = 'fldt_multiple_position'//ct//'.dat',position='append')      
	    endif

        do i = 1, icc
	  	   write(11,*) fsave_t(i)         !stime  
           write(11,1102) fsave(n,:,i)  !ebfld(n,:)
        enddo
		close(11)
     endif
  enddo  
     icc = 0
  endif

  deallocate(ebfld,ijk)


  1102  format(5e16.5)

  return
end subroutine diagfine_fld_position




