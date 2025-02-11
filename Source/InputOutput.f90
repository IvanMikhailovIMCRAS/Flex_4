!######################################################################################
module InputOutput
!######################################################################################
    use common_param
	use common_var
    use errorModule
    implicit none
		
Contains

!**************************************************************************************
subroutine dataInput()
!**************************************************************************************

    open(n_input, file = name_input, status = 'old', iostat=ioer)
        if (ioer.ne.0) call ERROR(1)
        read(n_input,*,iostat=ioer) N1
        if (ioer.ne.0) call ERROR(2)
        if (N1.lt.1) call ERROR(3)
        read(n_input,*,iostat=ioer) sigma1
        if (ioer.ne.0) call ERROR(2)
        if (sigma1.lt.0.0) call ERROR(4)
        read(n_input,*,iostat=ioer) p1
        if (ioer.ne.0) call ERROR(2)
        if (p1.lt.1.0) call ERROR(5)
        read(n_input,*,iostat=ioer) chi1
        if (ioer.ne.0) call ERROR(12)
        read(n_input,*,iostat=ioer) zp1
        if (ioer.ne.0) call ERROR(19)
        read(n_input,*,iostat=ioer) N2
        if (ioer.ne.0) call ERROR(2)
        if (N2.lt.1) call ERROR(6)
        read(n_input,*,iostat=ioer) sigma2
        if (ioer.ne.0) call ERROR(2)
        if (sigma2.lt.0.0) call ERROR(7)
        read(n_input,*,iostat=ioer) p2
        if (ioer.ne.0) call ERROR(2)
        if (p2.lt.1.0) call ERROR(8)
        read(n_input,*,iostat=ioer) chi2
        if (ioer.ne.0) call ERROR(13)
        read(n_input,*,iostat=ioer) zp2
        if (ioer.ne.0) call ERROR(20)
        read(n_input,*,iostat=ioer) chi12
        if (ioer.ne.0) call ERROR(14)
        read(n_input,*,iostat=ioer) n_layer
        if (ioer.ne.0) call ERROR(21)
        read(n_input,*,iostat=ioer) left_w
        if (ioer.ne.0) call ERROR(22)
        read(n_input,*,iostat=ioer) right_w
        if (ioer.ne.0) call ERROR(23)
        if (sigma1+sigma2.ge.1.0) call ERROR(9)
        read(n_input,*,iostat=ioer) eta
        if (ioer.ne.0) call ERROR(11)
        read(n_input,*,iostat=ioer) ksi
        if (ioer.ne.0) call ERROR(17)
        read(n_input,*,iostat=ioer) nfree
        if (ioer.ne.0) call ERROR(15)
        read(n_input,*,iostat=ioer) swpro
        if (ioer.ne.0) call ERROR(16)
    close(n_input)
    if (zp1.gt.n_layer+1.or.zp1.lt.0) call ERROR(24)
    if (zp2.gt.n_layer+1.or.zp2.lt.0) call ERROR(25)
    return
end subroutine dataInput
!**************************************************************************************

!**************************************************************************************
subroutine read_initial_guess()
!**************************************************************************************
	implicit none
	integer(4) :: i
	logical :: ISexist
	real(8) :: value(1:4)
	
	INQUIRE(file = name_ig, exist = ISexist)
	if (ISexist) then
		open(n_ig, file = name_ig, iostat = ioer)
	else
		return
	endif
	i = 0
	if (ioer.eq.0) then
		read(n_ig,*,iostat=ioer) eta
	endif
	do while (ioer.eq.0)
		if (i.gt.n_layer+1) then
			exit
		endif
		read(n_ig,*,iostat=ioer) value
		if (ioer.ne.0) exit
		alpha(i) = value(1)
		U1(i) = value(2) 
		U2(i) = value(3)
		Us(i) = value(4)
		i = i + 1
	enddo
	close(n_ig)
		
	return
end subroutine read_initial_guess
!**************************************************************************************

!**************************************************************************************
subroutine print_initial_guess()
!**************************************************************************************
	implicit none
	integer(4) :: i
	
	open(n_ig, file = name_ig)
	write(n_ig,*) eta
	do i = 0, n_layer+1
		write(n_ig,*) alpha(i), U1(i), U2(i), Us(i)
	enddo
	close(n_ig) 
	
	return
end subroutine print_initial_guess
!**************************************************************************************

!**************************************************************************************
subroutine print_phi_profiles()
!**************************************************************************************
	implicit none
	integer(4) z
	real(8) :: phi_p1, phi_p2, ne1, ne2
	
	open(n_pro, file=name_pro)
        write(n_pro,'(8x,a,8x,a,16x,a,15x,a,15x,a)') 'z','phi1','phi2','end1','end2'
        do z = 1, n_layer
        	if (phi_pol1(z).lt.1e-99_8) then
        		phi_p1 = 0.0
        	else
        		phi_p1 = phi_pol1(z)
        	endif
        	if (phi_pol2(z).lt.1e-99_8) then
        		phi_p2 = 0.0
        	else
        		phi_p2 = phi_pol2(z)
        	endif
        	if (n_end1(z).lt.1e-99_8) then
        		ne1 = 0.0
        	else
        		ne1 = n_end1(z)
        	endif
        	if (n_end2(z).lt.1e-99_8) then
        		ne2 = 0.0
        	else
        		ne2 = n_end2(z)
        	endif
            write(n_pro,'(I9,E20.10,E20.10,E20.10,E20.10)') &
                    & z, phi_p1, phi_p2, ne1, ne2
        enddo
    close(n_pro)
end subroutine print_phi_profiles
!**************************************************************************************

!**************************************************************************************
subroutine print_data(m1, m2, F)
!**************************************************************************************
	implicit none
	real(8), intent(in) :: m1, m2, F
	open(n_data, file=name_data)
    write(n_data,'(8x,a,16x,a,15x,a)') 'H_1','H_2','F'  
    write(n_data,'(E20.10,E20.10,E20.10)') m1, m2, F
    close(n_data)
end subroutine print_data
!**************************************************************************************

!**************************************************************************************
subroutine Print_infor_title()
!**************************************************************************************
	implicit none
	
	open(unit=n_infor,file=name_infor,form='formatted',status='replace')
	
	write(n_infor,'(a)') '============================================================'
	write(n_infor,'(9x,a,1x,a,4x,a,1x,a)') 'PROGRAM',name_program,'ver',version 
	write(n_infor,'(a)') '============================================================'
	write(n_infor,'(a,2x,a)') '(c)', copyright
	write(n_infor,'(5x,a)')          authors
	write(n_infor,'(5x,a,a,a)')     date,', ',year
	write(n_infor,'(a)') '============================================================'
	close(n_infor)
	
	return
end subroutine Print_infor_title
!**************************************************************************************

!**************************************************************************************
subroutine Print_infor_finish()
!**************************************************************************************
	implicit none
	
	open(unit=n_infor,file=name_infor,&
	&    form='formatted',access='append',status='unknown')
	
	write(n_infor,'(3x,a,I0,2x,a,E13.7,2x,a,E13.7)') &
        &         'iter = ', iter, 'dev = ', deviation, 'eta = ', eta	
	write(n_infor,'(a)') '============================================================'
	write(n_infor,'(3x,a)') 'Awesome! Everything was working as it should!' 
	write(n_infor,'(a)') '============================================================'
	write(n_infor,'(3x,a,F16.5)') 'Lead time (sec): ', lead_time
	write(n_infor,'(a)') '============================================================'
	close(n_infor)
	
	return
end subroutine Print_infor_finish
!**************************************************************************************

!######################################################################################
end module InputOutput
!######################################################################################
