!######################################################################################
module errorModule
!######################################################################################
	use common_param
	use common_var
	implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
Contains

!**************************************************************************************
subroutine ERROR(n)
!**************************************************************************************
    integer(4) :: n
	
	open(unit=n_infor,file=name_infor,&
	&    form='formatted',access='append',status='unknown')
	
    select case (n)
        case(1)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. Input file is not found !'
        case(2)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. Input file is not correct !'
        case(3)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. N1 < 1 !'
        case(4)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. sigma1 < 0 !'
        case(5)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. p1 < 1 !'
        case(6)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. N2 < 1 !'
        case(7)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. sigma2 < 0 !'
        case(8)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. p2 < 1 !'
        case(9)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. sigma1+sigma2 >= 1 !'
        case(11)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. eta is not found !'
        case(12)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. chi1 is not found !'
        case(13)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. chi2 is not found !'
        case(14)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. chi12 is not found !'
        case(15)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. nfree is not found !'
        case(16)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. swpro is not found !'
        case(17)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. eta < tolerance !'
        case(18)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. ksi is not found !'
        case(19)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. zp1 is not found !'
        case(20)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. zp2 is not found !'
        case(21)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. n_layer is not found !'
        case(22)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. left_w is not found !'
        case(23)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. right_w is not found !'
        case(24)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. zp1 value is invalid !'
        case(25)
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. zp2 value is invalid !'
        case default
            write(n_infor,'(3x,a,I0,a)') 'ERROR ',n,'. Unknown error !'
    end select
    
    write(n_infor,'(a)') '============================================================'
	write(n_infor,'(3x,a)') 'Sorry, something went wrong there...' 
	write(n_infor,'(a)') '============================================================'
	write(n_infor,'(3x,a,F16.5)') 'Lead time (sec): ', lead_time
	write(n_infor,'(a)') '============================================================'
    close(n_infor)

	call Close_all_files()
	call Deallocate_all_variables()
    stop
end subroutine ERROR
!**************************************************************************************

!**************************************************************************************
subroutine Close_all_files()
!**************************************************************************************
	integer(4) :: i
	logical    :: ISopen
	do i = 1, size(list_file)
		INQUIRE(unit = list_file(i),opened = ISopen)
		if(ISopen) close(list_file(i))
	enddo
end subroutine Close_all_files
!**************************************************************************************


!######################################################################################
end module errorModule
!######################################################################################
