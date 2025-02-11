!######################################################################################
module Lib
!######################################################################################
	use common_param
	use common_var
    use errorModule
	implicit none
	
Contains

!******************************************************************************
subroutine forward_propagate(N,lambda_b,lambda_p,lambda_s,WB,Gf,zp)
!******************************************************************************
    integer(4), intent(in) :: N, zp
    real(8), intent(in) :: lambda_b, lambda_p, lambda_s
    real(8), intent(in) :: WB(0:n_layer+1)
    real(8), intent(inout) :: Gf(0:n_layer+1,0:N,-1:1)
    integer(4) s, z

	Gf(zp,0,1) = WB(1)
	Gf(zp,0,0) = WB(1)
	Gf(zp,0,-1) = WB(1)
	    
    if (zp.eq.0) then
	    Gf(zp,0,1) = WB(1)
	    Gf(zp,0,0) = 0.0
	    Gf(zp,0,-1) = 0.0
    endif
    
    if (zp.eq.n_layer) then
	    Gf(zp,0,1) = 0.0
	    Gf(zp,0,0) = 0.0
	    Gf(zp,0,-1) = WB(1)
    endif

    do s = 1, N
        do z = 1, n_layer
            Gf(z,s,-1) = WB(z)*(lambda_b*Gf(z-1,s-1,1)  &
               & + 4.0*lambda_p*Gf(z,s-1,0)  &
               & + lambda_s*Gf(z+1,s-1,-1))
            
               !!!
               Gf(z,s,0) = WB(z)*(lambda_p*Gf(z-1,s-1,1)  &
               & + (2.0*lambda_p + lambda_s + lambda_b)*Gf(z,s-1,0)  &
               & + lambda_p*Gf(z+1,s-1,-1))
            
               !!!
               Gf(z,s,1) = WB(z)*(lambda_s*Gf(z-1,s-1,1)  &
               & + 4.0*lambda_p*Gf(z,s-1,0)  &
               & + lambda_b*Gf(z+1,s-1,-1))
            
           enddo
        enddo


return
end subroutine forward_propagate
!******************************************************************************

!******************************************************************************
subroutine back_propagate(N,lambda_b,lambda_p,lambda_s,WB,Gb)
!******************************************************************************
    integer(4), intent(in) :: N
    real(8), intent(in) :: lambda_b, lambda_p, lambda_s
    real(8), intent(in) :: WB(0:n_layer+1)
    real(8), intent(inout) :: Gb(0:n_layer+1,1:N,-1:1)
    integer(4) s, z
    
     Gb(1:n_layer,N,1) = WB(1:n_layer)
     Gb(1:n_layer,N,0)  = WB(1:n_layer)
     Gb(1:n_layer,N,-1) = WB(1:n_layer)
    
    do s = 2, N
        do z = 1, n_layer
            Gb(z,N-s+1,-1) = WB(z) * ( &
               &   lambda_b*Gb(z+1,N-s+2,1) &
               & + 4.0*lambda_p*Gb(z+1,N-s+2,0) &
               & + lambda_s*Gb(z+1,N-s+2,-1) )
               !!!
               Gb(z,N-s+1,0) = WB(z) * ( &
               &   lambda_p*Gb(z,N-s+2,1) &
               & + (2.0*lambda_p + lambda_s + lambda_b)*Gb(z,N-s+2,0) &
               & + lambda_p*Gb(z,N-s+2,-1) )
               !!!
               Gb(z,N-s+1,1) = WB(z) * ( &
               &   lambda_s*Gb(z-1,N-s+2,1) &
               & + 4.0*lambda_p*Gb(z-1,N-s+2,0) &
               & + lambda_b*Gb(z-1,N-s+2,-1) )
           enddo
   enddo
        
return
end subroutine back_propagate
!******************************************************************************
!**************************************************************************************
subroutine calc_trans_probabilities(p,lambda)
!**************************************************************************************
    implicit none
    real(8), intent(in) :: p 
    real(8), intent(out) :: lambda(-1:1) 
    real(8) exp_K 
    ! p = (1 + <cos(ang)>) / (1 - <cos(ang)>),                                      (1)
    ! where <cos(ang)>  is angle between two segment-vectors
    ! <cos(ang)> = cos(pi)*lambda(-1) + cos(pi/2)*lambda(0) + cos(0)*lambda(+1)     (2)
    ! U(ang) = K * (1 - cos(ang)) is angle potential                                (3)
    ! lambda(-1) = A*exp(-U(pi))                                                    (4)
    ! lambda( 0) = A*exp(-U(pi/2))                                                  (5)
    ! lambda(+1) = A*exp(-U(0))                                                     (6)
    ! from (1-6) obtain:

    exp_K = sqrt((1.0-p)**2+p) + p - 1.0
    
    lambda(0)  = 1.0/(exp_K + 1.0/exp_K + 4.0)
    lambda(-1) = lambda(0) / exp_K
    lambda(+1) = lambda(0) * exp_K
    lambda(0)  = (1.0 - lambda(-1) - lambda(+1)) / 4.0
    
return
end subroutine calc_trans_probabilities
!**************************************************************************************

!**************************************************************************************
subroutine calc_U_Flory(chi, phi, phi_bulk, U)
!**************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), intent(in) :: chi, phi_bulk, phi(0:n_layer+1)
    real(8), intent(inout) :: U(0:n_layer+1)
    integer(4) z
    real(8) mean_phi

    mean_phi = 0.0      
    do z = 1, n_layer
        mean_phi = (phi(z-1)+4.0*phi(z)+phi(z+1))/6.0
        U(z) = U(z) + chi*(mean_phi-phi_bulk)
    enddo
    
return
end subroutine calc_U_Flory
!**************************************************************************************

!**************************************************************************************
subroutine calc_U_int(chi, phi, phi_bulk, phi_base, U_int)
!**************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), intent(in) :: chi, phi_bulk, phi(0:n_layer+1), phi_base(0:n_layer+1)
    real(8), intent(inout) :: U_int
    integer(4) z
    real(8) mean_phi

    mean_phi = 0.0      
    do z = 1, n_layer
        mean_phi = (phi(z-1)+4.0*phi(z)+phi(z+1))/6.0
        U_int = U_int + chi*(mean_phi-phi_bulk)*phi_base(z)
    enddo
    
return
end subroutine calc_U_int
!**************************************************************************************

!**************************************************************************************
subroutine check_deviation()
!**************************************************************************************

    eta = eta / golden_ratio
    deviation = 1000000000.0
    deviation_old = deviation + 1.0
    open(unit=n_infor,file=name_infor,&
	&    form='formatted',access='append',status='unknown')
    write(n_infor,'(3x,a,I0,2x,a,E13.7,2x,a,E13.7)') &
        &         'iter = ', iter, 'dev = ', deviation, 'eta = ', eta	
    close(n_infor)
    iter = 0.0
    if (eta.lt.tolerance) call ERROR(17)
    
    return
end subroutine check_deviation
!**************************************************************************************

end module Lib

