!######################################################################################
!! DESCRIPTION: SCF-simulation of two semi-flexible chains grafted onto 
!!              a flat surface under selective solvent condition
!! AUTHOR : I.V. Mikhailov, Institute of Macromolecular Compounds IMC RAS (c)
!######################################################################################
program FLEX
	use common_param
	use common_var
	use InputOutput
	use Lib
	use errorModule
    implicit none
    real(8) :: time(1:3) 
    
    call cpu_time(time(1)); call cpu_time(time(2))	
    call Print_infor_title()
    ! reading all parameters from input-file
    call dataInput()
    ! allocate memory for arrays
    call Allocate_all_variables()
    ! reading initial guess for Lagrange field if it is possible
	call read_initial_guess() 
    ! running the program engine
    call main()
    ! printing Lagrange field into external file "initial_guess.in"
    call print_initial_guess()

	call Deallocate_all_variables()
	
	call cpu_time(time(3))
	lead_time = time(3) - time(2) - (time(2) - time(1))	
	call Print_infor_finish()	

    stop
!######################################################################################
end program FLEX
!######################################################################################

!**************************************************************************************
subroutine main()
!**************************************************************************************
	use common_param
	use common_var
	use InputOutput
	use Lib
	use errorModule
	implicit none
    integer(4) s, z
    real(8) F, Uint1, Uint2, Uints
    real(8) sum1, sum2, m1_1, m1_2
    real(8) sum_n_r1, sum_n_r2
    real(8) lambda(-1:+1)
    
    !!!! transition probabilities
    lambda(:) = 0.0
    call calc_trans_probabilities(p1,lambda)  
    lambda_1p = lambda(0)  ! perpendicular kink
    lambda_1s = lambda(+1) ! straight transition (open angle)
    lambda_1b = lambda(-1) ! backfold fracture
    !!!!
    lambda(:) = 0.0
    call calc_trans_probabilities(p2,lambda)
    !!!
    lambda_2p = lambda(0)  ! perpendicular kink
    lambda_2s = lambda(+1) ! straight transition (open angle)
    lambda_2b = lambda(-1) ! backfold fracture
    !!!!
    iter = 0
    deviation = 1000000000.0
    counter_free = 0
    phi_solv(n_layer+1) = 1.0_8
!!!!!!!!!!!!!!!!!!!! ENGINE OF PROGRAM  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do while (deviation.gt.tolerance)
        iter = iter + 1
        if (iter.ge.max_iter) then
            write(*,'(a,I0)') 'Iteration steps > ', max_iter
            stop
        endif
        !!!! Boltzmann's weights
        WB1(:)=exp(-alpha(:)-U1(:))
        WB2(:)=exp(-alpha(:)-U2(:))
        WBs(1:n_layer)=exp(-alpha(1:n_layer)-Us(1:n_layer))
        !!!! Boundary conditions
        if (left_w.ne.0) then
           WB1(0) = WB1(1)
           WB2(0) = WB2(1)
           WBs(0) = WBs(1)
        endif
        if (right_w.ne.0) then
           WB1(n_layer+1) = WB1(n_layer)
           WB2(n_layer+1) = WB2(n_layer)
           WBs(n_layer+1) = WBs(n_layer)
        endif        
        !!!! Volume fraction profile of solvent
        phi_solv(1:n_layer) = WBs(1:n_layer)
		!!!! zero-initialization of propogators
        Gf1(:,:,:) = 0.0; Gb1(:,:,:) = 0.0; Gf2(:,:,:) = 0.0; Gb2(:,:,:) = 0.0
        !!!! Now we iterate over all Green's functions (propagators) 
        call forward_propagate(N1,lambda_1b,lambda_1p,lambda_1s,WB1,Gf1,zp1)
        call back_propagate(N1,lambda_1b,lambda_1p,lambda_1s,WB1,Gb1)
        call forward_propagate(N2,lambda_2b,lambda_2p,lambda_2s,WB2,Gf2,zp2) 
        call back_propagate(N2,lambda_2b,lambda_2p,lambda_2s,WB2,Gb2)
        ! Here we calculate profiles phi1 and phi2 according to composition law
        phi_pol1(:) = 0.0; phi1(:,:) = 0.0; phi_pol2(:) = 0.0; phi2(:,:) = 0.0
        do s = 1, N1
            do z = 1, n_layer
                phi1(z,0)  = phi1(z,0)  + Gf1(z,s,0)*Gb1(z,s,0)/WB1(z)*4.0
                phi1(z,-1) = phi1(z,-1) + Gf1(z,s,-1)*Gb1(z,s,1)/WB1(z)
                phi1(z,1)  = phi1(z,1)  + Gf1(z,s,1)*Gb1(z,s,-1)/WB1(z)
            enddo
        enddo
        !!!!
        do s = 1, N2
            do z = 1, n_layer
                phi2(z,0)  = phi2(z,0)  + Gf2(z,s,0)*Gb2(z,s,0)/WB2(z)*4.0
                phi2(z,-1) = phi2(z,-1) + Gf2(z,s,-1)*Gb2(z,s,1)/WB2(z)
                phi2(z,1)  = phi2(z,1)  + Gf2(z,s,1)*Gb2(z,s,-1)/WB2(z)
            enddo
        enddo
        !!!!
        phi_pol1(:) = phi1(:,-1) + phi1(:,0) + phi1(:,1)
        part_fun1 = sum(phi_pol1(1:n_layer))
        phi_pol2(:) = phi2(:,-1) + phi2(:,0) + phi2(:,1)
        part_fun2 = sum(phi_pol2(1:n_layer))
        !!!!
        phi1(:,-1) = sigma1*phi1(:,-1)*dble(N1)/part_fun1
        phi1(:,0) = sigma1*phi1(:,0)*dble(N1)/part_fun1
        phi1(:,1) = sigma1*phi1(:,1)*dble(N1)/part_fun1
        phi_pol1(:) = phi1(:,-1) + phi1(:,0) + phi1(:,1)
        !!!!
        phi2(:,-1) = sigma2*phi2(:,-1)*dble(N2)/part_fun2
        phi2(:,0) = sigma2*phi2(:,0)*dble(N2)/part_fun2
        phi2(:,1) = sigma2*phi2(:,1)*dble(N2)/part_fun2      
        phi_pol2(:) = phi2(:,-1) + phi2(:,0) + phi2(:,1)
        !!!! Flory-Huggins interactions 
        W1(:) = U1(:); W2(:) = U2(:); Ws(:) = Us(:)
        U1(:) = 0.0;   U2(:) = 0.0;   Us(:) = 0.0   
        call calc_U_Flory(chi1, phi_solv, 1.0_8, U1)
        call calc_U_Flory(chi2, phi_solv, 1.0_8, U2)
        call calc_U_Flory(chi1, phi_pol1, 0.0_8, Us)
        call calc_U_Flory(chi2, phi_pol2, 0.0_8, Us)       
        call calc_U_Flory(chi12, phi_pol1, 0.0_8, U2)
        call calc_U_Flory(chi12, phi_pol2, 0.0_8, U1)
      
        U1(:) = ksi*U1(:) + (1.0 - ksi)*W1(:)
        U2(:) = ksi*U2(:) + (1.0 - ksi)*W2(:)
        Us(:) = ksi*Us(:) + (1.0 - ksi)*Ws(:)
        !!!! Gradient descent for Lagrange field
        alpha(1:n_layer) = alpha(1:n_layer) + eta*( &
        &            phi1(1:n_layer,-1) + phi1(1:n_layer,0) + phi1(1:n_layer,1) + &
        &            phi2(1:n_layer,-1) + phi2(1:n_layer,0) + phi2(1:n_layer,1) + &
        &            phi_solv(1:n_layer) - 1.0)
        !!!! estimating the deviations from the incompressibility conditions, etc.
        deviation_old = deviation
        deviation = sum((phi1(1:n_layer,-1) + phi1(1:n_layer,0) + phi1(1:n_layer,1) + &
        &                phi2(1:n_layer,-1) + phi2(1:n_layer,0) + phi2(1:n_layer,1) + &
        &                phi_solv(1:n_layer) - 1.0)**2 + &
        &                (U1(1:n_layer)-W1(1:n_layer))**2 + &
        &                (U2(1:n_layer)-W2(1:n_layer))**2 + &
        &                (Us(1:n_layer)-Ws(1:n_layer))**2)
        deviation = sqrt(deviation)
        !!!! automatic selection of step size $eta$ and $ksi$
        if (deviation.ge.deviation_old) then
            counter_free = counter_free + 1
        endif      
        if (deviation.ne.deviation.or.counter_free.ge.nfree) then
            call check_deviation()
            alpha(:) = 0.0; U1(:) = 0.0; U2(:) = 0.0; Us(:) = 0.0
            counter_free = 0
            ksi = eta
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !write(*,*) iter, deviation
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! calculating ends-profiles and their characteristics
    do z = 1, n_layer
        n_end1(z)  =  Gf1(z,N1,0)*Gb1(z,N1,0)/WB1(z) * 4.0 + &
        &             Gf1(z,N1,-1)*Gb1(z,N1,1)/WB1(z) + &
        &             Gf1(z,N1,1)*Gb1(z,N1,-1)/WB1(z)
        !!!!
        n_end2(z)  =  Gf2(z,N2,0)*Gb2(z,N2,0)/WB2(z) * 4.0 + &
        &             Gf2(z,N2,-1)*Gb2(z,N2,1)/WB2(z) + &
        &             Gf2(z,N2,1)*Gb2(z,N2,-1)/WB2(z)
    enddo
   	do z = 1, n_layer
    	n_end1(z) = n_end1(z)*dble(N1)/part_fun1
    	n_end2(z) = n_end2(z)*dble(N2)/part_fun2
   	enddo
   	!!!! printing profiles
    if (swpro.ne.0) call print_phi_profiles()
    !!!! first moments
    sum_n_r1 = 0.0;  sum_n_r2 = 0.0
    do z = 1, n_layer
        sum_n_r1 = sum_n_r1 + n_end1(z)*(dble(z)-0.5)
        sum_n_r2 = sum_n_r2 + n_end2(z)*(dble(z)-0.5)
    enddo
    sum1 = 0.0; sum2 = 0.0; m1_1 = 0.0; m1_2 = 0.0
    do z = 1, n_layer
	    sum1 = sum1 + n_end1(z)
	    sum2 = sum2 + n_end2(z)
	    m1_1 = m1_1 + n_end1(z)*(dble(z)-0.5)
	    m1_2 = m1_2 + n_end2(z)*(dble(z)-0.5)
	enddo
    m1_1 = m1_1/sum1
    m1_2 = m1_2/sum2
    !!!! calculating Flory's interaction free energy
    Uint1 = 0.0; Uint2 = 0.0; Uints = 0.0    
    call calc_U_int(chi1, phi_solv, 1.0_8, phi_pol1, Uint1)
    call calc_U_int(chi2, phi_solv, 1.0_8, phi_pol2, Uint2)
    call calc_U_int(chi1, phi_pol1, 0.0_8, phi_solv, Uints)
    call calc_U_int(chi2, phi_pol2, 0.0_8, phi_solv, Uints)       
    call calc_U_int(chi12, phi_pol1, 0.0_8, phi_pol2, Uint2)
    call calc_U_int(chi12, phi_pol2, 0.0_8, phi_pol1, Uint1)
    !!!! calculating total free energy
    F = 0.0
    if (sigma1.gt.0.0) F = F + sigma1*log(dble(N1)*sigma1/part_fun1)
    if (sigma2.gt.0.0) F = F + sigma2*log(dble(N2)*sigma2/part_fun2)
    F = F + Uint1 + Uint2 + Uints  &
    & - sum(alpha(1:n_layer)) &
    & - sum(U1(1:n_layer)*phi_pol1(1:n_layer)) &
    & - sum(U2(1:n_layer)*phi_pol2(1:n_layer)) &
    & - sum(Us(1:n_layer)*phi_solv(1:n_layer)) &
    & - dble(N1)*sigma1*(chi1+chi12-1.0) - dble(N2)*sigma2**(chi2+chi12-1.0) 
    !!!!
    call print_data(m1_1, m1_2, F)
    
    return
end subroutine main
!######################################################################################
