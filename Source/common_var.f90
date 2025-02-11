module common_var		
	! integer common variables		
	integer(4) :: N1           ! - polymerization degree (chain 1)
	integer(4) :: N2           ! - polymerization degree (chain 2)
	integer(4) :: counter_free ! - counter for "free steps" (gradient descent)
	integer(4) :: ioer         ! - input/output error label
	integer(4) :: iter         ! - number of algorithm iterations
	integer(4) :: n_layer      ! - number of layer in z-direction 
	integer(4) :: nfree        ! - number of "free steps" at gradient descent
	integer(4) :: swpro        ! - if swpro=0: switch off print of profiles  
	integer(4) :: zp1 	   ! - z-coordinate of the grafting point (A-chain)
	integer(4) :: zp2 	   ! - z-coordinate of the grafting point (B-chain)
	integer(4) :: left_w       ! - type of left wall  (0 - solid, 1 - mirror)
	integer(4) :: right_w 	   ! - type of right wall (0 - solid, 1 - mirror)
	! real common variables	
	real(8) :: chi1          ! - Flory's parameter for polymer-solvent (chain 1)
	real(8) :: chi2          ! - Flory's parameter for polymer-solvent (chain 2) 
	real(8) :: chi12         ! - Flory's parameter for polymer-polymer
	real(8) :: deviation     ! - deviations from the incompressibility conditions
	real(8) :: deviation_old ! - old value (at iter-1) of deviation
	real(8) :: eta           ! - step size of gradient descent
	real(8) :: ksi           ! - step size of gradient descent for interaction potential
	real(8) :: lambda_1b     ! - transition probability at backfold fracture (chain 1)
	real(8) :: lambda_1p     ! - transition probability at perpendicular kink (chain 1)
	real(8) :: lambda_1s     ! - transition probability at straight transition (chain 1)
	real(8) :: lambda_2b     ! - transition probability at backfold fracture (chain 2)
	real(8) :: lambda_2p     ! - transition probability at perpendicular kink (chain 2)
	real(8) :: lambda_2s     ! - transition probability at straight transition (chain 2)
	real(8) :: lead_time     ! - algorithm lead time
	real(8) :: p1            ! - Kuhn length (chain 1)
	real(8) :: p2            ! - Kuhn length (chain 2)
	real(8) :: part_fun1     ! - partition function (chain 1)
	real(8) :: part_fun2     ! - partition function (chain 2)
	real(8) :: sigma1        ! - grafting density (chain 1)
	real(8) :: sigma2        ! - grafting density (chain 2)
	! one-dimension allocatable arrays	
	real(8), allocatable, dimension(:) :: alpha    ! - Lagrange field
	real(8), allocatable, dimension(:) :: phi_pol1 ! - volume fraction profile (chain 1)
	real(8), allocatable, dimension(:) :: phi_pol2 ! - volume fraction profile (chain 2)
	real(8), allocatable, dimension(:) :: phi_solv ! - volume fraction profile (solvent)
	real(8), allocatable, dimension(:) :: WB1      ! - Boltzmann's weight (chain 1)
	real(8), allocatable, dimension(:) :: WB2      ! - Boltzmann's weight (chain 2)
	real(8), allocatable, dimension(:) :: WBs      ! - Boltzmann's weight (solvent)
	real(8), allocatable, dimension(:) :: n_end1   ! - ends distribution (chain 1)
	real(8), allocatable, dimension(:) :: n_end2   ! - ends distribution (chain 2)
	real(8), allocatable, dimension(:) :: U1       ! - interaction potential (chain 1)
	real(8), allocatable, dimension(:) :: U2       ! - interaction potential (chain 2)
	real(8), allocatable, dimension(:) :: Us       ! - interaction potential (solvent)
	real(8), allocatable, dimension(:) :: W1       ! - old value (at iter-1) of U1
	real(8), allocatable, dimension(:) :: W2       ! - old value (at iter-1) of U2
	real(8), allocatable, dimension(:) :: Ws       ! - old value (at iter-1) of Us
	! two-dimension allocatable arrays	
	real(8), allocatable, dimension(:,:) :: phi1   ! - phi_pol1 separated to orientation
	real(8), allocatable, dimension(:,:) :: phi2   ! - phi_pol2 separated to orientation
	! three-dimension allocatable arrays	
	real(8), allocatable, dimension(:,:,:) :: Gf1  ! - forward propagator (chain 1)
	real(8), allocatable, dimension(:,:,:) :: Gf2  ! - forward propagator (chain 2)
	real(8), allocatable, dimension(:,:,:) :: Gb1  ! - back propagator (chain 1)
	real(8), allocatable, dimension(:,:,:) :: Gb2  ! - back propagator (chain 2)
    
Contains

!******************************************************************************
subroutine Allocate_all_variables()
!******************************************************************************

	allocate(alpha(0:n_layer+1))
    allocate(phi_solv(0:n_layer+1))
    allocate(phi_pol1(0:n_layer+1))
    allocate(phi_pol2(0:n_layer+1))
    allocate(WB1(0:n_layer+1))
    allocate(WB2(0:n_layer+1))
    allocate(WBs(0:n_layer+1))
    allocate(n_end1(1:n_layer))
    allocate(n_end2(1:n_layer))
    allocate(phi1(0:n_layer+1,-1:1))
    allocate(phi2(0:n_layer+1,-1:1))
    allocate(U1(0:n_layer+1))
    allocate(U2(0:n_layer+1))
    allocate(Us(0:n_layer+1))
    allocate(W1(0:n_layer+1))
    allocate(W2(0:n_layer+1))
    allocate(Ws(0:n_layer+1))
    allocate(Gf1(0:n_layer+1,0:N1,-1:1))
    allocate(Gb1(0:n_layer+1,1:N1,-1:1))
    allocate(Gf2(0:n_layer+1,0:N2,-1:1))
    allocate(Gb2(0:n_layer+1,1:N2,-1:1))
    !!!! zero-initialization of arrays
    alpha(:) = 0.0; WB1(:) = 0.0; WB2(:) = 0.0; WBs(:) = 0.0; phi_solv(:) = 0.0
    phi1(:,:) = 0.0; phi_pol1(:) = 0.0;  Gf1(:,:,:) = 0.0; Gb1(:,:,:) = 0.0
    n_end1(:) = 0.0; n_end2(:) = 0.0; phi2(:,:) = 0.0; phi_pol2(:) = 0.0
    Gf2(:,:,:) = 0.0; Gb2(:,:,:) = 0.0
    U1(:) = 0.0; U2(:) = 0.0; Us(:) = 0.0
    W1(:) = 0.0; W2(:) = 0.0; Ws(:) = 0.0
	
	return
end subroutine
!******************************************************************************

!******************************************************************************
subroutine Deallocate_all_variables()
!******************************************************************************
	
	if (allocated(alpha))    deallocate(alpha)
	if (allocated(phi_solv)) deallocate(phi_solv)
	if (allocated(phi_pol1)) deallocate(phi_pol1)
	if (allocated(phi_pol2)) deallocate(phi_pol2)
	if (allocated(WB1))      deallocate(WB1)
	if (allocated(WB2))      deallocate(WB2)
	if (allocated(WBs))      deallocate(WBs)
	if (allocated(n_end1))   deallocate(n_end1)
	if (allocated(n_end2))   deallocate(n_end2)
	if (allocated(phi1))     deallocate(phi1)
	if (allocated(phi2))     deallocate(phi2)
	if (allocated(U1))       deallocate(U1)
	if (allocated(U2))       deallocate(U2)
	if (allocated(Us))       deallocate(Us)
	if (allocated(W1))       deallocate(W1)
	if (allocated(W2))       deallocate(W2)
	if (allocated(Ws))       deallocate(Ws)
	if (allocated(Gf1))      deallocate(Gf1)
	if (allocated(Gf2))      deallocate(Gf2)
	if (allocated(Gb1))      deallocate(Gb1)
	if (allocated(Gb2))      deallocate(Gb2)
	
	return
end subroutine
!******************************************************************************

end module common_var
