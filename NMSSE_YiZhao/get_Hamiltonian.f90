! The first two terms in eq 20 and 22 (Hel + F(t))
subroutine get_Hamiltonian(t, Hamiltonian)

  use constants
  use spectral_density

  double precision, intent(in) :: t
  double complex, intent(out) :: Hamiltonian(N_basis, N_basis)
  double precision :: x_operator(2, 2)
  double precision :: omega
  double complex :: force
  
  force = 0.0d0
  x_operator = 0.0d0
  x_operator(1, 1) = 1.0d0
  x_operator(2, 2) = -1.0d0
  
  do i = 1, N_omega
    omega = i * interval_omega
    force = force + C(1, 1, i) * (dsqrt(n_therm(i) * 0.5d0 + 0.5d0) * (phi(1, i, 1) + ci * phi(1, i, 2)) * cdexp(ci * omega * t) + dsqrt(n_therm(i) * 0.5d0) * (phi(1, i, 1) - ci * phi(1, i, 2)) * cdexp(-ci * omega * t))
  end do
  
  Hamiltonian = H0 + force * x_operator
  
end subroutine

!! Ohmic spectral density
subroutine discretization()

  use constants
  use spectral_density
  
  double precision :: x_operator(2, 2)
  double precision :: omega, SP, pC
  
  x_operator = 0.0d0
  x_operator(1, 1) = 1.0d0
  x_operator(2, 2) = -1.0d0
  
  do i = 1, N_omega
    omega = i * interval_omega
    SP = pi * 0.5d0 * alpha * omega * dexp(- omega / omega_c)
    n_therm(i) = 1.0d0 / (dexp(omega * hbar * beta) - 1.0d0)
    pC = dsqrt(SP * interval_omega * hbar / pi)
    C(:, :, i) = pC * x_operator
!    stop
  end do

end subroutine

subroutine thirdterm_NM()

  use mpi
  use constants
  use spectral_density
  use time_evolution
  
  integer :: ierr, num_procs, my_id, n_time_step_per_para, i, i_total, i_other
  integer, allocatable :: stat(:)
  double precision :: time, omega
  double complex, allocatable :: dthirdterm(:, :, :)
  double complex, allocatable :: expiHt_plus(:, :, :), expiHt_minus(:, :, :)
  double complex :: tmp1(N_basis, N_basis), tmp2(N_basis, N_basis)
  double complex :: alpha0, beta0
  
  allocate(stat(N_basis))
  third_term_NM = 0.0d0
  alpha0 = dcmplx(1.0d0, 0.0d0)
  beta0 = dcmplx(0.0d0, 0.0d0)
  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  n_time_step_per_para = (time_steps - 1) / num_procs + 1
  allocate(dthirdterm(N_basis, N_basis, n_time_step_per_para * num_procs))
  allocate(expiHt_plus(N_basis, N_basis, n_time_step_per_para * num_procs), expiHt_minus(N_basis, N_basis, n_time_step_per_para * num_procs))
  
  do i = 1, n_time_step_per_para * num_procs
    time = (i - 1) * interval_time
    call expansion_Hamiltonian_complex(N_basis, H0, time, expiHt_plus(:, :, i))
    call expansion_Hamiltonian_complex(N_basis, H0, -time, expiHt_minus(:, :, i))
  end do
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  
  dthirdterm = 0.0d0
  do i = 1, n_time_step_per_para
    i_total = n_time_step_per_para * my_id + i
    time = (i_total - 1) * interval_time
!    write(*, *) i
!    write(*, '(4f20.7)') expiHt_plus(1, :, i)
!    write(*, '(4f20.7)') expiHt_plus(2, :, i)
!    write(*, *)
!    if (i == 5) stop
    do j = 1, N_omega
      omega = j * interval_omega
!      if(j < 3) then
!      write(*, '(2f20.7)') C(1, :, j)
!      write(*, '(2f20.7)') C(2, :, j)
!      write(*, *)
!      end if
      call dzgemm('N', 'N', N_basis, N_basis, N_basis, alpha0, C(:, :, j), N_basis, expiHt_plus(:, :, i_total), N_basis, beta0, tmp1, N_basis)
      call zgemm('N', 'N', N_basis, N_basis, N_basis, alpha0, expiHt_minus(:, :, i_total), N_basis, tmp1, N_basis, beta0, tmp2, N_basis)
!      if(j < 3) then
!      write(*, '(4f20.7)') tmp2(1, :)
!      write(*, '(4f20.7)') tmp2(2, :)
!      write(*, *)
!      end if
!      stop
      tmp2 = - tmp2 * dcmplx(dtanh(hbar * omega * beta * 0.25d0) * dcos(omega * time), -dsin(omega * time)) * ci / hbar
      call dzgemm('N', 'N', N_basis, N_basis, N_basis, alpha0, C(:, :, j), N_basis, tmp2, N_basis, alpha0, dthirdterm(:, :, i_total), N_basis)
    end do
    do j = 0, num_procs - 1
      if(j /= my_id) then
        call MPI_SEND(dthirdterm(:, :, i_total), N_basis * N_basis, MPI_DOUBLE_COMPLEX, j, i, MPI_COMM_WORLD, ierr)
        i_other = n_time_step_per_para * j + i
        call MPI_RECV(dthirdterm(:, :, i_other), N_basis * N_basis, MPI_DOUBLE_COMPLEX, j, i, MPI_COMM_WORLD, stat, ierr)
      end if
    end do
  end do
  
  do i = 1, time_steps
    if(i == 1) then
      third_term_NM(:, :, i) = 0.0d0
    else
      third_term_NM(:, :, i) = third_term_NM(:, :, i - 1) + (dthirdterm(:, :, i - 1) + dthirdterm(:, :, i)) * 0.5d0 * interval_time
    end if
  end do
  
end subroutine

!subroutine thirdterm_NM()
!
!  use constants
!  use spectral_density
!  use time_evolution
!  
!  double precision :: time, omega
!  double complex :: e(N_basis, N_basis)
!  double complex :: thirdterm
!  double complex :: exp_omega(N_omega)
!  double complex :: alpha0, beta0
!  
!  third_term_NM = 0.0d0
!  e = dcmplx(0.0d0, 0.0d0)
!  do i = 1, N_basis
!    e(i, i) = dcmplx(1.0d0, 0.0d0)
!  end do
!  alpha0 = dcmplx(1.0d0, 0.0d0)
!  beta0 = dcmplx(0.0d0, 0.0d0)
!  
!  do i = 1, time_steps + 1
!    time = (i - 1) * interval_time
!    do j = 1, N_omega
!      omega = j * interval_omega
!      exp_omega(j) = (cdexp(dcmplx(0.0d0, -omega * time))  - 1.0d0) / (hbar * omega)
!    end do
!    thirdterm = 0.0d0
!    call dzgemm('T', 'N', 1, 1, N_omega, alpha0, C(1, 1, :), N_omega, exp_omega, N_omega, beta0, thirdterm, 1)
!    third_term_NM(:, :, i) = thirdterm * e
!!    write(22, '(f10.5, 4ES20.10)') time, third_term_NM(1, 1, i), third_term_NM(2, 2, i)
!  end do
!!  stop
!  
!end subroutine
