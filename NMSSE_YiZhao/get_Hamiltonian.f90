! The first two terms in eq 20 and 22 (Hel + F(t))
subroutine get_Hamiltonian(t, Hamiltonian)

  use constants
  use spectral_density

  double precision, intent(in) :: t
  double complex, intent(out) :: Hamiltonian(N_basis, N_basis)
  double precision :: omega
  double complex :: force
  
  Hamiltonian = H0
  
  force = 0.0d0
  
!  open(22, file = 'H_t.dat', access = 'append')
!  write(22, '(f10.5)', advance = 'no') t
  do i = 1, N_basis
    force = 0.0d0
    do j = 1, N_omega
      omega = j * interval_omega
      force = force + h(j) * (dsqrt(n_therm(j) * 0.5d0 + 0.5d0) * (phi(i, j, 1) + ci * phi(i, j, 2)) * cdexp(ci * omega * t) + dsqrt(n_therm(j) * 0.5d0) * (phi(i, j, 1) - ci * phi(i, j, 2)) * cdexp(-ci * omega * t))
    end do
    Hamiltonian(i, i) = Hamiltonian(i, i) + force
!    write(22, '(2f14.7)', advance = 'no') Hamiltonian(i, i)
  end do
!  write(22, *)
!  close(22)
  
end subroutine

!! Debye spectral density
subroutine discretization()

  use constants
  use spectral_density
  
  double precision :: omega, SP
  
  open(11, file = 'results_ft.dat')
  read(11, *)
  do i = 1, N_omega
    read(11, *) omega, SP
    omega = i * interval_omega
    if(SP < 0.0d0) SP = 0.0d0
    n_therm(i) = 1.0d0 / (dexp(beta * hbar * omega) - 1.0d0)
    h(i) = dsqrt(2.0d0 * SP * hbar * interval_omega)
  end do
  close(11)

end subroutine

subroutine thirdterm_NM()

  use mpi
  use constants
  use spectral_density
  use time_evolution
  
  double precision :: time, omega
  double complex, allocatable :: dthirdterm(:, :, :)
  double complex, allocatable :: expiHt(:, :), tmp(:, :)
  double complex, allocatable :: U(:, :), e(:, :)
  double complex :: alpha0, beta0
  double complex :: alpha_c
  
  allocate(dthirdterm(N_basis, N_basis, 2))
  allocate(expiHt(N_basis, N_basis), tmp(N_basis, N_basis))
  allocate(U(N_basis, N_basis), e(N_basis, N_basis))
  third_term_NM = 0.0d0
  alpha0 = dcmplx(1.0d0, 0.0d0)
  beta0 = dcmplx(0.0d0, 0.0d0)
  e = dcmplx(0.0d0, 0.0d0)
  do i = 1, N_basis
    e(i, i) = dcmplx(1.0d0, 0.0d0)
  end do
  call expansion_Hamiltonian_real(N_basis, H0, interval_time, U)
  
  dthirdterm = 0.0d0
  do i = 1, time_steps
    time = (i - 1) * interval_time
    if (i == 1) then
      expiHt = e
    else
      call zgemm('N', 'N', N_basis, N_basis, N_basis, alpha0, U, N_basis, tmp, N_basis, beta0, expiHt, N_basis)
    end if
    tmp = expiHt
    alpha_c = 0.0d0
    do j = 1, N_omega
      omega = j * interval_omega
      alpha_c = alpha_c - ci / hbar * dcmplx(dtanh(hbar * omega * beta * 0.25d0) * dcos(omega * time), -dsin(omega * time)) * h(j) ** 2
    end do
    do k1 = 1, N_basis; do k2 = 1, N_basis
        dthirdterm(k1, k2, 2) = alpha_c * expiHt(k1, k2) * dconjg(expiHt(k1, k1))
    end do; end do
    if(i == 1) then
      third_term_NM(:, :, i) = 0.0d0
    else
      third_term_NM(:, :, i) = third_term_NM(:, :, i - 1) + (dthirdterm(:, :, 1) + dthirdterm(:, :, 2)) * 0.5d0 * interval_time
    end if
    dthirdterm(:, :, 1) = dthirdterm(:, :, 2)
    dthirdterm(:, :, 2) = 0.0d0
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
