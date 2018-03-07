! time evolution based on Markovian stochastic Schrodinger equation (Eq. 22) using Chebyshev polynomials methods
subroutine solve_MSSE(n_time, t_total, n_eq, psi0, psi)

  use constants

  integer, intent(in) :: n_time, n_eq
  double precision, intent(in) :: t_total
  double complex, intent(in) :: psi0(n_eq)
  double complex, intent(out) :: psi(n_time, n_eq)
  double precision :: time, interval
  double complex :: Hamiltonian(n_eq, n_eq)
  double complex :: third_term_NM(n_eq, n_eq, n_time)
  double complex :: U(n_eq, n_eq)

!  call thirdterm_NM(n_time, t_total, third_term_NM)

  interval = t_total / n_time
  
  do i = 1, n_time
    time = i * interval
    call get_Hamiltonian(time, Hamiltonian)
!    write(*, '(4f20.7)') Hamiltonian(1, :)
!    write(*, '(4f20.7)') Hamiltonian(2, :)
!    write(*, *)
!    Hamiltonian = Hamiltonian + third_term_NM(:, :, i)
!    write(*, '(4f20.7)') Hamiltonian(1, :)
!    write(*, '(4f20.7)') Hamiltonian(2, :)
!    write(*, *)
!    stop
    call expansion_Hamiltonian(n_eq, Hamiltonian, -interval, U)
!    write(*, '(4f20.7)') U(1, :)
!    write(*, '(4f20.7)') U(2, :)
!    write(*, *)
!    if(i == 20) stop
!    write(*, '(4f20.7)') U(1, :)
!    write(*, '(4f20.7)') U(2, :)
!    write(*, *)
    if(i == 1) then
      call zgemm('N', 'N', n_eq, 1, n_eq, 1.0d0, U, n_eq, psi0, n_eq, 0.0d0, psi(i, :), n_eq)
    else
      call zgemm('N', 'N', n_eq, 1, n_eq, 1.0d0, U, n_eq, psi(i - 1, :), n_eq, 0.0d0, psi(i, :), n_eq)
    end if
  end do

end subroutine

subroutine thirdterm_NM(time_steps, total_time, third_term_NM)

  use constants
  use spectral_density
  
  integer, intent(in) :: time_steps
  double precision, intent(in) :: total_time
  double complex, intent(out) :: third_term_NM(N_basis, N_basis, time_steps)
  double complex :: dthirdterm(N_basis, N_basis, time_steps)
  double precision :: time, omega, interval
  double complex :: expiHt_plus(N_basis, N_basis, time_steps), expiHt_minus(N_basis, N_basis, time_steps)
  double complex :: U_plus(N_basis, N_basis), U_minus(N_basis, N_basis)
  double complex :: tmp1(N_basis, N_basis), tmp2(N_basis, N_basis), ddthirdterm(N_basis, N_basis)
  
  interval = total_time / time_steps
  dthirdterm = 0.0d0
  third_term_NM = 0.0d0

  call expansion_Hamiltonian_real(N_basis, H0, interval, U_plus)
  call expansion_Hamiltonian_real(N_basis, H0, -interval, U_minus)

  do i = 1, time_steps
    time = i * interval
    if(i == 1) then
      expiHt_plus(:, :, i) = U_plus
      expiHt_minus(:, :, i) = U_minus
    else
      call zgemm('N', 'N', N_basis, N_basis, N_basis, 1.0d0, U_plus, N_basis, expiHt_plus(:, :, i - 1), N_basis, 0.0d0, expiHt_plus(:, :, i), N_basis)
      call zgemm('N', 'N', N_basis, N_basis, N_basis, 1.0d0, U_minus, N_basis, expiHt_minus(:, :, i - 1), N_basis, 0.0d0, expiHt_minus(:, :, i), N_basis)
    end if
!    write(*, *) i
!    write(*, '(4f20.7)') expiHt_plus(1, :, i)
!    write(*, '(4f20.7)') expiHt_plus(2, :, i)
!    write(*, *)
!    if (i == 3) stop
    do j = 1, N_omega
      omega = j * interval_omega
      call dzgemm('N', 'N', N_basis, N_basis, N_basis, 1.0d0, C(:, :, j), N_basis, expiHt_plus, N_basis, 0.0d0, tmp1, N_basis)
      call zgemm('N', 'N', N_basis, N_basis, N_basis, 1.0d0, expiHt_minus, N_basis, tmp1, N_basis, 0.0d0, tmp2, N_basis)
      tmp2 = tmp2 * cdexp(dcmplx(0.0d0, -omega * time))
      call dzgemm('N', 'N', N_basis, N_basis, N_basis, 1.0d0, C(:, :, j), N_basis, tmp2, N_basis, 0.0d0, ddthirdterm, N_basis)
      ddthirdterm = ddthirdterm * dcmplx(0.0d0, -1.0d0)
      dthirdterm(:, :, i) = dthirdterm(:, :, i) + ddthirdterm * interval
    end do
    if(i == 1) then
      third_term_NM(:, :, i) = dthirdterm(:, :, i)
    else
      third_term_NM(:, :, i) = third_term_NM(:, :, i - 1) + dthirdterm(:, :, i)
    end if
  end do

end subroutine
