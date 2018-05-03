! time evolution based on Markovian stochastic Schrodinger equation (Eq. 22) using Chebyshev polynomials methods
subroutine solve_MSSE(n_time, t_total, n_eq, psi0, third_term, psi)

  use constants

  integer, intent(in) :: n_time, n_eq
  double precision, intent(in) :: t_total
  double complex, intent(in) :: psi0(n_eq)
  double complex, intent(in) :: third_term(n_eq, n_eq, n_time)
  double complex, intent(out) :: psi(n_time, n_eq)
  double precision :: time, interval
  double complex :: Hamiltonian(n_eq, n_eq)
  double complex :: U(n_eq, n_eq)
  double complex :: alpha0, beta0
  
  interval = t_total / n_time
  alpha0 = dcmplx(1.0d0, 0.0d0)
  beta0 = dcmplx(0.0d0, 0.0d0)

  do i = 1, n_time
    time = (i - 1) * interval
    call get_Hamiltonian(time, Hamiltonian)
!    write(*, '(4f20.7)') Hamiltonian(1, :)
!    write(*, '(4f20.7)') Hamiltonian(2, :)
!    write(*, *)
!    write(22, '(5f14.7)') time, Hamiltonian(1, 1), Hamiltonian(2, 2) 
!    do j1 = 1, n_eq; do j2 = 1, n_eq
!        write(*, '(2f20.7)', advance = 'no') Hamiltonian(j1, j2)
!      end do
!      write(*, *)
!    end do
!    stop
!    write(*, *)
    Hamiltonian = Hamiltonian + third_term(:, :, i)
!    write(33, '(5f10.5)') time, Hamiltonian(1, 1), Hamiltonian(2, 2)
!    write(*, '(4f20.7)') Hamiltonian(1, :)
!    write(*, '(4f20.7)') Hamiltonian(2, :)
!    write(*, *)
!    write(*, '(4f20.7)') Hamiltonian(1, :)
!    write(*, '(4f20.7)') Hamiltonian(2, :)
!    write(*, *)
!    stop
!    call expansion_Hamiltonian_test(n_eq, Hamiltonian, -interval, U)
!    write(*, '(4f20.10)') U(1, :)
!    write(*, '(4f20.10)') U(2, :)
!    write(*, *)
    call expansion_Hamiltonian_1(n_eq, Hamiltonian, -interval, U)
!    write(*, '(4f20.10)') U(1, :)
!    write(*, '(4f20.10)') U(2, :)
!    write(*, *)
!    stop
!    if(i == 20) stop
!    write(*, '(4f20.7)') U(1, :)
!    write(*, '(4f20.7)') U(2, :)
!    write(*, *)
    if(i == 1) then
      call zgemm('N', 'N', n_eq, 1, n_eq, alpha0, U, n_eq, psi0, n_eq, beta0, psi(i, :), n_eq)
    else
      call zgemm('N', 'N', n_eq, 1, n_eq, alpha0, U, n_eq, psi(i - 1, :), n_eq, beta0, psi(i, :), n_eq)
    end if
  end do
  
end subroutine