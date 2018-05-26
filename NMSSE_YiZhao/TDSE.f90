! time evolution based on Markovian stochastic Schrodinger equation (Eq. 22) using Chebyshev polynomials methods
subroutine solve_MSSE(n_time, t_total, n_eq, psi0, third_term, psi)

  use constants

  integer, intent(in) :: n_time, n_eq
  double precision, intent(in) :: t_total
  double complex, intent(in) :: psi0(n_eq)
  double complex, intent(in) :: third_term(n_eq, n_eq, n_time)
  double complex, intent(out) :: psi(n_time, n_eq)
  double precision :: time, interval
  double complex, allocatable:: Hamiltonian(:, :)
  double complex, allocatable :: U(:, :)
  double complex :: alpha0, beta0
  
  interval = t_total / n_time
  alpha0 = dcmplx(1.0d0, 0.0d0)
  beta0 = dcmplx(0.0d0, 0.0d0)
  allocate(Hamiltonian(n_eq, n_eq))
  allocate(U(n_eq, n_eq))

  do i = 1, n_time
    time = (i - 1) * interval
    call get_Hamiltonian(time, Hamiltonian)
    Hamiltonian = Hamiltonian + third_term(:, :, i)
    call expansion_Hamiltonian_1(n_eq, Hamiltonian, -interval, U)
    if(i == 1) then
      call zgemm('N', 'N', n_eq, 1, n_eq, alpha0, U, n_eq, psi0, n_eq, beta0, psi(i, :), n_eq)
    else
      call zgemm('N', 'N', n_eq, 1, n_eq, alpha0, U, n_eq, psi(i - 1, :), n_eq, beta0, psi(i, :), n_eq)
    end if
  end do
  
  deallocate(Hamiltonian)
  deallocate(U)
  
end subroutine