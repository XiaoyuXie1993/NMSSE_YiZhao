! time-dependent Schrodinger equation (Markovian stochastic Schrodinger equation from eq 22)
subroutine MSSE(n_eq, dpsi, psi, t)
  
  use constants

  integer, intent(in) :: n_eq
  double complex, intent(in) :: psi(n_eq)
  double precision, intent(in) :: t
  double complex, intent(out) :: dpsi(n_eq)
  double complex :: Hamiltonian(n_eq, n_eq)
  double complex :: third_term(n_eq, n_eq)
  
  call get_Hamiltonian(t, Hamiltonian)
  
  call zgemm('N', 'N', n_eq, n_eq, n_eq, 1.0d0, Hamiltonian, n_eq, y, n_eq, 0.0d0, dy, n_eq)
  
  dy = 1 / hbar * cmplx(0.0d0, -1.0d0) * dy

end subroutine
