module constants

!! constants
  double precision :: pi, hbar, kB
  double complex :: ci
  parameter(pi = 3.1415926535897932)
! reduced Planck constant with eV * fs
  parameter (hbar = 0.658211951440)
! Boltzmann constant with eV * K-1
  parameter (kB = 0.000086173304)
  parameter (ci = dcmplx(0.0d0, 1.0d0))
  
end module

module Hamiltonian_electronic

! number of molecules on a/b direction (2D system)
  integer :: Na, Nb
! parameters for model Hamiltonian
  double precision :: epsilon, Va, Vb, Vab
! basis set of electronic states
  integer :: N_basis
! Hamiltonian elements
  double precision, allocatable :: H0(:, :)

end module

module spectral_density

  use Hamiltonian_electronic

  double precision :: beta
!! parameters for discretization of spectral density
  integer :: N_omega
  double precision :: interval_omega
  double precision, allocatable :: n_therm(:), h(:)
  double precision, allocatable :: C(:, :, :)
  double precision, allocatable :: phi(:, :, :)

end module

module time_evolution

  use Hamiltonian_electronic

!! parameters of initial electronic states and time-dependent simulation
  double complex, allocatable :: psi0(:)
! NM term
  double complex, allocatable :: third_term_NM(:, :, :)
  integer :: time_steps
  double precision :: interval_time, total_time
!! parameters of statistic average
  integer :: N_statistic

end module