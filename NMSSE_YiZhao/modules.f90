module constants

!! constants
  double precision :: pi, hbar
  parameter(pi = 3.1415926535897932)
  parameter(hbar = 1.0d0)

end module

module Hamiltonian_electronic

!! parameters of electronic Hamiltonian (two level system with spin-boson model)
! basis set of electronic states
  integer :: N_basis
! Hamiltonian elements
  double complex, allocatable :: H0(:, :)

end module

module spectral_density

  use Hamiltonian_electronic

!! parameters of Ohmic spectral density J(omega) = pi/2 * alpha * exp(-omega / omega_c)
  double precision :: alpha, omega_c
  double precision :: beta
!! parameters for discretization of spectral density
  integer :: N_omega
  double precision :: interval_omega
  double precision, allocatable :: SP(:, :, :)
  double precision, allocatable :: C(:, :, :)
!  double precision, allocatable :: S(:, :, :)
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