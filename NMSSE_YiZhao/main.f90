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
  parameter(N_basis = 2)
! Hamiltonian elements
  double complex :: H0(n_basis, n_basis)
  parameter(H0 = (/0.0d0, 1.0d0, 1.0d0, 0.0d0/))

end module

module spectral_density

  use Hamiltonian_electronic

!! parameters of Debye-Drude spectral density J(omega) = eta * omega * omega_c / (omega ^ 2 + omega_c ^ 2)
  double precision :: eta, omega_c
  double precision :: beta
!  parameter(alpha = 1.2d0)
!  parameter(omega_c = 2.5d0)
!  parameter(beta = 0.2d0)
!  parameter(beta = 5.0d0)
!  parameter(alpha = 0.1d0)
!  parameter(omega_c = 7.5d0)
!! parameters for discretization of spectral density
  integer :: N_omega
  parameter(N_omega = 10000)
  double precision :: interval_omega, omega_max
  parameter(omega_max = 100.0d0)
  double precision :: SP(N_basis, N_basis, N_omega)
  double precision :: phi(N_basis, N_basis, N_omega)
  double precision :: C(N_basis, N_basis, N_omega)

end module

module time_evolution

  use Hamiltonian_electronic

!! parameters of initial electronic states and time-dependent simulation
  double complex :: psi0(N_basis)
  parameter(psi0 = (/1.0d0, 0.0d0/))
  integer :: time_steps
!  parameter(time_steps = 1000)
  double precision :: interval_time, total_time
!  parameter(total_time = 15.0d0)
!! parameters of statistic average
  integer :: N_statistic
!  parameter(N_statistic = 500)
  double complex, allocatable :: U_minus(:, :, :), U_plus(:, :, :)

end module

!! time evolution of a spin-boson system using parallel programing via mpi (J. Chem. Phys. 2013, 138, 014111)
program spin_boson
  
  use time_evolution
  use mpi

  integer :: ierr, num_procs, my_id, n_traj_per_para
  integer :: stat(N_basis)
  double precision :: density
  double complex, allocatable :: psi(:, :, :)
  double precision, allocatable :: diff_density(:)
  
  call initial()
  allocate(diff_density(time_steps))
  call discretization()
  call init_random_seed()

!  write(*, '(i4, 2f14.7)') time_steps, interval_time, total_time

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  n_traj_per_para = (N_statistic - 1) / num_procs + 1
  allocate(psi(n_traj_per_para, time_steps, N_basis))
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  diff_density = 0.0d0
  do i = 1, n_traj_per_para
    call initialphi()
    call Runge_Kutta(time_steps, total_time, N_basis, psi0, psi(i, :, :))
    if(my_id /= 0) then
      call MPI_SEND(psi(i, :, :), time_steps * N_basis, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, stat, ierr)
    else
      do k = 1, time_steps
        diff_density(k) = diff_density(k) + (density(psi(i, k, 1)) - density(psi(i, k, 2))) / (n_traj_per_para * num_procs)
      end do
      do j = 1, num_procs - 1
        call MPI_RECV(psi(i, :, :), time_steps * N_basis, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, stat, ierr)
        do k = 1, time_steps
          diff_density(k) = diff_density(k) + (density(psi(i, k, 1)) - density(psi(i, k, 2))) / (n_traj_per_para * num_procs)
        end do
      end do
    end if
  end do

  call MPI_FINALIZE(ierr)

  if(my_id == 0) then
    time = 0.0d0
    do i = 1, time_steps
      write(22, '(2f14.7)') time, diff_density(i)
      time = time + interval_time
    end do
  end if

end program 

subroutine initial()

  use spectral_density
  use time_evolution
  
  character*12 :: ch

  open(11, file = 'input')
! paramters in spectral_density
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) eta
  read(11, '(A)', advance = 'no') ch
  read(11, *) omega_c
  read(11, '(A)', advance = 'no') ch
  read(11, *) beta
!  write(*, '(3f10.5, l5)') eta, omega_c, beta
!  stop
! paramters in time_evolution
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) total_time
  interval_time = total_time / time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) N_statistic
!  write(*, '(i7, f10.5, i7)') time_steps, total_time, N_statistic
  close(11)
  allocate(U_minus(N_basis, N_basis, N_statistic), U_plus(N_basis, N_basis, N_statistic))

end subroutine
  
subroutine initialphi()

  use constants
  use spectral_density
  
  do i1 = 1, N_basis; i2 = 1, N_basis
    do j = 1, N_omega
      call random_number(phi(i1, i2, j))
      phi(i1, 12, j) = phi(i1, i2, j) * 2.0d0 * pi
    end do
  end do
  
end subroutine

! initial of random seed from fortran manual  
          SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
          END SUBROUTINE

! density (population) of a electronic state  
double precision function density(coeff)
  
  double complex, intent(in) :: coeff

  density = real(coeff) ** 2 + imag(coeff) ** 2

end function