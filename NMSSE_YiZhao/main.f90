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
  double precision :: H0(N_basis, N_basis)
  parameter(H0 = (/0.0d0, 1.0d0, 1.0d0, 0.0d0/))

end module

module spectral_density

  use Hamiltonian_electronic

!! parameters of Debye-Drude spectral density J(omega) = eta * omega * omega_c / (omega ^ 2 + omega_c ^ 2)
  double precision :: eta, omega_c
  double precision :: beta
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
!  call thirdterm_NM()

!  write(*, '(i4, 2f14.7)') time_steps, interval_time, total_time

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  n_traj_per_para = (N_statistic - 1) / num_procs + 1
  allocate(psi(n_traj_per_para, time_steps, N_basis))
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  diff_density = 0.0d0
  do i = 1, n_traj_per_para
    call init_random_seed(my_id)
    call initialphi()
    call solve_MSSE(time_steps, total_time, N_basis, psi0, psi(i, :, :))
    if(my_id /= 0) then
      call MPI_SEND(psi(i, :, :), time_steps * N_basis, MPI_DOUBLE_COMPLEX, 0, i, MPI_COMM_WORLD, stat, ierr)
    else
      do k = 1, time_steps
        diff_density(k) = diff_density(k) + (density(psi(i, k, 1)) + density(psi(i, k, 2))) / (n_traj_per_para * num_procs)
      end do
      do j = 1, num_procs - 1
        call MPI_RECV(psi(i, :, :), time_steps * N_basis, MPI_DOUBLE_COMPLEX, j, i, MPI_COMM_WORLD, stat, ierr)
        do k = 1, time_steps
          diff_density(k) = diff_density(k) + (density(psi(i, k, 1)) + density(psi(i, k, 2))) / (n_traj_per_para * num_procs)
        end do
      end do
    end if
  end do

  call MPI_FINALIZE(ierr)

  if(my_id == 0) then
    open(22, file = 'result.dat')
    write(22, '(2f14.7)') 0.0d0, density(psi0(1)) + density(psi0(2))
    do i = 1, time_steps
      time = interval_time * i
      write(22, '(2f14.7)') time, diff_density(i)
    end do
    close(22)
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

end subroutine
  
subroutine initialphi()

  use constants
  use spectral_density

  do i1 = 1, N_basis; do i2 = 1, N_basis
    do j = 1, N_omega
      call random_number(phi(i1, i2, j))
      phi(i1, i2, j) = phi(i1, i2, j) * 2.0d0 * pi
    end do
  end do; end do

end subroutine

!! initial random seed from fortran manual
          subroutine init_random_seed(pid)
            use iso_fortran_env, only: int64
            
            implicit none
            
            integer, intent(in) :: pid
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8)
            integer(int64) :: t
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
!            open(newunit=un, file="/dev/urandom", access="stream", &
!                 form="unformatted", action="read", status="old", iostat=istat)
!            if (istat == 0) then
!               read(un) seed
!               close(un)
!            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
!            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
        end subroutine init_random_seed


! density (population) of a electronic state  
double precision function density(coeff)

  double complex, intent(in) :: coeff

  density = dreal(coeff) ** 2 + dimag(coeff) ** 2

end function