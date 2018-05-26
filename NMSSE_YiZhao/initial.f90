subroutine initial()

  use constants
  use Hamiltonian_electronic
  use spectral_density
  use time_evolution
  
  integer :: Hamiltonian_index, center_index1, center_index2
  double precision :: temperature
  character*20 :: ch

  open(11, file = 'input')
! Hamiltonian
  read(11, *)

! paramters in molecule
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) Na
  read(11, '(A)', advance = 'no') ch
  read(11, *) Nb
!  write(*, '(2i5)') Na, Nb
!  stop
! paramters in electronic
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) epsilon
  read(11, '(A)', advance = 'no') ch
  read(11, *) Va
  read(11, '(A)', advance = 'no') ch
  read(11, *) Vb
  read(11, '(A)', advance = 'no') ch
  read(11, *) Vab
  N_basis = Na * Nb * 2
  allocate(H0(N_basis, N_basis))
! paramters in spectral_density
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) N_omega
  read(11, '(A)', advance = 'no') ch
  read(11, *) interval_omega
  read(11, '(A)', advance = 'no') ch
  read(11, *) temperature
  beta = 1.0d0 / (kB * temperature)
  allocate(n_therm(N_omega), h(N_omega))
  allocate(phi(N_basis, N_omega, 2))
! paramters in time_evolution
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) total_time
  interval_time = total_time / time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) N_statistic
  allocate(psi0(N_basis))
  center_index1 = Hamiltonian_index(Na / 2 + 1, Nb / 2 + 1, 1)
  center_index2 = Hamiltonian_index((Na + 1) / 2, (Nb + 1) / 2, 2)
  psi0 = 0.0d0
!  psi0(center_index1) = dcmplx(1.0d0, 0.0d0)
!  psi0(center_index2) = dcmplx(1.0d0, 0.0d0)
  psi0(center_index1) = dcmplx(dsqrt(0.5d0), 0.0d0)
  psi0(center_index2) = dcmplx(dsqrt(0.5d0), 0.0d0)
!  write(*, *) Hamiltonian_index(Na / 2, Nb / 2, 1)
!  stop
  allocate(third_term_NM(N_basis, N_basis, time_steps))
!  write(*, '(i7, 2f10.5, i7)') time_steps, total_time, interval_time, N_statistic
  close(11)
  
end subroutine

! model Hamiltonian based on molecular calculation
subroutine initialHamiltonian()
  
  use Hamiltonian_electronic
  
  integer :: a1, a2, b1, b2, k1, k2
  integer :: index1, index2, Hamiltonian_index

  H0 = 0.0d0
  do a1 = 1, Na; do b1 = 1, Nb; do k1 = 1, 2
    index1 = Hamiltonian_index(a1, b1, k1)
    do a2 = 1, Na; do b2 = 1, Nb; do k2 = 1, 2
      index2 = Hamiltonian_index(a2, b2, k2)
      if(a1 == a2 .and. b1 == b2 .and. k1 == k2) then
        H0(index1, index2) = epsilon * 0.001d0
      else if(abs(a1 - a2) == 1 .and. b1 == b2 .and. k1 == k2) then
        H0(index1, index2) = Va * 0.001d0
      else if(a1 == a2 .and. abs(b1 - b2) == 1 .and. k1 == k2) then
        H0(index1, index2) = Vb * 0.001d0
      else if(k1 == 1 .and. k2 == 2) then
        if((a1 == a2 .or. a1 - a2 == 1) .and. (b1 == b2 .or. b1 - b2 == 1)) then
          H0(index1, index2) = Vab * 0.001d0
        end if
      else if(k1 == 2 .and. k2 == 1) then
        if((a1 == a2 .or. a2 - a1 == 1) .and. (b1 == b2 .or. b2 - b1 == 1)) then
          H0(index1, index2) = Vab * 0.001d0
        end if
      end if
    end do; end do; end do
  end do; end do; end do
  
!  open(22, file = 'H.dat')
!  do i = 1, N_basis; do j = 1, N_basis
!    write(22, '(f10.5)', advance = 'no') H0(i, j)
!    end do
!    write(22, *)
!  end do
!  close(22)
  
end subroutine

subroutine initialphi()

  use constants
  use spectral_density
  
  double precision, allocatable :: U1(:), U2(:), tmp1(:), tmp2(:)
  
  allocate(U1(N_omega), U2(N_omega), tmp1(N_omega), tmp2(N_omega))
  
  do i = 1, N_basis
    call random_number(U1)
    call random_number(U2)
    tmp1 = dsqrt(-2.0d0 * dlog(U1))
    tmp2 = 2.0d0 * pi * U2
    phi(i, :, 1) = tmp1 * dcos(tmp2)
    phi(i, :, 2) = tmp1 * dsin(tmp2)
  end do

  deallocate(U1, U2, tmp1, tmp2)
  
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
   open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
      read(un) seed
      close(un)
   else
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
   end if
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

integer function Hamiltonian_index(ia, ib, ik)

  use Hamiltonian_electronic
  
  integer, intent(in) :: ia, ib, ik
  
  Hamiltonian_index = ik + (ia - 1) * 2 + (ib - 1) * 2 * Na

end function
