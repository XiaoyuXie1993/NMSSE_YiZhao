program plot

!! constants
  double precision :: pi, hbar
  parameter(pi = 3.1415926535897932)
  parameter(hbar = 1.0d0)
!! parameters of Debye-Drude spectral density J(omega) = eta * omega * omega_c / (omega ^ 2 + omega_c ^ 2)
  double precision :: eta, omega_c
  parameter (eta = 0.5d0)
  parameter (omega_c = 5.0d0)
!! discretization of spectral density
  integer :: N_omega
  parameter (N_omega = 100000)
  double precision :: interval_omega
  parameter (interval_omega = 0.1d0)
  double precision :: SP
  double precision :: C(N_omega)
!! parameters of thirdterm
  integer :: N_time
  parameter (N_time = 150)
  double precision :: interval_time
  parameter (interval_time = 0.1d0)
  double complex :: exp_omega(N_omega)
!  double complex :: dthirdterm(N_time + 1)
  double complex :: thirdterm(N_time + 1)
  double precision :: omega, time
  double complex :: alpha, beta

  alpha = 1.0d0
  beta = 0.0d0
  C = 0.0d0
  do i = 1, N_omega
    omega = i * interval_omega
    SP = eta * omega * omega_c / (omega ** 2 + omega_c ** 2)
    C(i) = SP * interval_omega * hbar / pi
!    write(22, '(2f10.5)') omega, SP, C(i)
  end do
!stop

  thirdterm = 0.0d0
  do i = 1, N_time + 1
    time = (i - 1) * interval_time
    do j = 1, N_omega
      omega = j * interval_omega
      exp_omega(j) = (cdexp(dcmplx(0.0d0, -omega * time))  - 1.0d0) / (hbar * omega)
!      dthirdterm(i) = dthirdterm(i) + C(j) ** 2 * zexp(dcmplx(0.0d0, -omega * time)) * dcmplx(0.0d0, -1.0d0) / hbar
    end do
    call dzgemm('T', 'N', 1, 1, N_omega, alpha, C, N_omega, exp_omega, N_omega, beta, thirdterm(i), 1)
!    if(i == 1) then
!      thirdterm(i) = 0.0d0
!    else
!      thirdterm(i) = thirdterm(i - 1) + (dthirdterm(i - 1) + dthirdterm(i)) / 2.0d0 * interval_time
!    end if
    write(22, '(f10.5, 2ES20.10)') time, thirdterm(i)!, thirdterm(i)
  end do

end program