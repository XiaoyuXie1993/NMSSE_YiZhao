! The first two terms in eq 20 and 22 (Hel + F(t))
subroutine get_Hamiltonian(t, Hamiltonian)

  use spectral_density

  double precision, intent(in) :: t
  double complex, intent(out) :: Hamiltonian(N_basis, N_basis)
  double complex :: force(N_basis, N_basis)
!  double complex :: H0_complex
  double precision :: omega, real_f, imag_f

  force = 0.0d0

  do i1 = 1, N_basis; do i2 = 1, i1
    do j = 1, N_omega
      omega = j * interval_omega
      force(i1, i2) = force(i1, i2) + dsqrt(C(i1, i2, j)) * dcmplx(real_f(beta, omega, phi(i1, i2, j), t), imag_f(beta, omega, phi(i1, i2, j), t))
    end do
  end do; end do

  Hamiltonian = H0 + force
  
  do i1 = 1, N_basis; do i2 = i1 + 1, N_basis
    Hamiltonian(i1, i2) = dconjg(Hamiltonian(i2, i1))
  end do; end do
  
!  write(22, '(5f10.5)') t, Hamiltonian(1, 1), Hamiltonian(2, 2)
  
end subroutine

!! Debye-Drude spectral density
subroutine discretization()

  use constants
  use spectral_density

  double precision :: omega, time

  SP = 0.0d0

  do i = 1, N_basis
    do j = 1, N_omega
      omega = j * interval_omega
      SP(i, i, j) = eta * omega * omega_c / (omega ** 2 + omega_c ** 2)
!      write(22, '(2f10.5)') omega, SP(i, i, j)      
    end do
!    stop
  end do

  C = 0.0d0
  do i1 = 1, N_basis; do i2 = 1, N_basis
    do j = 1, N_omega
!        omega = j * interval_omega
        C(i1, i2, j) = SP(i1, i2, j) * interval_omega * hbar / pi
    end do
  end do; end do

end subroutine

!subroutine thirdterm_NM()
!
!  use constants
!  use spectral_density
!  use time_evolution
!  
!  double precision :: time, omega
!  double complex :: e(N_basis, N_basis)
!  double complex :: dthirdterm(N_basis, N_basis, time_steps + 1)
!  double complex :: expiHt_plus(N_basis, N_basis, time_steps + 1), expiHt_minus(N_basis, N_basis, time_steps + 1)
!  double complex :: U_plus(N_basis, N_basis), U_minus(N_basis, N_basis)
!  double complex :: tmp1(N_basis, N_basis), tmp2(N_basis, N_basis)
!  double complex :: alpha0, beta0
!  
!  third_term_NM = 0.0d0
!  e = dcmplx(0.0d0, 0.0d0)
!  do i = 1, N_basis
!    e(i, i) = dcmplx(1.0d0, 0.0d0)
!  end do
!  alpha0 = dcmplx(1.0d0, 0.0d0)
!  beta0 = dcmplx(0.0d0, 0.0d0)
!
!  call expansion_Hamiltonian_real(N_basis, H0, interval_time, U_plus)
!  call expansion_Hamiltonian_real(N_basis, H0, -interval_time, U_minus)
!  
!  dthirdterm = 0.0d0
!  do i = 1, time_steps + 1
!    time = (i - 1) * interval_time
!    if(i == 1) then
!      expiHt_plus(:, :, i) = e
!      expiHt_minus(:, :, i) = e
!    else
!      call zgemm('N', 'N', N_basis, N_basis, N_basis, alpha0, U_plus, N_basis, expiHt_plus(:, :, i - 1), N_basis, beta0, expiHt_plus(:, :, i), N_basis)
!      call zgemm('N', 'N', N_basis, N_basis, N_basis, alpha0, U_minus, N_basis, expiHt_minus(:, :, i - 1), N_basis, beta0, expiHt_minus(:, :, i), N_basis)
!    end if
!!    write(*, *) i
!!    write(*, '(4f20.7)') expiHt_plus(1, :, i)
!!    write(*, '(4f20.7)') expiHt_plus(2, :, i)
!!    write(*, *)
!!    if (i == 3) stop
!    do j = 1, N_omega
!      omega = j * interval_omega
!!      write(*, '(2f20.7)') C(1, :, j)
!!      write(*, '(2f20.7)') C(2, :, j)
!!      write(*, *)
!      call dzgemm('N', 'N', N_basis, N_basis, N_basis, alpha0, C(:, :, j), N_basis, expiHt_plus(:, :, i), N_basis, beta0, tmp1, N_basis)
!      call zgemm('N', 'N', N_basis, N_basis, N_basis, alpha0, expiHt_minus(:, :, i), N_basis, tmp1, N_basis, beta0, tmp2, N_basis)
!!      write(*, '(4f20.7)') tmp2(1, :)
!!      write(*, '(4f20.7)') tmp2(2, :)
!!      write(*, *)
!!      stop
!      tmp2 = tmp2 * zexp(dcmplx(0.0d0, -omega * time)) * dcmplx(0.0d0, -1.0d0) / hbar
!      call dzgemm('N', 'N', N_basis, N_basis, N_basis, alpha0, C(:, :, j), N_basis, tmp2, N_basis, alpha0, dthirdterm(:, :, i), N_basis)
!    end do
!    if(i == 1) then
!      third_term_NM(:, :, i) = dcmplx(0.0d0, 0.0d0)
!    else
!      third_term_NM(:, :, i) = third_term_NM(:, :, i - 1) + (dthirdterm(:, :, i - 1) +  dthirdterm(:, :, i)) * interval_time / 2.0d0
!    end if
!    write(22, '(5ES20.10)') time, dthirdterm(1, 1, i), third_term_NM(1, 1, i)
!  end do
!  stop
!  
!end subroutine

subroutine thirdterm_NM()

  use constants
  use spectral_density
  use time_evolution
  
  double precision :: time, omega
  double complex :: e(N_basis, N_basis)
  double complex :: thirdterm
  double complex :: exp_omega(N_omega)
  double complex :: alpha0, beta0
  
  third_term_NM = 0.0d0
  e = dcmplx(0.0d0, 0.0d0)
  do i = 1, N_basis
    e(i, i) = dcmplx(1.0d0, 0.0d0)
  end do
  alpha0 = dcmplx(1.0d0, 0.0d0)
  beta0 = dcmplx(0.0d0, 0.0d0)
  
  do i = 1, time_steps + 1
    time = (i - 1) * interval_time
    do j = 1, N_omega
      omega = j * interval_omega
      exp_omega(j) = (cdexp(dcmplx(0.0d0, -omega * time))  - 1.0d0) / (hbar * omega)
    end do
    thirdterm = 0.0d0
    call dzgemm('T', 'N', 1, 1, N_omega, alpha0, C(1, 1, :), N_omega, exp_omega, N_omega, beta0, thirdterm, 1)
    third_term_NM(:, :, i) = thirdterm * e
!    write(22, '(f10.5, 4ES20.10)') time, third_term_NM(1, 1, i), third_term_NM(2, 2, i)
  end do
!  stop
  
end subroutine

double precision function real_f(beta, omega, phi, t)

  use constants
  
  double precision, intent(in) :: beta, omega, phi, t
  
  real_f = dsqrt(1.0d0 / dtanh(omega * hbar * beta * 0.5d0) + 1.0d0 / dsinh(omega * hbar * beta * 0.5d0)) * dcos(omega * t + phi)

end function

double precision function imag_f(beta, omega, phi, t)

  use constants
  
  double precision, intent(in) :: beta, omega, phi, t
  
  imag_f = dsqrt(1.0d0 / dtanh(omega * hbar * beta * 0.5d0) - 1.0d0 / dsinh(omega * hbar * beta * 0.5d0)) * dsin(omega * t + phi)

end function