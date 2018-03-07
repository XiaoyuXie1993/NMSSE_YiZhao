! The first two terms in eq 20 and 22 (Hel + F(t))
subroutine get_Hamiltonian(t, Hamiltonian)

  use constants
  use spectral_density

  double precision, intent(in) :: t
  double complex, intent(out) :: Hamiltonian(n_basis, n_basis)
  double complex :: force(n_basis, n_basis)
!  double complex :: H0_complex
  double precision :: omega, real_f, imag_f

  force = 0.0d0

  do i1 = 1, N_basis; do i2 = 1, N_basis
    do j = 1, N_omega
      omega = j * interval_omega
      force(i1, i2) = force(i1, i2) + C(i1, i2, j) *  &
                      dcmplx(real_f(beta, omega, phi(i1, i2, j), t), imag_f(beta, omega, phi(i1, i2, j), t))
    end do
  end do; end do

  Hamiltonian = H0 + force

end subroutine

!! Debye-Drude spectral density
subroutine discretization()

  use constants
  use spectral_density

  double precision :: omega, time

  interval_omega = omega_max / N_omega
  SP = 0.0d0

  do i = 1, N_basis
    do j = 1, N_omega
      omega = j * interval_omega
      SP(i, i, j) = eta * omega * omega_c / (omega ** 2 + omega_c ** 2)
    end do
  end do

  C = 0.0d0
  do i1 = 1, N_basis; do i2 = 1, N_basis
    do j = 1, N_omega
!        omega = j * interval_omega
        C(i1, i2, j) = dsqrt(SP(i1, i2, j) * interval_omega * hbar / pi)
    end do
  end do; end do

end subroutine

double precision function real_f(beta, omega, phi, t)

  use constants
  
  double precision, intent(in) :: beta, omega, phi, t
  
  real_f = dsqrt(1.0d0 / dtanh(omega * hbar * beta * 0.5d0) + 1.0d0 / dsinh(omega * hbar * beta * 0.5d0)) * &
           dcos(omega * t + phi)

end function

double precision function imag_f(beta, omega, phi, t)

  use constants
  
  double precision, intent(in) :: beta, omega, phi, t
  
  imag_f = dsqrt(1.0d0 / dtanh(omega * hbar * beta * 0.5d0) - 1.0d0 / dsinh(omega * hbar * beta * 0.5d0)) * &
           dsin(omega * t + phi)

end function