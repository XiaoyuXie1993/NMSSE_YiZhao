!! calculate exp(i * H * t / hbar) using Chebyshev polynomials methods for complex H (not Hermitian but Hij = Hji* for i /= j)
!! using exp(i * M * t) = exp(i * (H + i * H_diagonal) * t) ~ exp(- H_diaganol * t / 2)exp(i * H * t)exp(- H_diaganol * t / 2) + O(t ^ 3)
subroutine expansion_Hamiltonian_1(n_matrix, H, t, expiHt)

  use constants

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: t
  double complex, intent(in) :: H(n_matrix, n_matrix)
  double complex, intent(out) :: expiHt(n_matrix, n_matrix)
  double complex, allocatable :: H_Hermitian(:, :)
  double complex, allocatable :: expHermitian(:, :)
  double complex, allocatable :: tmp(:, :)
  double precision, allocatable :: H_diag(:, :)
  double complex, allocatable :: expdiag(:, :)
  double complex :: alpha0, beta0

!  write(*, '(4f14.7)') H(1, :)
!  write(*, '(4f14.7)') H(2, :)
!  write(*, *)
  
  allocate(H_diag(n_matrix, n_matrix), expdiag(n_matrix, n_matrix))
  H_diag = 0.0d0
  expdiag = 0.0d0
  do i = 1, n_matrix
    H_diag(i, i) = dimag(H(i, i))
    expdiag(i, i) = dcmplx(dexp(-H_diag(i, i) * 0.5 * t / hbar), 0.0d0)
  end do
  allocate(H_Hermitian(n_matrix, n_matrix))
  H_Hermitian = H - H_diag * dcmplx(0.0d0, 1.0d0)
  deallocate(H_diag)

!  write(*, '(4f14.7)') H_Hermitian(1, :)
!  write(*, '(4f14.7)') H_Hermitian(2, :)
!  write(*, *)
!  write(*, '(2f14.7)') H_diag(1, :)
!  write(*, '(2f14.7)') H_diag(2, :)
!  write(*, *)
!  write(*, '(4f14.7)') H_commutor(1, :)
!  write(*, '(4f14.7)') H_commutor(2, :)
!  write(*, *)
!stop

  allocate(expHermitian(n_matrix, n_matrix))
  call expansion_Hamiltonian_complex(n_matrix, H_Hermitian, t, expHermitian)
  deallocate(H_Hermitian)

!  write(*, '(4f14.7)') expHermitian(1, :)
!  write(*, '(4f14.7)') expHermitian(2, :)
!  write(*, *)
!  write(*, '(2f14.7)') expdiag(1, :)
!  write(*, '(2f14.7)') expdiag(2, :)
!  write(*, *)
!  write(*, '(4f14.7)') expcommutor(1, :)
!  write(*, '(4f14.7)') expcommutor(2, :)
!  write(*, *)
!  stop

  allocate(tmp(n_matrix, n_matrix))
  alpha0 = dcmplx(1.0d0, 0.0d0)
  beta0 = dcmplx(0.0d0, 0.0d0)
  call zgemm('N', 'N', n_matrix, n_matrix, n_matrix, alpha0, expdiag, n_matrix, expHermitian, n_matrix, beta0, tmp, n_matrix)
  
  deallocate(expHermitian)
  
!  write(*, '(4f14.7)') tmp(1, :)
!  write(*, '(4f14.7)') tmp(2, :)
!  write(*, *)
  
  call zgemm('N', 'N', n_matrix, n_matrix, n_matrix, alpha0, tmp, n_matrix, expdiag, n_matrix, beta0, expiHt, n_matrix)
  
  deallocate(expdiag, tmp)
  
!  write(*, '(4f14.7)') expiHt(1, :)
!  write(*, '(4f14.7)') expiHt(2, :)
!  write(*, *)
!  stop
  
end subroutine

!! calculate exp(i * H * t / hbar) using Chebyshev polynomials methods for complex H (not Hermitian but Hij = Hji* for i /= j)
!! using exp(i * M * t) = exp(i * (H + i * H_diagonal) * t) ~ exp(i * H * t / 2)exp(- H_diaganol * t)exp(i * H * t / 2) + O(t ^ 3)
subroutine expansion_Hamiltonian_2(n_matrix, H, t, expiHt)

  use constants

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: t
  double complex, intent(in) :: H(n_matrix, n_matrix)
  double complex, intent(out) :: expiHt(n_matrix, n_matrix)
  double complex, allocatable :: H_Hermitian(:, :)
  double complex, allocatable :: expHermitian(:, :)
  double complex, allocatable :: tmp(:, :)
  double precision, allocatable :: H_diag(:, :)
  double complex, allocatable :: expdiag(:, :)
  double complex :: alpha0, beta0

!  write(*, '(4f14.7)') H(1, :)
!  write(*, '(4f14.7)') H(2, :)
!  write(*, *)
  
  allocate(H_diag(n_matrix, n_matrix), expdiag(n_matrix, n_matrix))
  H_diag = 0.0d0
  expdiag = 0.0d0
  do i = 1, n_matrix
    H_diag(i, i) = dimag(H(i, i))
    expdiag(i, i) = dcmplx(dexp(-H_diag(i, i) * t / hbar), 0.0d0)
  end do
  allocate(H_Hermitian(n_matrix, n_matrix))
  H_Hermitian = H - H_diag * ci
  deallocate(H_diag)

!  write(*, '(4f14.7)') H_Hermitian(1, :)
!  write(*, '(4f14.7)') H_Hermitian(2, :)
!  write(*, *)
!  write(*, '(2f14.7)') H_diag(1, :)
!  write(*, '(2f14.7)') H_diag(2, :)
!  write(*, *)
!  write(*, '(4f14.7)') H_commutor(1, :)
!  write(*, '(4f14.7)') H_commutor(2, :)
!  write(*, *)
!stop

  allocate(expHermitian(n_matrix, n_matrix))
  call expansion_Hamiltonian_complex(n_matrix, H_Hermitian, 0.5d0 *t, expHermitian)
  deallocate(H_Hermitian)

!  write(*, '(4f14.7)') expHermitian(1, :)
!  write(*, '(4f14.7)') expHermitian(2, :)
!  write(*, *)
!  write(*, '(2f14.7)') expdiag(1, :)
!  write(*, '(2f14.7)') expdiag(2, :)
!  write(*, *)
!  write(*, '(4f14.7)') expcommutor(1, :)
!  write(*, '(4f14.7)') expcommutor(2, :)
!  write(*, *)
!  stop

  allocate(tmp(n_matrix, n_matrix))
  alpha0 = dcmplx(1.0d0, 0.0d0)
  beta0 = dcmplx(0.0d0, 0.0d0)
  call zgemm('N', 'N', n_matrix, n_matrix, n_matrix, alpha0, expdiag, n_matrix, expHermitian, n_matrix, beta0, tmp, n_matrix)
  
  deallocate(expdiag)
  
!  write(*, '(4f14.7)') tmp(1, :)
!  write(*, '(4f14.7)') tmp(2, :)
!  write(*, *)
  
  call zgemm('N', 'N', n_matrix, n_matrix, n_matrix, alpha0, expHermitian, n_matrix, tmp, n_matrix, beta0, expiHt, n_matrix)
  
  deallocate(expHermitian, tmp)
  
!  write(*, '(4f14.7)') expiHt(1, :)
!  write(*, '(4f14.7)') expiHt(2, :)
!  write(*, *)
!  stop
  
end subroutine

!! calculate exp(i * H * t / hbar) using Chebyshev polynomials methods for real Hamiltonian
subroutine expansion_Hamiltonian_real(n_matrix, H, t, expiHt)

  use constants

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: t
  double precision, intent(in) :: H(n_matrix, n_matrix)
  double complex, intent(out) :: expiHt(n_matrix, n_matrix)
  double precision :: Ha, Hb
  double precision :: tH(n_matrix, n_matrix)
  double precision :: jt

  call diagonal_real(n_matrix, H, Ha, Hb)
  jt = Hb * t / hbar
  tH = H
  do i = 1, n_matrix
    tH(i, i) = tH(i, i) - Ha
  end do
  tH = tH / Hb
  call expansion_real(n_matrix, tH, jt, expiHt)
  expiHt = expiHt * cdexp(dcmplx(0.0d0, 1.0d0) * t * Ha / hbar)

end subroutine

!! calculate exp(i * H * t / hbar) using Chebyshev polynomials methods for complex Hamiltonian
subroutine expansion_Hamiltonian_complex(n_matrix, H, t, expiHt)

  use constants

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: t
  double complex, intent(in) :: H(n_matrix, n_matrix)
  double complex, intent(out) :: expiHt(n_matrix, n_matrix)
  double precision :: Ha, Hb
  double complex :: tH(n_matrix, n_matrix)
  double precision :: jt

  call diagonal_complex(n_matrix, H, Ha, Hb)
!  write(*, '(2f14.7)') Ha, Hb
!  write(*, *)
  jt = Hb * t / hbar
  tH = H
  do i = 1, n_Matrix
    tH(i, i) = tH(i, i) - Ha
  end do
  tH = tH / Hb
!  write(*, '(4f14.7)') tH(1, :)
!  write(*, '(4f14.7)') tH(2, :)
!  write(*, *)
  call expansion_complex(n_Matrix, tH, jt, expiHt)
  expiHt = expiHt * cdexp(dcmplx(0.0d0, 1.0d0) * t * Ha / hbar)

end subroutine 