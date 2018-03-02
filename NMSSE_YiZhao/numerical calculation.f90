subroutine Runge_Kutta(n_time, t_total, n_eq, y0, y)

  integer, intent(in) :: n_time, n_eq
  double precision, intent(in) :: t_total
  double complex, intent(in) :: y0(n_eq)
  double complex, intent(out) :: y(n_time, n_eq)
  double precision :: time, interval, density

  interval = t_total / n_time

  time = 0.0d0
  do i = 1, n_time
    if(i == 1) then
      y(i, :) = y0
    else
      call solve(n_eq, y(i - 1, :), time, interval, y(i, :))
    end if
    time = time + interval
!    write(*, *) i, density(y(i, 1)) + density(y(i, 2))
  end do
!stop

end subroutine

!!using fourth-order Runge-Kutta methods (for MSSE)
subroutine solve(n_eq, y, t, interval, newy)

  integer, intent(in) :: n_eq
  double complex, intent(in) :: y(n_eq)
  double precision, intent(in) :: t, interval
  double complex, intent(out) :: newy(n_eq)
  double complex :: K1(n_eq), K2(n_eq), K3(n_eq), K4(n_eq)

  call MTDSE(n_eq, K1, y, t)
  call MTDSE(n_eq, K2, y + 0.5d0 * interval * K1, t + 0.5d0 * interval)
  call MTDSE(n_eq, K3, y + 0.5d0 * interval * K2, t + 0.5d0 * interval)
  call MTDSE(n_eq, K4, y + interval * K3, t + interval)

  newy = y + interval / 6.0d0 * (K1 + 2.0d0 * K2 + 2.0d0 * K3 + K4)

end subroutine

subroutine integral(n_eq, t0, t, dt, y, inty)

  integer, intent(in) :: n_eq, t0, t
  double precision, intent(in) :: dt
  double complex, intent(in) :: y(n_eq, n_eq, t - t0 + 1)
  double complex, intent(out) :: inty(n_eq, n_eq)
  
  inty = 0.0d0

  do i = 1, t - t0 + 1
    do j1 = 1, n_eq; do j2 = 1, n_eq
    inty(j1, j2) = inty(j1, j2) + dt * y(j1, j2, i)  
  end do

end subroutine