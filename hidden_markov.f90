! A(i,j) == Pr( hidden i -> hidden j)
! B(i,j) == Pr( observed j | hidden i)

subroutine alphapass(alpha,sigma,A,B,pi,observed,N,T,S)
      integer :: N,T,S
      real, dimension(N,T) :: alpha
      real, dimension(N,N) :: A
      real, dimension(N,S) :: B
      real, dimension(N) :: pi
      integer, dimension(T) :: observed
      real, dimension(T) :: sigma
      integer :: i
!f2py intent(out) :: alpha, sigma
!f2py intent(hide), depend(A)  :: N = size(A,1)
!f2py intent(hide), depend(observed) :: T = size(observed)
!f2py intent(hide), depend(B) :: S = size(B,2)

! initial iteration
      alpha(:,1) = pi
      alpha(:,1) = alpha(:,1)*B(:,observed(1))
      sigma(1) = sum(alpha(:,1))
      alpha(:,1) = alpha(:,1)/sigma(1)

! now do the remainder
      do i = 2,T
         alpha(:,i) = matmul(alpha(:,i-1),A)
         alpha(:,i) = alpha(:,i)*B(:,observed(i))
         sigma(i) = sum(alpha(:,i))
         alpha(:,i) = alpha(:,i)/sigma(i)
      end do
      end subroutine

! A(i,j) == Pr( hidden i -> hidden j)
! B(i,j) == Pr( observed j | hidden i)
subroutine betapass(beta,sigma,A,B,observed,N,T,S)
      integer :: N,T,S
      real, dimension(N,T) :: beta
      real, dimension(N,N) :: A
      real, dimension(N,S) :: B
      integer, dimension(T) :: observed
      real, dimension(T) :: sigma
      integer :: i
      real, dimension(N,N) :: At

!f2py intent(out) :: beta
!f2py intent(hide), depend(A)  :: N = size(A,1)
!f2py intent(hide), depend(observed) :: T = size(observed)
!f2py intent(hide), depend(B) :: S = size(B,2)

      At = transpose(A)

! initial iteration
      beta(:,T) = 1.

! now do the remainder
      do i = T,2,-1
         beta(:,i-1) = matmul(beta(:,i)*B(:,observed(i)),At)/sigma(i)
      end do
      end subroutine

subroutine gammapass(alpha,beta,gamma,N,T)
      integer :: N,T
      real, dimension(N,T) :: alpha,beta,gamma
!      integer :: i
!      real totalsum
!f2py intent(out) :: gamma
!f2py intent(hide), depend(alpha)  :: N = size(alpha,1)
!f2py intent(hide), depend(alpha)  :: T = size(alpha,2)

      gamma = alpha*beta
!      do i = 1, T
!         totalsum= sum(gamma(:,i))
!         gamma(:,i) = gamma(:,i)/totalsum
!      end do
      end subroutine

! A(i,j) == Pr( hidden i -> hidden j)
! B(i,j) == Pr( observed j | hidden i)
subroutine generate_data(observed,hidden,pi,A,B,N,S,T)
      integer :: N, S, T
      real, dimension(N):: pi
      real, dimension(N,N):: A
      real, dimension(N,S):: B
      integer, dimension(T) :: observed,hidden
      integer :: i

!f2py intent(out) :: observed, hidden
!f2py intent(hide), depend(pi) :: N = size(pi)
!f2py intent(hide), depend(B)  :: S = size(B,1)
      interface
         subroutine grn(X,R)
         real, dimension(:) :: X
         integer :: R
         end subroutine grn
      end interface
      call random_seed

      ! initial hidden state
      call grn(pi,hidden(1))

      ! initial observed state
      call grn(B(hidden(1),:),observed(1))

      ! now do the rest of the sequence
      do i = 1, T
         call grn(A(hidden(i-1),:),hidden(i))
         call grn(B(hidden(i),:),observed(i))
      end do
      end subroutine

subroutine grn(X, R)
! generates a random number with value 1...len(X)
      real, dimension(:) :: X
      real, dimension(:), allocatable :: cdf
      integer :: N, R, i
      real harvest
!f2py intent(out) :: R

      N = size(X)
      allocate(cdf(N))

      cdf(1) = X(1)
      do i = 2, N
         cdf(i) = cdf(i-1) + X(i)
      end do

! Now generate a number and determine it's location.
      call random_number(harvest)
      do i = 1, N
         if (harvest <= cdf(i)) then
            R = i
            exit
         end if
      end do
      deallocate(cdf)
      return
      end subroutine

