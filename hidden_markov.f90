subroutine alphapass(alpha,sigma,A,B,pi,observed,N,T,S)
      integer :: N,T,S
      real, dimension(0:N-1,0:T-1) :: alpha
      real, dimension(0:N-1,0:N-1) :: A
      real, dimension(0:S-1,0:N-1) :: B
      real, dimension(0:N-1) :: pi
      integer, dimension(0:T-1) :: observed
      real, dimension(0:T-1) :: sigma
      integer :: i
!f2py intent(out) :: alpha, sigma
!f2py intent(hide), depend(A)  :: N = size(A,1)
!f2py intent(hide), depend(observed) :: T = size(observed)
!f2py intent(hide), depend(B) :: S = size(B,1)

! initial iteration
      alpha(:,0) = pi
      alpha(:,0) = alpha(:,0)*B(observed(0),:)
      sigma(0) = sum(alpha(:,0))
      alpha(:,0) = alpha(:,0)/sigma(0)

! now do the remainder
      do i = 1,T-1
         alpha(:,i) = matmul(A,alpha(:,i-1))
         alpha(:,i) = alpha(:,i)*B(observed(i),:)
         sigma(i) = sum(alpha(:,i))
         alpha(:,i) = alpha(:,i)/sigma(i)
      end do

      end subroutine
subroutine betapass(beta,sigma,A,B,observed,N,T,S)
      integer :: N,T,S
      real, dimension(0:N-1,0:T-1) :: beta
      real, dimension(0:N-1,0:N-1) :: A
      real, dimension(0:S-1,0:N-1) :: B
      integer, dimension(0:T-1) :: observed
      real, dimension(0:T-1) :: sigma
      integer :: i
      real, dimension(0:N-1,0:N-1) :: At

!f2py intent(out) :: beta
!f2py intent(hide), depend(A)  :: N = size(A,1)
!f2py intent(hide), depend(observed) :: T = size(observed)
!f2py intent(hide), depend(B) :: S = size(B,1)

      At = transpose(A)

! initial iteration
      beta(:,T-1) = 1.

! now do the remainder
      do i = T-1,1,-1
         beta(:,i-1) = matmul(At,beta(:,i)*B(observed(i),:))/sigma(i)
      end do
      end subroutine

subroutine gammapass(alpha,beta,gamma,N,T)
      integer :: N,T
      real, dimension(0:N-1,0:T-1) :: alpha,beta,gamma
      integer :: i
      real totalsum
!f2py intent(out) :: gamma
!f2py intent(hide), depend(alpha)  :: N = size(alpha,1)
!f2py intent(hide), depend(alpha)  :: T = size(alpha,2)

      gamma = alpha*beta
      do i = 0, T-1
         totalsum= sum(gamma(:,i))
         gamma(:,i) = gamma(:,i)/totalsum
      end do
      end subroutine

subroutine generate_data(observed,hidden,pi,A,B,N,S,T)
      integer :: N, S, T
      real, dimension(0:N-1):: pi
      real, dimension(0:N-1,0:N-1):: A
      real, dimension(0:S-1,0:N-1):: B
      integer, dimension(0:T-1) :: observed,hidden
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
      call grn(pi,hidden(0))

      ! initial observed state
      call grn(B(:,hidden(0)),observed(0))

      ! now do the rest of the sequence
      do i = 1, T-1
         call grn(A(:,hidden(i-1)),hidden(i))
         call grn(B(:,hidden(i)),observed(i))
      end do
      end subroutine

subroutine grn(X, R)
! generates a random number with value 0...len(X)-1
      real, dimension(:) :: X
      real, dimension(:), allocatable :: cdf
      integer :: N, R, i
      real harvest
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
            R = i-1
            exit
         end if
      end do
      deallocate(cdf)
      return
      end subroutine

