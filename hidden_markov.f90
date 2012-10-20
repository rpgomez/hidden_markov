! alpha(i,t) = Pr(hidden_t = i | obs_1,...,obs_t)
! beta(i,t)  = Pr(obs_t+1,...,obs_T| hidden_t = i)*Pr(o_1,...,o_t)/Pr(o_1,...,o_T)
! gamma(i,t)  = Pr(hidden_t = i | obs_1,...,obs_T)
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

! Reestimates Pr(i -> j)
subroutine digammapass(alpha,beta,observed,digamma,A,B,N,T,S)
! Pr(x_t=i,x_t+1=j|o_1,...,o_T) = alpha_t^(i)A_ijbeta_^t+1(j)b(j,o_t+1)

      integer :: N,T,S
      real, dimension(N,T) :: alpha,beta
      real, dimension(N,N):: A,digamma
      real, dimension(N,S):: B
      integer, dimension(T) ::  observed
      integer :: i,j,k
      real, dimension(N) :: vecj
      digamma = 0.
      do k = 1,T-1
         do j = 1,N
            digamma(:,j) = digamma(:,j) + alpha(:,k)*A(:,j)*B(j,observed(k+1))*beta(j,k+1)
         end do
      end do

      digamma = digamma/T-1
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

! generates a random pi,A,B triplet for me.
subroutine make_random(pi,A,B,N,S)
      integer :: N
      integer :: S
      real, dimension(N) :: pi
      real, dimension(N,N) :: A
      real, dimension(N,S) :: B
      integer :: x,y

      call random_seed
      do x = 1, N
         call random_number(pi(x))
      end do
      pi = pi + 0.00000001
      pi = pi/sum(pi)

      do x = 1, N
         do y = 1, N
            call random_number(A(x,y))
         end do
         A(x,:) = A(x,:) +  0.00000001
         A(x,:) = A(x,:)/sum(A(x,:))
      end do

      do x = 1, N
         do y = 1, S
            call random_number(B(x,y))
         end do
         B(x,:) = B(x,:) +  0.00000001
         B(x,:) = B(x,:)/sum(B(x,:))
      end do
      end subroutine

!     does a single pass for me
subroutine single_pass(A,B,pi,observed,gamma,digamma,N,T,S)
      integer :: N
      integer :: T
      integer :: S

      real, dimension(N) :: pi
      real, dimension(N,N) :: A,digamma
      real, dimension(N,S) :: B
      real, dimension(N,T) :: gamma
      integer, dimension(T) :: observed

      real, dimension(N,T) :: alpha,beta
      real, dimension(T) :: sigma

      integer i
      call alphapass(alpha,sigma,A,B,pi,observed,N,T,S)
      call betapass(beta,sigma,A,B,observed,N,T,S)
      gamma = alpha*beta

      call digammapass(alpha,beta,observed,digamma,&
      & A,B,N,T,S)
      end subroutine

! re-estimates pi, A, and B from (gamma,observed,digamma)
subroutine reestimate_parameters(A,B,pi,observed,gamma,digamma,N,T,S)
      integer :: N
      integer :: S
      integer :: T

      real, dimension(N) :: pi
      real, dimension(N,N) :: A
      real, dimension(N,S) :: B
      integer, dimension(T) :: observed
      real, dimension(N,T) :: gamma
      real, dimension(N,N) :: digamma

      integer i,j
      pi = gamma(:,1)
      A = digamma
      B = 0.

      do i = 1, T
            B(:,observed(i)) = B(:,observed(i)) + gamma(:,i)
      end do

      do i = 1, N
         if (sum(B(i,:)) /= 0.) then
            B(i,:) = B(i,:)/sum(B(i,:))
         end if
      end do
    end subroutine reestimate_parameters

