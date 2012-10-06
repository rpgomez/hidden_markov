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
      do 10 i = 1,T-1
         alpha(:,i) = matmul(A,alpha(:,i-1))
         alpha(:,i) = alpha(:,i)*B(observed(i),:)
         sigma(i) = sum(alpha(:,i))
         alpha(:,i) = alpha(:,i)/sigma(i)
 10   continue

end
