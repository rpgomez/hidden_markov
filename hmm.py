import numpy as np
import numba

def alphapass(alpha,sigma,A,B,pi,observed):
    """
    alpha(i,t) = Pr(hidden_t = i | obs_1,...,obs_t)
    beta(i,t)  = Pr(obs_t+1,...,obs_T| hidden_t = i)*Pr(o_1,...,o_t)/Pr(o_1,...,o_T)
    gamma(i,t)  = Pr(hidden_t = i | obs_1,...,obs_T)
    A(j,i) == Pr( hidden i -> hidden j)
    B(i,j) == Pr( observed j | hidden i)
    """

    N = A.shape[0]
    T = observed.shape[0]
    S = B.shape[1]

    # initial iteration
    alpha[:,0] = pi
    alpha[:,0] = alpha[:,0]*B[:,observed[0]]
    sigma[0] = alpha[:,0].sum()
    alpha[:,0] /= sigma[0]

    # now do the remainder
    for t in range(1,T):
        alpha[:,t] = A.dot(alpha[:,t-1])
        alpha[:,t] *= B[:,observed[t]]
        sigma[t] = alpha[:,t].sum()
        alpha[:,t] /= sigma[t]



def betapass(beta,sigma,A,B,observed):
    """
    A(j,i) == Pr( hidden i -> hidden j)
    B(i,j) == Pr( observed j | hidden i)
    """

    N = A.shape[0]
    T = observed.shape[0]
    S = B.shape[1]

    At = A.transpose()

    # initial iteration
    beta[:,-1] = 1.

    # now do the remainder
    for t in range(T-1,0,-1):
        beta[:,t-1] = At.dot(beta[:,t]*B[:,observed[t]])/sigma[t]


def digammapass(alpha,beta,gamma,observed,sigma,digamma,A,B):
    """
    Reestimates Pr(i -> j)
    Pr(x_t=i,x_t+1=j|o_1,...,o_T) = alpha_t^(i)A_jibeta_^t+1(j)b(j,o_t+1)
    """
      N,T = alpha.shape
      S = B.shape[1]
      digamma[:] = 0
      digammatemp = np.zeros((N,N))
      for t in range(0,T-1):
         digammatemp[:] = 0
         for j in range(N):
            digammatemp[:,j] = alpha[:,t]*A[j,:]*B[j,observed[t+1]]*beta[j,t+1]

         digammatemp = digammatemp/sigma[t+1]
         digamma[:] += digammatemp

      for n in range(N):
         digamma[:,n] = digamma[:,n]/gamma[n,:T-1].sum()

def gammapass(alpha,beta,gamma):
    N,T = alpha.shape

    gamma[:,:] = alpha*beta


def gammaoneshot(gamma,pi,A,B,observed,sigma):
    """ computes gamma directly for me. """

    N = A.shape[0]
    S = B.shape[1]
    T = observed.shape[0]

    alphapass(alpha,sigma,A,B,pi,observed)
    betapass(beta,sigma,A,B,observed)
    gammapass(alpha,beta,gamma)

def make_random(N,S):
    """generates a random pi,A,B triplet for me."""

      real, dimension(N) :: pi
      real, dimension(N,N) :: A
      real, dimension(N,S) :: B
      integer :: x,y
      pi = np.random.(np.ones(N))
      A  = np.random.(np.ones(N),N).transpose()
      B  = np.random.(np.ones(S),N)

      return pi,A,B
  
def single_pass(A,B,pi,observed,gamma,digamma):
    """computes alpha, beta, gamma, digamma passes """

    gammaoneshot(gamma,pi,A,B,observed,sigma)
    digammapass(alpha,beta,gamma,observed,sigma,digamma)


def reestimate_parameters(A,B,pi,observed,halting_criteria=1e-6):
    """re-estimates pi, A, and B from observed"""


      do tt = 1, iters
         call single_pass(A,B,pi,observed,gamma,digamma,N,T,S)
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
      end do
    end subroutine reestimate_parameters
