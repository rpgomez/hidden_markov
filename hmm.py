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
    alpha(t,i) = Pr(hidden_t = i | obs_1,...,obs_t)
    beta(t,i)  = Pr(obs_t+1,...,obs_T| hidden_t = i)*Pr(o_1,...,o_t)/Pr(o_1,...,o_T)
    gamma(t,i)  = Pr(hidden_t = i | obs_1,...,obs_T)
    A(j,i) == Pr( hidden i -> hidden j)
    B(i,j) == Pr( observed j | hidden i)

    Reestimates Pr(i -> j)
    Pr(x_t=i,x_t+1=j|o_1,...,o_T) = alpha_t(i)A_jibeta_{t+1}(j)b(j,o_t+1)
    """

    N,T = alpha.shape
    S = B.shape[1]
    digamma[:] = 0
    digammatemp = np.zeros((N,N))
    for t in range(0,T-1):
        digammatemp[:] = 0
        for j in range(N):
            digammatemp[j,:] = alpha[:,t]*A[j,:]*B[j,observed[t+1]]*beta[j,t+1]

        digammatemp = digammatemp/sigma[t+1]
        digamma[:] += digammatemp

    digamma = digamma/digamma.sum(axis=0).reshape(1,-1)
    #for n in range(N):
    #    digamma[:,n] = digamma[:,n]/gamma[n,:T-1].sum()

def gammapass(alpha,beta,gamma):
    N,T = alpha.shape

    gamma[:,:] = alpha*beta

def update_B(B,gamma, observed):
    """ Updates Belief on B """
    N, S = B.shape
    T = observed.shape[0]

    B_new = np.zeros((N,S))

    for t in range(T):
        B_new[:,observed[t]] += gamma[:,t]

    B_new = B_new/B_new.sum(axis=1).reshape(-1,1)
    B[:] = B_new

def update_pi(pi,gamma):
    """ Updates belief on pi """
    pi[:] = gamma[:,0]
    
def update_parameters(alpha,beta,gamma, pi, A,B,observed,sigma):
    """ Updates pi, A,B from observed data. """

    N,S = B.shape
    digamma = np.zeros((N,N))
    digammapass(alpha, beta, gamma, observed,sigma,digamma, A,B)

    update_B(B,gamma,observed)
    update_pi(pi,gamma)

    A[:] = digamma[:]
    
def gammaoneshot(gamma,pi,A,B,observed,sigma):
    """ computes gamma directly for me. """

    alphapass(alpha,sigma,A,B,pi,observed)
    betapass(beta,sigma,A,B,observed)
    gammapass(alpha,beta,gamma)

def make_random(N,S):
    """generates a random pi,A,B triplet for me."""

    pi = np.random.dirichlet(np.ones(N))
    A  = np.random.dirichlet(np.ones(N),N).transpose()
    B  = np.random.dirichlet(np.ones(S),N).transpose()

    return pi,A,B
  
def single_pass(A,B,pi,observed,gamma,digamma):
    """computes alpha, beta, gamma, digamma passes and returns log
    likelihood."""

    gammaoneshot(gamma,pi,A,B,observed,sigma)
    digammapass(alpha,beta,gamma,observed,sigma,digamma)

def compute_likelihood(sigma):
    return log(sigma).sum()

def reestimate_parameters(A,B,pi,observed,halting_criteria=1e-6,debug=False):
    """re-estimates pi, A, and B from observed"""

    N, S = B.shape
    T = observed.shape[0]
    alpha = zeros((N,T))
    beta  = zeros((N,T))
    gamma = zeros((N,T))
    sigma = zeros(T)

    current_score = -np.inf
    gammaoneshot(gamma,pi,A,B,observed,sigma)
    new_score = compute_likelihood(sigma)
    diff = new_score - current_score
    while diff > halting_criteria:
        current_score = new_score
        update_parameters(alpha,beta,gamma, pi, A,B,observed,sigma)
        gammaoneshot(gamma,pi,A,B,observed,sigma)
        new_score = compute_likelihood(sigma)
        diff = new_score - current_score
        if debug:
            print("diff = ", diff)

    
