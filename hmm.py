import numpy as np
import numba
import scipy.stats as st


def make_random(N,S):
    """generates a random pi,A,B triplet for me."""

    pi = np.random.dirichlet(np.ones(N))
    A  = np.random.dirichlet(np.ones(N),N).transpose()
    B  = np.random.dirichlet(np.ones(S),N)

    return pi,A,B


def generate_data(pi,A,B,T = 100):
    """Generates an observed sequence of length T"""
    N, S = B.shape

    # Generate hidden path
    X_t = np.zeros(T,dtype=np.int32)
    X_t[0] = np.argwhere(np.random.multinomial(1,pi,size=1).flatten()==1).flatten()[0]
    for t in range(1,T):
        X_t[t] = np.argwhere(np.random.multinomial(1,A[:,X_t[t-1]],size=1).flatten()==1).flatten()[0]

    # Generate observed sequence
    Y_t = np.zeros(T,dtype=np.int32)
    for t in range(T):
        Y_t[t] = np.argwhere(np.random.multinomial(1,B[X_t[t]],size=1).flatten()==1).flatten()[0]

    return X_t,Y_t


@numba.njit()
def alphapass(alpha,sigma,A,B,pi,observed):
    """
    alpha(t,i) = Pr(hidden_t = i | obs_1,...,obs_t)
    beta(t,i)  = Pr(obs_t+1,...,obs_T| hidden_t = i)*Pr(o_1,...,o_t)/Pr(o_1,...,o_T)
    gamma(t,i)  = Pr(hidden_t = i | obs_1,...,obs_T)
    A(j,i) == Pr( hidden i -> hidden j)
    B(i,j) == Pr( observed j | hidden i)
    """

    N,S = B.shape
    T = observed.shape[0]

    # initial iteration
    alpha[0] = pi
    alpha[0] = alpha[0]*B[:,observed[0]]
    sigma[0] = alpha[0].sum()
    alpha[0] /= sigma[0]

    # now do the remainder
    for t in range(1,T):
        alpha[t] = np.dot(A,alpha[t-1])
        alpha[t] *= B[:,observed[t]]
        sigma[t] = alpha[t].sum()
        alpha[t] /= sigma[t]


@numba.njit()
def betapass(beta,sigma,A,B,observed):
    """
    A(j,i) == Pr( hidden i -> hidden j)
    B(i,j) == Pr( observed j | hidden i)
    """

    N,S = B.shape
    T = observed.shape[0]

    At = A.transpose()

    # initial iteration
    beta[-1] = 1.

    # now do the remainder
    for t in range(T-1,0,-1):
        beta[t-1] = np.dot(At,beta[t]*B[:,observed[t]])/sigma[t]

@numba.njit()
def gammapass(alpha,beta,gamma):
    gamma[:,:] = alpha*beta


@numba.njit()
def update_B(B,gamma, observed):
    """ Updates belief on B """
    N, S = B.shape
    T = observed.shape[0]

    B_new = np.zeros((N,S))

    for t in range(T):
        B_new[:,observed[t]] += gamma[t]

    B_new = B_new/B_new.sum(axis=1).reshape(-1,1)
    B[:] = B_new

def update_pi(pi,gamma):
    """ Updates belief on pi """
    pi[:] = gamma[0]

@numba.njit()
def digammapass(alpha,beta,gamma,observed,sigma,digamma,A,B):
    """
    Computes an updated version of A based on the observed sequence of data.

    alpha(t,i) = Pr(hidden_t = i | obs_1,...,obs_t)
    beta(t,i)  = Pr(obs_t+1,...,obs_T| hidden_t = i)*Pr(o_1,...,o_t)/Pr(o_1,...,o_T)
    gamma(t,i)  = Pr(hidden_t = i | obs_1,...,obs_T)
    A(j,i) = Pr( hidden i -> hidden j)
    B(i,j) = Pr( observed j | hidden i)

    Reestimates Pr(X=i -> X=j) via averaging
    Pr(x_t=i,x_t+1=j|o_1,...,o_T) = alpha_t(i)A_jibeta_{t+1}(j)B(j,o_{t+1})
    """

    T = alpha.shape[0]
    N,S = B.shape

    digamma[:] = 0
    digammatemp = np.zeros((N,N))
    for t in range(0,T-1):
        temp_i = (alpha[t]/sigma[t+1]).reshape(1,-1)*A
        temp_j = (B[:,observed[t+1]]*beta[t+1]).reshape(-1,1)
        digammatemp = temp_i*temp_j
        digamma[:] += digammatemp

    # rescale columns as \sum_j A_{j,i} == 1
    digamma[:] = digamma/digamma.sum(axis=0).reshape(1,-1)

def update_parameters(alpha,beta,gamma, pi, A,B,observed,sigma):
    """ Updates pi, A,B from observed data. """

    N,S = B.shape
    digamma = np.zeros((N,N))
    digammapass(alpha, beta, gamma, observed,sigma,digamma, A,B)

    update_B(B,gamma,observed)
    update_pi(pi,gamma)

    A[:] = digamma[:]

@numba.njit()
def gammaoneshot(alpha,beta,gamma,pi,A,B,observed,sigma):
    """ computes alpha, beta, and gamma directly for me. """

    alphapass(alpha,sigma,A,B,pi,observed)
    betapass(beta,sigma,A,B,observed)
    gammapass(alpha,beta,gamma)

def compute_loglikelihood(sigma):
    """Computes Pr(Y_0,...,Y_T| pi, A, B) from (sigma_t)"""
    return np.log(sigma).sum()

def compute_loglikelihood_parameters(pi,A,B,observed):
    """Computes Pr(Y_0,...,Y_T| pi, A, B) from (pi, A, B)"""
    N, S = B.shape
    T = observed.shape[0]
    alpha = np.zeros((T,N))
    beta  = np.zeros((T,N))
    gamma = np.zeros((T,N))
    sigma = np.zeros(T)

    gammaoneshot(alpha,beta,gamma,pi,A,B,observed,sigma)
    score = compute_loglikelihood(sigma)
    return score

def reestimate_parameters(A,B,pi,observed,
                          halting_criteria=1e-6,debug=False):
    """re-estimates pi, A, and B from the observed sequence."""
    N, S = B.shape
    T = observed.shape[0]
    alpha = np.zeros((T,N))
    beta  = np.zeros((T,N))
    gamma = np.zeros((T,N))
    sigma = np.zeros(T)

    current_score = -np.inf
    gammaoneshot(alpha,beta,gamma,pi,A,B,observed,sigma)
    new_score = compute_loglikelihood(sigma)
    diff = new_score - current_score
    while diff > halting_criteria:
        current_score = new_score
        update_parameters(alpha,beta,gamma, pi, A,B,observed,sigma)
        gammaoneshot(alpha,beta,gamma,pi,A,B,observed,sigma)
        new_score = compute_loglikelihood(sigma)
        diff = new_score - current_score
        if debug:
            print("diff = ", diff)

def determine_number_of_hidden_states(Y_t,halting_criteria=1e-6,verbose=False):
    """Uses Wilks' theorem to determine the number of hidden states based on
    the observed data (Y_t). It returns the number of hidden states"""
    S = Y_t.max() + 1
    N = 1
    my_pi, my_A, my_B = make_random(N,S)

    reestimate_parameters(my_A,my_B, my_pi,Y_t,halting_criteria=halting_criteria)
    S_alternative =compute_loglikelihood_parameters(my_pi, my_A, my_B, Y_t)

    p = 0.0
    while p < 0.01:
        df = 2*N + S
        N += 1
        my_pi, my_A, my_B = make_random(N,S)

        S_null = S_alternative
        reestimate_parameters(my_A,my_B, my_pi,Y_t,halting_criteria=halting_criteria)
        S_alternative =compute_loglikelihood_parameters(my_pi, my_A, my_B, Y_t)

        model =st.chi2(df = df)
        lam = 2*(S_alternative - S_null)
        p = model.sf(lam)
        if verbose:
            print("N=",N,"df: ", df, "Lamda: ", lam, "p: ", p)
    return N-1

class hmm():
    """This class is intended to mimic the machine learning algorithms implemented in
    Scikit-Learn. It's intended to make it as easy as possible to recover the
    parameters pi, A, B of the HMM model as well as to project the observed data into the
    hidden state space using the gamma distribution.

    Call the fit() method with observed data to recover the pi, A, B parameters from
    an observed sequence of values. """

    def __init__(self, num_hiddenstates = None, halting_criteria=1e-6,
                 pi = None, A = None, B = None):
        """My initializer. If you know the number of hidden states, pass it
        off via num_hiddenstates.

        halting_criteria is the parameter for specifying when to halt updating the values of pi, A, B.
        """
        self.N = num_hiddenstates
        self.halting_criteria = halting_criteria
        self._called_fit = False

        self.pi = pi
        self.A  = A
        self.B  = B
        if not self.N is None:
            if self.A is None:
                self.A = np.random.dirichlet(np.ones(self.N),size=self.N).transpose()

            if self.pi is None:
                self.pi = np.random.dirichlet(np.ones(self.N),size=1).flatten()

    def fit(self, observed, verbose=False):
        """ The method to call to find the values of pi, A, B from the observed data.
        observed should be a list of observed values that are hashable.

        To map back from the index B[i,j] = Pr(Y=j | X= i) to the observed value Y,
        use the dictionary self.inv_map[j] = observed value"""

        unique = list(set([y for y in observed]))
        unique.sort()
        self.map = dict([(y,t) for t,y in enumerate(unique)])
        self.inv_map = dict([(t,y) for t,y in enumerate(unique)])
        self.Y_t = np.array([self.map[y] for y in observed])
        self.S = len(self.map)
        self.T = len(self.Y_t)

        if self.N is None:
            self.N =  determine_number_of_hidden_states(self.Y_t,
                                                        halting_criteria=self.halting_criteria,
                                                        verbose=verbose)

            self.pi, self.A, self.B = make_random(self.N,self.S)
        else:
            if self.B is None:
                self.B = np.random.dirichlet(np.ones(self.S),size=self.N)

        reestimate_parameters(self.A,self.B,self.pi,self.Y_t,
                              halting_criteria=self.halting_criteria)
        self._called_fit = True
    def transform(self):
        """Takes the observed data mapped to (Y_t) and then projected to the
        gamma sequence. Also computes the log likelihood of the observed data.
        Returns the gamma sequence. """

        if not self._called_fit:
            raise Exception("You need to call the fit method prior to calling the transform method.")

        alpha = np.zeros((self.T,self.N))
        beta  = np.zeros((self.T,self.N))
        gamma = np.zeros((self.T,self.N))
        sigma = np.zeros(self.T)
        gammaoneshot(alpha,beta,gamma,self.pi,self.A,self.B,self.Y_t,sigma)
        self.X_t = gamma
        self.logp = compute_loglikelihood(sigma)
        return gamma
    def fit_transform(self,observed,verbose=False):
        """Fits the model to the observed data and returns the projection of the observed
        data into the hidden state space via the gamma sequence, Y_t -> gamma_t"""

        self.fit(observed,verbose=verbose)
        return self.transform()

    def compute_loglikelihood(self):
        """ Computes log Pr(Y_0,...,Y_T| pi, A, B) """
        self.logp = compute_loglikelihood_parameters(self.pi,self.A,self.B,self.Y_t)
        return self.logp


