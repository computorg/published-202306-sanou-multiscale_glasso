import sys

# install the developmental version of pylearn-parsimony forked and available on 
# https://github.com/desanou/pylearn-parsimony/tree/mglasso_integration
# pip install git+https://github.com/desanou/pylearn-parsimony.git@mglasso_integration

import parsimony.algorithms as algorithms
import parsimony.estimators as estimators
import parsimony.functions as functions
import parsimony.functions.nesterov.tv as tv
from parsimony.algorithms.utils import Info
import random
import numpy as np
import sklearn.preprocessing
from scipy import sparse
from scipy.linalg import block_diag

from parsimony.functions.losses import MGLassoRegression
from parsimony.algorithms.subgradient import SubGradientDescent
from parsimony.algorithms.utils import NonSumDimStepSize

def Ak_from_pairs(k,p):
    Ak = sparse.lil_matrix((int(p*(p-1)/2),p*p))
    ij=0

    for i in range(0,p-1):
        for j in range(i+1,p):
            if i==k:
                Ak[ij,i*p+k]=1
                Ak[ij,j*p+j]=-1
            elif j==k:
                Ak[ij,i*p+j]=1
                Ak[ij,j*p+i]=-1
            else:
                Ak[ij,i*p+k]=1
                Ak[ij,j*p+k]=-1
            ij=ij+1
    to_keep = list(set(range(Ak.shape[1]))-set(range(0,p*p,p+1)))
    Aknew = sparse.lil_matrix(sparse.csr_matrix(Ak)[:,to_keep])
    return(Aknew)


def linear_operator_from_num_variables(num_variables, type_, W):

    """Generates the linear operator for the TV lasso Nesterov function
    from number of variables.

    Parameters:
    ----------
    num_variables : Integer. The total number of variables, including the
            intercept variable(s).

    """
    A = list()
    for k in range(0,num_variables):
        Ak = Ak_from_pairs(k,num_variables)
        A.append(Ak.tocsr())
    return A

def beta2Beta(beta,p):
    Beta=np.zeros((p,p))
    for j in range(0,(p-1)):
        for i in range(0,p):
            k=i
            l=j
            if j>=i:
                l=j+1
            Beta[k,l]=beta[i*(p-1)+j]
    return(Beta)

def precision2regression(K):
    p=K.shape[0]
    M=np.zeros((p,p))
    for i in range(0,p):
        for j in range(0,p):
            if i!=j:
                M[i,j]= - K[i,j]/K[i,i]
    return(M)

def admm_algo(X, lam1, lam2, max_iter, prec, rho = 1, mean=False, mu = 0, simulation = False, tau_incr = 2, tau_decr = 2):
    X = np.array(X)
    n = X.shape[0]
    p = X.shape[1]

    X = sklearn.preprocessing.scale(X)
    y = X.reshape(n * p, 1, order='F')
    Xvec = np.delete(np.kron(np.identity(p), X), range(0, p * p, p + 1), axis=1)
    A_ = linear_operator_from_num_variables(p, "initial", None)

    info = [Info.converged, Info.num_iter, Info.time, Info.fvalue, Info.ok, Info.admm_iter, Info.fista_iter]

    x_ = np.zeros((A_[0].shape[1], 1))
    Ax = [0.0] * len(A_)
    for i in range(len(A_)):
        Ax[i] = A_[i].dot(x_)
    Ax = np.vstack(Ax)
    dimr = Ax.shape[0]

    algorithm = algorithms.proximal.ADMMglasso(info=info, eps=prec, max_iter=max_iter, A=A_, mu = mu, simulation = simulation, tau_incr = tau_incr, tau_decr = tau_decr, rho = rho)
    function = functions.combinedfunctions.AugmentedLinearRegressionMglasso(Xvec, y, lam1, lam2, A=A_,
                                                                                      mean=mean, rho = rho)

    pp = Xvec.shape[1]
    xa = [np.zeros((pp, 1)), np.zeros((dimr, 1))]
    beta = algorithm.run(function, xa, A_, lam2)
    info_admm = algorithm.info_get()

    out = {
        "method": "ADMM",
        "time" : info_admm[Info.time],
        "time_cumu": np.cumsum(info_admm[Info.time]),
        "niter_cumu": np.cumsum(info_admm[Info.fista_iter]),
        "fvalue": info_admm[Info.fvalue],
        "converged": info_admm[Info.converged],
        "beta": beta2Beta(beta,p),
        "niter": info_admm[Info.num_iter],
        "admm_iters": info_admm[Info.admm_iter]
    }

    return (out)

def conesta_algo(X, lam1, lam2, max_iter, prec, mean = False):
    X=np.array(X)
    n=X.shape[0]
    p=X.shape[1]

    y=X.reshape(n*p,1,order='F')
    Xvec=np.delete(np.kron(np.identity(p),X),range(0,p*p,p+1),axis=1)
    A_=linear_operator_from_num_variables(p, "initial", None)

    info = [Info.converged, Info.num_iter, Info.time, Info.fvalue, Info.ok, Info.func_val, Info.continuations, Info.fista_iter]

    hgmm_conesta = estimators.LinearRegressionL1L2TV(l1 = lam1, l2 = 0.0, tv = lam2, A = A_,
                                           algorithm=algorithms.proximal.CONESTA(info = info, eps=prec, max_iter=max_iter), mean=mean)

    res_conesta = hgmm_conesta.fit(Xvec, y)
    info_conesta = res_conesta.get_info()

    out = {
        "method": "CONESTA",
        "time" : info_conesta[Info.time],
        "time_cumu": np.cumsum(info_conesta[Info.time]),
        "niter_cumu": np.cumsum(info_conesta[Info.fista_iter]),
        "fvalue": info_conesta[Info.fvalue],
        "converged": info_conesta[Info.converged],
        "beta": beta2Beta(res_conesta.beta,p),
        "niter": info_conesta[Info.num_iter],
        "continuations": info_conesta[Info.continuations]
    }

    return(out)

def sgd_algo(X, lam1, lam2, max_iter, prec):
    X=np.array(X)
    n=X.shape[0]
    p=X.shape[1]

    X=sklearn.preprocessing.scale(X)
    y=X.reshape(n*p,1,order='F')
    Xvec=np.delete(np.kron(np.identity(p),X),range(0,p*p,p+1),axis=1)
    A_=linear_operator_from_num_variables(p, "initial", None)

    info_sgd = [Info.converged, Info.num_iter, Info.time, Info.ok, Info.func_val]

    sgd = SubGradientDescent(max_iter=max_iter, step_size=NonSumDimStepSize(a=0.1),
                         use_gradient=False, info=info_sgd, eps=prec, use_best_f=True)

    pp = Xvec.shape[1]
    function = MGLassoRegression(Xvec, y, l1 = lam1, tv = lam2, A = A_, mean=False)
    beta = sgd.run(function, np.zeros((pp, 1)))

    info_sgd = sgd.info_get()
    niter = info_sgd[Info.num_iter]

    out = {
        "method": "SGD",
        "time" : info_sgd[Info.time],
        "time_cumu": np.cumsum(info_sgd[Info.time]),
        "niter_cumu": np.cumsum(np.repeat(1, niter)),
        "fvalue": info_sgd[Info.func_val],
        "converged": info_sgd[Info.converged],
        "beta": beta2Beta(beta, p),
        "niter": niter
    }

    return out
