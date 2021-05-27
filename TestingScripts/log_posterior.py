from numba import njit
from math import ceil

@njit
def func_invKgYinvKt(Ng,Nt,lt,invKgY,invKgYinvKt):

    u = np.exp(-1.0/lt)
    u2 = u*u
    oneplusu2 = 1.0 + u2
    oneoveroneminusu2 = 1.0 / ( 1.0 - u2 )
    
    for ig in range(Ng):
        invKgYinvKt[ig,0] = (invKgY[ig,0] - u * invKgY[ig,1])*oneoveroneminusu2
        invKgYinvKt[ig,-1] = (invKgY[ig,-1] - u * invKgY[ig,-2])*oneoveroneminusu2
        for it in range(1,Nt-1):
            invKgYinvKt[ig,it] = ( oneplusu2 * invKgY[ig,it] - u * ( invKgY[ig,it-1] + invKgY[ig,it+1] ) )*oneoveroneminusu2
            
            
@njit
def func_YJt(Ng,Nt,lt,Y,result,log_eps=np.log(1e-16),safety_margin=10):
    
    u = np.exp(-1/lt)
    u2 = u*u
    
    M = safety_margin+ceil(-lt*log_eps)
    
    power_u = np.power(u,np.arange(M))
    
    lt2 = lt*lt
    
    for i in range(Ng):
        
        for l in range(Nt):
            
            res =  - (1.0+u2)*Y[i,l]/(1.0-u2)
            
            if l < M:
                res += u2*Y[i,0]*power_u[l]/(1.0-u2)
                
            if Nt-l-1 < M:
                res += u2*Y[i,-1]*power_u[Nt-l-1]/(1.0-u2)
            
            #res = (u2*( Y[i,0]*power_u[min(M,l)] + Y[i,-1]*power_u[min(M,Nt-l-1)] ) - (1.0+u2)*Y[i,l])/(1.0-u2)
            
            for j in range(max(0,l-M),min(Nt,l+M)):
                res += Y[i,j]*power_u[abs(l-j)]
                
            result[i,l] = res/lt2
    

def log_prior_x(x,mu,lt,lg,sigma2,G,Ng,Nt):
    
    # X is a (Ng,Nt) matrix
    # mu is a (Ng,) vector
    # sigma2, lt and lg are scalars
    # G2 is a (Ng,Ng) matrix
    
    Y = x.reshape((Ng,Nt)) - mu.reshape((Ng,1))
    
    G2 = np.square(G[:,np.newaxis]-G[np.newaxis,:])
    Kg = np.exp(-G2/lg/lg/2.0)
    dKg_dlg = G2*Kg/lg/lg/lg
    Jg = np.linalg.lstsq(np.dot(Kg.T,Kg),np.dot(Kg.T,dKg_dlg),rcond=None)[0]
    invKgY = np.linalg.lstsq(np.dot(Kg.T,Kg),np.dot(Kg.T,Y),rcond=None)[0]
    
    
    invKgYinvKt = np.zeros((Ng,Nt))
    func_invKgYinvKt(Ng,Nt,lt,invKgY,invKgYinvKt)
    
    u = np.exp(-1.0/lt)
    u2 = u*u
    logdetKg = np.linalg.slogdet(Kg)[1]
    logdetinvKt = -(Nt-1.0)*np.log( 1.0 - u2 )
    TrJg =  np.trace(Jg)
    TrJt = -2.0*(Nt-1.0)*u2/(1.0-u2)/lt/lt
    
    Y_invKgYinvKt = (Y*invKgYinvKt).sum()
    
    JgTY = np.dot(Jg.T,Y)
    JgTY_invKgYinvKt = (JgTY*invKgYinvKt).sum()
    
    YJt = np.zeros((Ng,Nt))
    func_YJt(Ng,Nt,lt,Y,YJt,log_eps=np.log(1e-16))
    YJt_invKgYinvKt = (YJt*invKgYinvKt).sum()
    
    lnP = -Ng*Nt*np.log(2.0*np.pi*sigma2)/2.0 + Ng*logdetinvKt/2.0 - Nt*logdetKg/2.0 - Y_invKgYinvKt/2.0/sigma2
    
    dlnP_dX = -invKgYinvKt/sigma2
    dlnP_dx = dlnP_dX.ravel()
    dlnP_dmu = -dlnP_dX.sum(axis=1)
    dlnP_dsigma2 = -Ng*Nt/2.0/sigma2 + Y_invKgYinvKt/2.0/sigma2/sigma2
    dlnP_dlt = -Ng*TrJt/2.0 + YJt_invKgYinvKt/2.0/sigma2
    dlnP_dlg = -Nt*TrJg/2.0 + JgTY_invKgYinvKt/2.0/sigma2
    #print(invKgY)
    #print(Kg,Y)
    
    return lnP, dlnP_dx, dlnP_dmu, dlnP_dlt, dlnP_dlg, dlnP_dsigma2

def log_prior_mu(mu,m,tau2):
    
    J = tau2*np.exp(-G_squared/2.0/m**2)
    J_inv = np.linalg.pinv(J)
    J_inv_mu = np.dot(J_inv,mu)
    
    lnQ = +0.5*np.linalg.slogdet(J_inv)[1]-0.5*np.dot(mu.T,J_inv_mu) - (M/2)*np.log(2.0*np.pi)
    
    dJdm = J*G2/m**3
    dJdtau2 = J/tau2

    dlnQdm = -0.5*np.trace(np.dot(J_inv,dJdm))+0.5*np.dot(J_inv_mu.T,np.dot(dJdm,J_inv_mu))
    dlnQdtau2 = -0.5*np.trace(np.dot(J_inv,dJdtau2))+0.5*np.dot(J_inv_mu.T,np.dot(dJdtau2,J_inv_mu))
    dlnQdmu = -J_inv_mu
    
    return lnQ, dlnQdmu, dlnQdm, dlnQdtau2

def log_prior_lengthscale(l):
    
    # InvGamma(1,2)
    
    lnA = np.log(2.0)-2.0*np.log(l)-2.0/l
    
    dlnAdl = 2.0*(1.0-l)/l/l
    
    return lnA, dlnAdl

def log_prior_variance(sigma2):
    
    # Gamma(1,1)
    
    lnS = - sigma2
    
    dlnSdsigma2 = - 1.0
    
    return lnS, dlnSsigma2

def log_likelihood(p):
    
    return 0.0, np.zeros(p.shape)

def log_posterior(args,Nh,Ng,Nt,G):
    
    # Unpack parameters
    ln_lt, ln_lg, ln_sigma2, ln_m, ln_tau2 = args[:Nh]
    mu = args[Nh:Nh+Ng]
    x = args[Nh+Ng:]
    
    # Transform parameters
    lt, lg, sigma2, m, tau2 = np.exp(ln_lg), np.exp(ln_lt), np.exp(ln_sigma2), np.exp(ln_m), np.exp(ln_tau2)
    p = 1.0/(1.0+np.exp(-x))
    
    ##### Priors
    
    ### Initialise values
    
    lnF = 0.0
    dlnF_dlt, dlnF_dlg, dlnF_dm = 0.0, 0.0
    dlnF_dsigma2, dlnF_dtau2 = 0.0, 0.0
    dlnF_dmu = np.zeros(Ng)
    dlnF_dx = np.zeros(Ng*Nt)
    dlnF_dp = np.zeros(Ng*Nt)
    
    
    ### Prior on x
    
    # Calculate quantities
    lnP, dlnP_dx, dlnP_dmu, dlnP_dlt, dlnP_dlg,dlnP_dsigma2 = log_prior_x(x,mu,lt,lg,sigma2,G,Ng,Nt)

    # Increment log posterior
    lnF += lnP

    # Increment gradients
    dlnF_dx += dlnP_dx
    dlnF_dmu += dlnP_dmu
    dlnF_dlt += dlnP_dlt
    dlnF_dlg += dlnP_dlg
    dlnF_dsigma2 += dlnP_dsigma2
    
    ### Prior on mu
    
    # Calculate quantities
    lnQ, dlnQ_dmu, dlnQ_dm, dlnQ_dtau2 = log_prior_mu(mu,m,tau2)
    
    # Increment log posterior
    lnF += lnQ
    
    # Increment gradients
    dlnF_dmu += dlnQ_dmu
    dlnF_dm += dlnQ_dm
    dlnF_dtau2 += dlnQ_dtau2
    
    ### Prior on lt, lg and m
    
    # Calculate quantities
    lnA, dlnA_dlt = log_prior_lengthscale(lt)
    lnB, dlnB_dlg = log_prior_lengthscale(lg)
    lnC, dlnC_dm  = log_prior_lengthscale(m)
    
    # Increment log posterior
    lnF += lnA
    lnF += lnB
    lnF += lnC
    
    # Increment gradients
    dlnF_dlt += dlnA_dlt
    dlnF_dlg += dlnB_dlg
    dlnF_dm  += dlnC_dm
    
    ### Prior on sigma2 and tau2
    
    # Calculate quantities
    lnS, dlnS_dsigma2 = log_prior_variance(sigma2)
    lnT, dlnT_dtau2 = log_prior_variance(tau2)
    
    # Increment log posterior
    lnF += lnS
    lnF += lnT
    
    # Increment gradients
    dlnF_dsigma2 += dlnS_dsigma2
    dlnF_dtau2 += dlnT_dtau2
    
    ##### Likelihood
    
    # Calculate quantities
    lnL, dlnL_dp = log_likelihood(p)
        
    # Increment log posterior
    lnF += lnL
    
    # Increment gradients
    dlnF_dp += dlnL_dp
    
    ##### Construct gradient with respect to input parameters
    
    # Initialise
    dlnF_dargs = np.zeros(Nh+Ng+Ng*Nt)
    
    # Correct for log-parameterisation
    dlnF_dln_lt = lt*dlnF_dlt
    dlnF_dln_lg = lg*dlnF_dlg
    dlnF_dln_m = m*dlnF_dm
    dlnF_dln_sigma2 = sigma2*dlnF_dsigma2
    dlnF_dln_tau2 = tau2*dlnF_dtau2
    
    # Correct for logit-parameterisation
    dlnF_dx += (1.0/p+1.0/(1.0-p))*dlnF_dp
    
    # Fill in values
    dlnF_dargs[0] = dlnF_dln_lt
    dlnF_dargs[1] = dlnF_dln_lg
    dlnF_dargs[2] = dlnF_dln_sigma2
    dlnF_dargs[3] = dlnF_dln_m
    dlnF_dargs[4] = dlnF_dln_tau2
    dlnF_dargs[Nh:Nh+M] = dlnF_dmu
    dlnF_dargs[Nh+M:] = dlnF_dx
    
    return lnF, dlnFdargs

