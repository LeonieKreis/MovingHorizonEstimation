from casadi import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import random
import time

def simulate_pendulum(N,T,L,p_true,s0,x,p,rhs_exp):
    x0 = s0 ## initial conditions from paper
    Res1 = s0
    dae = {'x': x, 'p': p, 'ode': rhs_exp}
    opts = {'tf': T/N}
    F = integrator('F', 'cvodes', dae, opts)
    for i in range(L*N):
        Fi = F(x0=s0, p=p_true)
        Xk_end = Fi['xf']
        # for k in range(1,M+1):
        Res1 = vertcat(Res1,Xk_end[:,-1])
        s0 = Xk_end[:,-1]
    return Res1
    

def generate_meas(x, p, xdot, T, sigma, p_true, M, x0):
    """
    x = an MX.sym, state of the ODE
    p = an MX.sym, parameters of the ODE
    xdot = an MX.sym, righthandside of the ODE
    T = Final Time 
    sigma = variance of the parameters, eg.g [1,1,2,1] 
    p_true = the true parameters of the problem
    M+1 = number of measurements
    x0 = array or list, initial value for the ODE 
    
    return: gaussian perturbed measurements (of the states!) for Multiple shooting nodes
    """
    n = len(x0) #x0.shape[0]
    #print(n)
    meas_times = np.linspace(0, T, num=M+1)
    
    # Computing measurements with true parameters
    dae = {'x': x, 'p': p, 'ode': xdot} 
    opts = {'grid' : meas_times}  
    F = integrator('F', 'cvodes', dae, opts)
    solution = F(x0=x0,p=p_true)['xf']
    x_0 = DM(x0)
    meas_exact = horzcat(x_0, solution)
    
    # adding noise to measurements
    noise = []
    for i in range(M+1):
        mean = 0#np.zeros(n)
        #sigma = [1,2,3,4]
        noise=horzcat(noise,np.random.normal(mean, sigma))#, (n,1))

    meas = meas_exact + noise
    return(meas)

def arrival_cost_values(x,p,xdot,T,N,x_opt, p_opt, last_y, last_P,last_V, last_W):
    #print('Arrival cost not implemented yet!')
    n = x.shape[0]
    m = p.shape[0]
    
    x_now = x_opt[0:n] # x*(t_L)
    #print('x_now',type(x_now))
    p_now = p_opt # p*
    #print('p_now',type(p_now))
    #compute x*(t_L+1)
    #print(type(x))
    #print(type(p))
    dae = {'x':x, 'p':p, 'ode':xdot} 
    opts = {'t0':0, 'tf':T/N}  
    F = integrator('F', 'cvodes', dae, opts)
    #print(F)
    x_plus1 = F(x0=x_now,p=p_now)['xf'] #x*(t_L+1)
    # compute X_x, X_p
    x00 = MX.sym('x00',n)
    p00 = MX.sym('p00',m)
    X_x_exp = jacobian(F(x0=x00,p=p00)['xf'],x00)
    X_p_exp = jacobian(F(x0=x00,p=p00)['xf'],p00)
    X_x = Function('X_x',[x00,p00],[X_x_exp])
    X_p = Function('X_p',[x00,p00],[X_p_exp]) # x00 later x_plus1 p00 later p_now
    H_x = X_x
    H_p = X_p
    #print('X_x',X_x(x_plus1,p_now).shape)
    #print('X_p',X_p(x_plus1,p_now).shape)
    # setup M
    zeros_or = np.zeros((n+m,n+m))
    zeros_mr = np.zeros((n,n+m))
    c2 = np.concatenate((np.concatenate((zeros_or,zeros_mr)),last_W))
    #print(c2.shape)
    m_l = -1*(horzcat(last_V@X_x(x_plus1,p_now),last_V@X_p(x_plus1,p_now)))
    #print(m_l.shape)
    b_row1 = np.concatenate((X_x(x_plus1,p_now),X_p(x_plus1,p_now)),axis=1)
    #print(b_row1.shape)
    ############b_row22 = np.concatenate((np.zeros(n),np.ones(m)))
    #print(b_row22.shape)
    ##########b_row2 = np.reshape(b_row22,(1,m+n))
    b_row2 = np.concatenate((np.zeros((m,n)),np.diag(np.ones(m))),axis=1)
    #print('brow2',b_row2.shape)
    blocks = np.concatenate((b_row1,b_row2))
    #print(blocks.shape)
    lo_l = -1*(last_W@blocks)
    c1 = vertcat(vertcat(last_P,m_l),lo_l)
    M_num = horzcat(c1,c2)
    M = MX(M_num)#SX(M_num)
    #MSX = SX(M_num)
    #print('M',M.shape)
    S = SX.sym('S',M.shape[0],M.shape[1])
    QQ, RR = qr(S)
    expR = RR
    QR = Function('QR',[S],[expR])
    R = QR(M)
    R2 = R[n+m:2*(n+m),n+m:2*(n+m)]
    #print('R2.shape',R2.shape)
    R2_inv = inv(R2)
    h_tilde = x_plus1- X_x(x_plus1,p_now)@x_now - X_p(x_plus1,p_now)@p_now
    rho2 = last_V@(last_y - h_tilde)
    x_tilde = h_tilde
    rho3 = -1*(last_W[0:n,0:n]@x_tilde)
    #print(rho3.shape)
    last_component = rho3[0:0+m,:]
    #print(last_component)
    rho2 = vertcat(rho2,last_component)
    P_Lplus1 = R2
    xandpbar = -R2_inv@rho2
    return P_Lplus1, xandpbar

def get_arrival_obj(n,m,P_Lplus1, xandpbar):
    S = MX.sym('S',n)
    pp = MX.sym('pp',m)
    arr_exp = P_Lplus1@(vertcat(S,pp)-xandpbar)
    arr = Function('arr',[vertcat(S,pp)],[arr_exp])
    return arr



def MS_functions_MHE(T, N, x, p, xdot, meas, sigma,sigma2,arr):
    """
    T = length of horizon
    N = Number of shooting intervals in the given horizon
    x = an MX.sym, state of the ODE
    p = an MX.sym, parameters for parameter estimation
    xdot = righthandside of the ODE
    p0 = initial guess of parameters,
    meas = given measurements 
    sigma = variance of meas errors e.g. sigma = [s1,s2,s3,s4]
    sigma2 = variance of states e.g. sigma2 = [s1,s2,s3,s4]
    
    new: arr function (or expression?) of the arrival cost part of the objective taking one input (eg vertcat(S0,p))
    not any more:
    x_opt= optimal end state of the last interval (computed in the last interval)
    p_opt = optimal parameter computed in the last interval
    last_P = P matrix from the last interval
    last_W = W matrix from the last interval 
    
    return: functions F1 and F2 for CNLLS-Problem
    """
    M = meas.shape[1]
    assert N > 0, 'N must be positive'
    assert T > 0, 'T must be positive'
    #assert M%N == 1, 'N must be a divisor of M' 
    
    #sigma2 = sigma
    n = x.shape[0] # dim states
    m = p.shape[0] # dim parameters
    
    y = MX.sym('y',n) # observations are states
    #arr1 = MX.sym('arr1', n)
    #arr2 #todo
    
    # build integrator for whole shooting interval
    dae = {'x':x, 'p':p, 'ode':xdot} 
    opts = {'t0':0, 'tf':T/N}  
    F = integrator('F', 'cvodes', dae, opts)
    
    ss = [1/sigma[i] for i in range(n)]
    d = np.diag(ss)
    S = MX(d)
    
    ss2 = [1/sigma2[i] for i in range(n)]
    dd = np.diag(ss2)
    S2 = MX(dd)
    
    # Pose the NLP
    # Set initial conditions for loop
    S0 = MX.sym('S_0',n)
    v0 = MX.sym('v_0',n)
    Sk = S0
    w1 = []
    w1 = vertcat(w1,Sk) #vertcat(p,Sk)
    
    #Setting objective
    obj = []
    obj_noi = []
    noi_var = []
    arrival_cost = []
    arrival_cost = vertcat(arrival_cost,arr(vertcat(S0,p)))
    # mc denotes the matching constraints
    mc=[]
    
    if True: #(M-1)/N == 1:
        obj = vertcat(obj,S@(meas[:,0] - Sk) ) 
        noi_var = vertcat(noi_var,v0)
        obj_noi = S2@v0
        # Iterate over all shooting intervals
        for k in range(1,N+1): #(1, M):  
            #integrate ode in kth interval and determine s^k
            Fk = F(x0=Sk, p=p) 
            Xk_end = Fk['xf']
            # new symbolic variable for initial condition in interval k
            Sk = MX.sym('S_' + str(k), n)
            
            w1 = vertcat(w1, Sk)
            if k == N: #M-1:
                obj = vertcat(obj,S@(y - Sk) )
                
            else:
                vk = MX.sym('v_'+str(k),n)
                obj = vertcat(obj,S@(meas[:,k] - Sk))
                obj_noi = vertcat(obj_noi,S2@vk)
                noi_var = vertcat(noi_var,vk)
             # Add k-th matching constraint
            mc = vertcat(mc, Xk_end + vk - Sk)
        obj = vertcat(obj,obj_noi)
        obj = vertcat(obj, arrival_cost)
        w1 = vertcat(w1,noi_var)
        w1 = vertcat(w1,p)
    else:
        print('Error: measurements dont match no. of shooting intervals!')
        #return 0
    
    #last_y = meas[:,M-1]
    #P, xandpbar = arrival_cost_values(x_opt, p_opt, last_y, last_P, last_W)
    #S0p = vertcat(S0,p)
    #arrival_cost = vertcat(arrival_cost,P@(S0p-xandpbar))
    #obj = vertcat(arrival_cost, obj)
      
    # convert casadi expressions into functions
    F1 = Function('F1', [w1,y] , [obj] )
    F2 = Function('F2', [w1,y] , [mc] )
    
    #print('w:', w1)
    #print(F1)
    #print(F2)
    return F1,F2

def cnlls_solver(F1, F2, w0, itmax=1, tol=1e-7, ggn = True, show_iteration = False):
    """
    F1: objective function of the CNLLS problem
    F2: constraint function of the CNLLS problem
    w0: initial guess
    itmax: maximal number of iterations
    tol: tolerance
    ggn: bool, True: Generalized Gauß Newton, False: ipopt 
    """
    if ggn:
        x = MX.sym('x', w0.shape[0]) 
        w = w0
        if show_iteration:
            parameter_iterates = []
        
        for i in range(itmax):
            F_1 = F1(w)        
            F_2 = F2(w)
            
            J_exp_1 = jacobian(F1(x), x)   
            J_fun_1 = Function('J_fun_1',[x],[J_exp_1]) 
            J_1 = J_fun_1(w)

            J_exp_2 = jacobian(F2(x), x)
            J_fun_2 = Function('J_fun_2',[x],[J_exp_2]) 
            J_2 = J_fun_2(w)
            
            #Setting qpsolver   
            qp = {'x': x, 'f': dot(F_1 + J_1@x, F_1 + J_1@x) ,'g': F_2 + J_2@x} 
            solver = qpsol('solver', 'qpoases', qp, {"printLevel": "low"}) #tabular, none, low, medium, high, debug_iter
            
            #setting bounds for constraints
            n = F2(w0).shape[0]
            eq_c_bounds = np.zeros(n)
            

            sol = solver(x0=w0, lbg = 0, ubg = 0) #eq_c_bounds, ubg=eq_c_bounds) #did not help
            dw = sol['x']
            #print('This are the dual solutions of the QP wrt g (shape first):', sol['lam_g'].shape[0],sol['lam_g'])
            #print('This are the dual solutions of the QP wrt p (shape first):', sol['lam_p'].shape[0],sol['lam_p'])
            #print('iteration ',i)
            
            if show_iteration:
                parameter_iterates.append(w)
            
            if norm_2(dw) <= tol:
                if show_iteration:
                    return parameter_iterates
                x_opt = w.full().flatten()
                return x_opt #[x_opt]
            
            #update iterate
            w += dw   
            
            if i == (itmax-1):
                #print('itmax reached')
                if show_iteration:
                    return parameter_iterates
                x_opt = w.full().flatten()
                return x_opt #[x_opt]
        
    else:  
        x = MX.sym('x', w0.shape[0]) 
        prob ={'x': x, 'f': dot(F1(x),F1(x)), 'g': F2(x)}
        opts = {}
        opts["ipopt"] = {"max_iter": itmax, "tol": tol,"print_level": 0}
        solver = nlpsol('solver', 'ipopt', prob, opts);
        
        #setting bounds for constraints
        n = F2(w0).shape[0]
        eq_c_bounds = np.zeros(n)
        
        # Solve the NLP
        sol = solver(x0=w0, lbg=eq_c_bounds, ubg=eq_c_bounds)
        x_opt = sol['x'].full().flatten()
    
    return x_opt



##if this does not work, the old version is in v3###
def MHE(P, r0, T, N, length_simulation, x, p, xdot, meas, sigma,sigma2, W, ggn = False):
    '''
    #x0bar: initial guess for arrival cost state
    #p0bar: initial guess for arrival cost parameter
    P: initial guess for arrival covariance matrix
    r0: initial guess for states and p
    
    T: time horizon 
    N: no of shooting intervals in time horizon
    length_simulation: time how long MHE runs
    x: state variable
    p: parameter variable
    xdot: rhs of ODE (must be autonomous and first order)
    meas: measuremnts (for this synthetic scenario)
    sigma = variance of meas errors e.g. sigma = [s1,s2,s3,s4]
    sigma2 = variance of states e.g. sigma2 = [s1,s2,s3,s4]
    W: covariance of x and p for arrival cost computation
    ggn: boolean, indicates whether to use ggn or ipopt
    
    returns: x_opt : sammlung der optimal states over time
             p_opt : sammlung des optimalen ps over time
    
    '''
    n = x.shape[0] # dim state
    m = p.shape[0] # dim parameter
    
    p_opt = []
    x_opt = []
    
    ss2 = [1/sigma2[i] for i in range(n)]
    dd = np.diag(ss2)
    last_V = MX(dd)
    last_P = P
    last_y = meas[:,0]
    last_W = W
    x_o = r0[0:n]
    p_o = r0[r0.shape[0]-m:]
    
    rk = r0
    len_w = rk.shape[0]
    x_opt = vertcat(x_opt,rk[0:n])
    
    dae = {'x': x, 'p': p, 'ode': xdot}
    opts = {'tf': T/N}
    F = integrator('F', 'cvodes', dae, opts)
    
    for k in range(length_simulation):
        toc = time.perf_counter()
        print('starting (time horizon) loop no. ', k+1)
        ## prep phase
        # update arrival cost (or in ms_fun)
        #later... see below in ms_fun
        #print(rk)
        cut = rk.shape[0]-m
        xk = rk[cut-((N+1)*n):cut-(N*n)] #rk[N*n:(N+1)*n]
        pk = rk[cut:]
        #print(xk.shape)
        #print(pk.shape)
        
        # integrate on next new shooting interval
        Res = F(x0=xk,p=pk)
        xkminus = Res['xf']
        
        # set rkminus and dkminus
        rkminus1 = rk[n:(N+1)*n]
        rkminus2 = xkminus
        rkminus3 = rk[(N+2)*n:(2*N+1)*n]
        rkminus4 = np.zeros(n)
        rkminus5 = pk
        rkminus = vertcat(vertcat(vertcat(vertcat(rkminus1,rkminus2),rkminus3),rkminus4),rkminus5)
        #print('rk- shape: ',rkminus.shape)
        #Dkminus_y = meas[:,] will be measured later, not available yet
        
        # setup ms setting with variable y
        last_W = W
        last_V = last_V
        x_o = rk[0:n]
        p_o = pk
        #print('k',k)
        #print('k+N+1',k+N+1)
        #print(meas[:,k:k+N+1].shape)
        #,k+N:k+2*N-1
        
        P_Lplus1, xandpbar = arrival_cost_values(x,p,xdot,T,N,x_o, p_o, last_y, last_P,last_V, last_W)
        arr = get_arrival_obj(n,m,P_Lplus1, xandpbar)
        FF1,FF2 = MS_functions_MHE(T, N, x, p, xdot, meas[:,k:k+N+1], sigma,sigma2, arr) 
        last_y = meas[:,k+1]
        last_P =P_Lplus1
        #compute vector components of QP (dep on variable y)
        # later..
        #compute matrix components of QP
        # later..
        
        ### we get the new measuremnet y
        ## est phase
        # set variable y in QP
        w = MX.sym('w', len_w)
        y = meas[:,k+N+1]
        F1 = Function('F1',[w],[FF1(w,y)])
        F2 = Function('F2',[w],[FF2(w,y)])
        
        # solve qp with initial value rkminus
        w0 = rkminus
        drk = cnlls_solver(F1, F2, w0, itmax=1, tol=1e-7, ggn=ggn, show_iteration = False)
        
        
        #set rk, xk, pk
        rk = drk #rkminus + drk -> here was an mistake! I updated twice!
        #print('rk new shape',rk.shape)
        xk = rk[N*n:(N+1)*n]
        #print('xk new shape',xk.shape)
        pk = rk[(2*N+1)*n:]
        #print('pk new',pk)
        p_opt = vertcat(p_opt,pk)
        if k == 0:
            x_opt = vertcat(x_opt,rk[0:(N+1)*n])
        else:
            x_opt = vertcat(x_opt,xk)
        tic = time.perf_counter()
        if k <=10:
            print('needed time for loop '+str(k),tic-toc)
    return x_opt, p_opt
