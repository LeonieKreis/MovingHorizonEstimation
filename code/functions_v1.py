from casadi import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import random

def generate_meas(x, p, xdot, T, sigma, p_true, M, x0):
    """
    x = an MX.sym, state of the ODE
    p = an MX.sym, parameters of the ODE
    xdot = an MX.sym, righthandside of the ODE
    T = Final Time 
    sigma = variance of the parameters 
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
    noise= np.random.normal(0.0, sigma, (n,M+1))
    meas = meas_exact + noise
    return(meas)

def arrival_cost_values(x_opt, p_opt, last_y, last_P, last_W):
    print('Arrival cost not implemented yet!')
    P = 0
    xandpbar = 0
    return P, xandpbar

def MS_functions_MHE(T, N, x, p, xdot, meas, sigma, x_opt, p_opt, last_P, last_W):
    """
    T = length of horizon
    N = Number of shooting intervals in the given horizon
    x = an MX.sym, state of the ODE
    p = an MX.sym, parameters for parameter estimation
    xdot = righthandside of the ODE
    p0 = initial guess of parameters,
    meas = given measurements 
    sigma = variance of parameters
    
    new:
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
    
    n = x.shape[0] # dim states
    m = p.shape[0] # dim parameters
    
    y = MX.sym('y',n) # observations are states
    #arr1 = MX.sym('arr1', n)
    #arr2 #todo
    
    # build integrator for whole shooting interval
    dae = {'x':x, 'p':p, 'ode':xdot} 
    opts = {'t0':0, 'tf':T/N}  
    F = integrator('F', 'cvodes', dae, opts)
    
    # Pose the NLP
    # Set initial conditions for loop
    S0 = MX.sym('S_0',n)
    Sk = S0
    w1 = []
    w1 = vertcat(w1,Sk) #vertcat(p,Sk)
    
    #Setting objective
    obj = []
    arrival_cost = []
    # mc denotes the matching constraints
    mc=[]
    
    if True: #(M-1)/N == 1:
        obj = vertcat(obj,(meas[:,0] - Sk) / sigma) 
        # Iterate over all shooting intervals
        for k in range(1,N+1): #(1, M):  
            #integrate ode in kth interval and determine s^k
            Fk = F(x0=Sk, p=p) 
            Xk_end = Fk['xf']
            # new symbolic variable for initial condition in interval k
            Sk = MX.sym('S_' + str(k), n)
            w1 = vertcat(w1, Sk)
            if k == N: #M-1:
                obj = vertcat(obj,(y - Sk) / sigma)
            else:
                obj = vertcat(obj,(meas[:,k] - Sk) / sigma)
             # Add k-th matching constraint
            mc = vertcat(mc, Xk_end - Sk)
        obj = vertcat(obj, )
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

def MHE(x0bar, p0bar, P, r0, T, N, length_simulation, x, p, xdot, meas, sigma, W, ggn = False):
    '''
    x0bar: initial guess for arrival cost state
    p0bar: initial guess for arrival cost parameter
    P: initial guess for arrival covariance matrix
    r0: initial guess for states and p
    
    T: time horizon 
    N: no of shooting intervals in time horizon
    length_simulation: time how long MHE runs
    x: state variable
    p: parameter variable
    xdot: rhs of ODE (must be autonomous and first order)
    meas: measuremnts (for this synthetic scenario)
    sigma: variance of iid measurement noise
    W: covariance of x and p for arrival cost computation
    
    returns: x_opt : sammlung der optimal states over time
             p_opt : sammlung des optimalen ps over time
    
    '''
    n = x.shape[0] # dim state
    m = p.shape[0] # dim parameter
    
    p_opt = []
    x_opt = []
    
    rk = r0
    len_w = rk.shape[0]
    x_opt = vertcat(x_opt,rk[0:n])
    
    dae = {'x': x, 'p': p, 'ode': xdot}
    opts = {'tf': T/N}
    F = integrator('F', 'cvodes', dae, opts)
    
    for k in range(length_simulation):
        print('starting (time horizon) loop no. ', k+1)
        ## prep phase
        # update arrival cost (or in ms_fun)
        #later... see below in ms_fun
        #print(rk)
        cut = rk.shape[0]-m
        xk = rk[cut-n:cut] #rk[N*n:(N+1)*n]
        pk = rk[cut:]
        #print(xk.shape)
        #print(pk.shape)
        
        # integrate on next new shooting interval
        Res = F(x0=xk,p=pk)
        xkminus = Res['xf']
        
        # set rkminus and dkminus
        rkminus1 = rk[n:(N+1)*n]
        rkminus2 = xkminus
        rkminus3 = pk
        rkminus = vertcat(vertcat(rkminus1,rkminus2),rkminus3)
        #print('rk- shape: ',rkminus.shape)
        #Dkminus_y = meas[:,] will be measured later, not available yet
        
        # setup ms setting with variable y
        last_P =0
        last_W = 0
        #print('k',k)
        #print('k+N+1',k+N+1)
        #print(meas[:,k:k+N+1].shape)
        #,k+N:k+2*N-1
        FF1,FF2 = MS_functions_MHE(T, N, x, p, xdot, meas[:,k:k+N+1], sigma, rk, rk, last_P, last_W) # not debugged yet!
        
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
        pk = rk[(N+1)*n:]
        #print('pk new',pk)
        p_opt = vertcat(p_opt,pk)
        if k == 0:
            x_opt = vertcat(x_opt,rk[0:(N+1)*n])
        else:
            x_opt = vertcat(x_opt,xk)
    
    return x_opt, p_opt
