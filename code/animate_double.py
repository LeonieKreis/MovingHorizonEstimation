# import ipywidgets as widgets
from casadi import *
import numpy as np
import matplotlib.pyplot as plt

# from matplotlib import rc #,animation
import matplotlib.animation as animation
from collections import deque
# from IPython.display import HTML

#parameters
g = 9.81
l1 = 0.323
l2 = 0.480
a1 = 0.2145
a2 = 0.223
m1 = 0.853
m2 = 0.510
J1 = 0.0126
J2 = 0.0185
#d1 = 0.005
#d2 = 0.005


def get_cart_x(Res):
    return Res[0]
def get_cart_y(Res):
    return 0
def get_ball1_x(Res):
    # return Res[0]-a1*np.sin(Res[2])  # this is the center of mass, not the ball
    return Res[0]-l1*np.sin(Res[2])  # this should be the ball
def get_ball1_y(Res):
    # return a1*np.cos(Res[2])  # this is the center of mass, not the ball
    return l1*np.cos(Res[2])
def get_ball2_x(Res):
    # return Res[0]-l1*np.sin(Res[2])-a2*np.sin(Res[4])  # this is the center of mass, not the ball
    return Res[0]-l1*np.sin(Res[2])-l2*np.sin(Res[4]) 
def get_ball2_y(Res):
    # return l1*np.cos(Res[2])+a2*np.cos(Res[4])  # this is the center of mass, not the ball
    return l1*np.cos(Res[2])+l2*np.cos(Res[4])


def animate_double_pendulum(gifname, ss,title= 'Double Pendulum', T=-1, N=-1):
# T=-1, N=-1 are placeholders in case we don't want to presribe T and N 
    
    s = np.array(ss).flatten()
    ll = int(s.shape[0]/6)
    #t_stop = 5  # how many seconds to simulate
    history_len = ll  # how many trajectory points to display

    # create a time array from 0..t_stop sampled at 0.02 second steps
    if T>0 and N>0:
        dt = T/N
    else:
        dt = 0.1  # default
    #t = np.arange(0, t_stop, dt)

    fig = plt.figure(figsize=(10, 8))
    fig.suptitle(title, fontsize=16)
    ax = fig.add_subplot(autoscale_on=False, xlim=(-1.5, 1.5), ylim=(-1, 1))
    ax.set_aspect('equal')
    ax.grid()

    line1, = ax.plot([], [], 'o-', lw=2)
    dot1, = ax.plot([],[],'bs', lw=3)
    line2, = ax.plot([], [], 'o-', lw=2)
    #dot2, = ax.plot([],[],'bs', lw=3)
    trace, = ax.plot([], [], '.-', lw=1, ms=2)

    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)

    x1 = np.array([s[6 * i] for i in range(ll)])
    y1 = np.zeros(ll)
    x2 = np.array([get_ball1_x(s[6*i:6*(i+1)]) for i in range(ll)])
    y2 = np.array([get_ball1_y(s[6*i:6*(i+1)]) for i in range(ll)])
    x3 = np.array([get_ball2_x(s[6*i:6*(i+1)]) for i in range(ll)])
    y3 = np.array([get_ball2_y(s[6*i:6*(i+1)]) for i in range(ll)])
    
    def animate(i):
        thisx1 = [x1[i], x2[i]]  
        thisy1 = [y1[i], y2[i]] 
        
        thisx2 = [x2[i], x3[i]]  
        thisy2 = [y2[i], y3[i]]

        if i == 0:
            history_x.clear()
            history_y.clear()

        history_x.appendleft(thisx2[1])
        history_y.appendleft(thisy2[1])

        line1.set_data(thisx1, thisy1)
        line2.set_data(thisx2, thisy2)
        trace.set_data(history_x, history_y)
        dot1.set_data(thisx1[0],thisy1[0])
        #dot2.set_data(thisx2[0],thisy2[0])
        time_text.set_text(time_template % (i*dt))
        return line1, line2, trace, dot1, time_text

    no_anim = ll
    #print(no_anim)

    if T>0 and N>0:
        interval = int(N/(T*10))
        fps = int(N*10/T)
    else:
        interval = 500  # your previous default
        fps = 6

    myAnimation = animation.FuncAnimation(
        fig, animate, no_anim, interval=interval, blit=True)
    myAnimation.save(gifname, writer="Pillow", fps=fps)
    plt.show()
    
    
def plot_widget(simulations): #, n):
    plt.figure(figsize=(10,8))
    plt.title("Comparison of the simulations:")
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.,1.)
    
    
    if 'true' in simulations:
        ss = Res1
        s = np.array(ss).flatten()
        ll = int(s.shape[0]/6)
        x1 = np.array([s[6 * i] for i in range(ll)])
        y1 = np.zeros(ll)
        x2 = np.array([get_ball1_x(s[6*i:6*(i+1)]) for i in range(ll)])
        y2 = np.array([get_ball1_y(s[6*i:6*(i+1)]) for i in range(ll)])
        x3 = np.array([get_ball2_x(s[6*i:6*(i+1)]) for i in range(ll)])
        y3 = np.array([get_ball2_y(s[6*i:6*(i+1)]) for i in range(ll)])
        plt.plot(x1,y1,marker = 'bs')#, label = "explicit")
        plt.plot(x2,y2,marker = 'o-')#, label = "explicit")
        plt.plot(x3,y3,marker = 'o')#, label = "explicit")
    if 'measured' in simulations:
        ss = meas1
        s = np.array(ss).flatten()
        ll = int(s.shape[0]/6)
        x1 = np.array([s[6 * i] for i in range(ll)])
        y1 = np.zeros(ll)
        x2 = np.array([get_ball1_x(s[6*i:6*(i+1)]) for i in range(ll)])
        y2 = np.array([get_ball1_y(s[6*i:6*(i+1)]) for i in range(ll)])
        x3 = np.array([get_ball2_x(s[6*i:6*(i+1)]) for i in range(ll)])
        y3 = np.array([get_ball2_y(s[6*i:6*(i+1)]) for i in range(ll)])
        plt.plot(x1,y1,marker = 'bs')#, label = "explicit")
        plt.plot(x2,y2,marker = 'o-')#, label = "explicit")
        plt.plot(x3,y3,marker = 'o')#, label = "explicit")
    if 'open-loop' in simulations:
        ss = s_opt
        s = np.array(ss).flatten()
        ll = int(s.shape[0]/6)
        x1 = np.array([s[6 * i] for i in range(ll)])
        y1 = np.zeros(ll)
        x2 = np.array([get_ball1_x(s[6*i:6*(i+1)]) for i in range(ll)])
        y2 = np.array([get_ball1_y(s[6*i:6*(i+1)]) for i in range(ll)])
        x3 = np.array([get_ball2_x(s[6*i:6*(i+1)]) for i in range(ll)])
        y3 = np.array([get_ball2_y(s[6*i:6*(i+1)]) for i in range(ll)])
        plt.plot(x1,y1,marker = 'bs')#, label = "explicit")
        plt.plot(x2,y2,marker = 'o-')#, label = "explicit")
        plt.plot(x3,y3,marker = 'o')#, label = "explicit")
        
    if 'mhe' in simulations:
        ss = xx
        s = np.array(ss).flatten()
        ll = int(s.shape[0]/6)
        x1 = np.array([s[6 * i] for i in range(ll)])
        y1 = np.zeros(ll)
        x2 = np.array([get_ball1_x(s[6*i:6*(i+1)]) for i in range(ll)])
        y2 = np.array([get_ball1_y(s[6*i:6*(i+1)]) for i in range(ll)])
        x3 = np.array([get_ball2_x(s[6*i:6*(i+1)]) for i in range(ll)])
        y3 = np.array([get_ball2_y(s[6*i:6*(i+1)]) for i in range(ll)])
        plt.plot(x1,y1,marker = 'bs')#, label = "explicit")
        plt.plot(x2,y2,marker = 'o-')#, label = "explicit")
        plt.plot(x3,y3,marker = 'o')#, label = "explicit")
        
    #plt.plot(0,1,marker = 'o',label = "starting point")
    plt.legend()
    plt.show()
    
    
    
    
def plot_widget_old(integrators, n):
    T = 2*np.pi
    p0 = 0
    q0 = 1
    del_t = T/n
    plt.figure(figsize=(7,7))
    plt.xlabel("p")
    plt.ylabel("q")
    plt.title("Comparison of the integrators:")
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.5,1.5)
    # exact solution
    N = 100
    p_ex = np.empty(N+1)
    q_ex = np.empty(N+1)
    p_ex[0] = p0
    q_ex[0] = q0
    for i in range(1,N+1):
        p_ex[i] = q0*np.cos(i*del_t + np.pi/2)
        q_ex[i] = q0*np.sin(i*del_t + np.pi/2)  
    plt.plot(p_ex,q_ex, label = "exact")
    if 'explicit' in integrators:
        p_e = np.empty(n+1)
        q_e = np.empty(n+1)
        p_e[0] = p0
        q_e[0] = q0
        for i in range(1,n+1):
            p_e[i] = p_e[i-1]- del_t*q_e[i-1]
            q_e[i] = q_e[i-1]+del_t *p_e[i-1]
        plt.plot(p_e,q_e,marker = 'o', label = "explicit")
    if 'implicit' in integrators:
        p_i = np.empty(n+1)
        q_i = np.empty(n+1)
        p_i[0] = p0
        q_i[0] = q0
        for i in range(1,n+1):
            q_i[i] = 1/(1+del_t*del_t)*(q_i[i-1]+del_t *p_i[i-1])
            p_i[i] = p_i[i-1]- del_t*q_i[i]
        plt.plot(p_i,q_i,marker = 'o', label = "implicit")
    if 'euler-B' in integrators:
        p_B = np.empty(n+1)
        q_B = np.empty(n+1)
        p_B[0] = p0
        q_B[0] = q0
        for i in range(1,n+1):
            p_B[i] = p_B[i-1]- del_t*q_B[i-1]
            q_B[i] = q_B[i-1]+ del_t*p_B[i]  
        plt.plot(p_B,q_B,marker = 'o', label = "Euler-B")
    plt.plot(0,1,marker = 'o',label = "starting point")
    plt.legend()
    plt.show()