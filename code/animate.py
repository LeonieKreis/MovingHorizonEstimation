import ipywidgets as widgets
from casadi import *
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc #,animation
import matplotlib.animation as animation
from collections import deque
from IPython.display import HTML

l = 1

def get_cart_x(Res):
    return Res[0]
def get_cart_y(Res):
    return 0
def get_ball_x(Res):
    return Res[0]-l*np.sin(Res[2])
def get_ball_y(Res):
    return -l*np.cos(Res[2])


def animate_pendulum(gifname, ss):
    
    s = np.array(ss).flatten()
    ll = int(s.shape[0]/4)
    #t_stop = 5  # how many seconds to simulate
    history_len = ll  # how many trajectory points to display

    # create a time array from 0..t_stop sampled at 0.02 second steps
    dt = 0.05
    #t = np.arange(0, t_stop, dt)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-1, 1), ylim=(-1.3, 1.3))
    ax.set_aspect('equal')
    ax.grid()

    line, = ax.plot([], [], 'o-', lw=2)
    trace, = ax.plot([], [], '.-', lw=1, ms=2)
    dot, = ax.plot([],[],'bs', lw=3)

    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)

    x1 = np.array([s[4 * i] for i in range(ll)])
    y1 = np.zeros(ll)
    x2 = np.array([get_ball_x(s[4*i:4*(i+1)]) for i in range(ll)])
    y2 = np.array([get_ball_y(s[4*i:4*(i+1)]) for i in range(ll)])

    def animate(i):
        thisx = [x1[i], x2[i]]  
        thisy = [y1[i], y2[i]] 

        if i == 0:
            history_x.clear()
            history_y.clear()

        history_x.appendleft(thisx[1])
        history_y.appendleft(thisy[1])

        line.set_data(thisx, thisy)
        trace.set_data(history_x, history_y)
        dot.set_data(thisx[0],thisy[0])
        time_text.set_text(time_template % (i*dt))
        return line, trace, dot, time_text

    no_anim = ll
    #print(no_anim)

    myAnimation = animation.FuncAnimation(
        fig, animate, no_anim, interval=500, blit=True)
    myAnimation.save(gifname, writer="Pillow", fps=6)
    plt.show()
    