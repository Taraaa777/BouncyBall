###########
# Warning #
# Code will run forever if variable != 'h' or variable != 'd'
###########
# Program plots initial parabola of ball dropped onto sphere with initial coordinates variable.
# Assumptions include elastic collisions, no wind res etc.
import numpy as np
import matplotlib.pyplot as plt

# Calculates a corresponding set of coords, given set of x- or y-coords
def circlise(x):
    y = []
    for xi in x:
        yi = np.sqrt(R**2 - xi**2)
        y.append(yi)
    return y

# Calculate maximum root of a quadratic equation
def root(a,b,c):
    det = np.sqrt(b*b - 4*a*c)
    #x = max( (-b + det)/(2*a), (-b - det)/(2*a))
    x = (-b - det)/(2*a)
    #if x < 0:
    #    print('max root is negative')
    return x

# ends list of coordinates when they enter sphere
def phasecheck(x,y):
    for i in range(len(x)):
        if x[i]**2+y[i]**2<R**2 and i != 0:
            x = x[0:i]
            y = y[0:i]
            return x, y
    return x, y

# Initialises diagram with sphere, labels etc.
def initialPlot(R,variable):
    sphere_x = np.linspace(0,R,250)
    sphere_y = circlise(sphere_x)
    
    plt.figure()
    # 11 too small for \columnwidth diagrams
    plt.title('Path of particle dependant on {0}'.format(variable),fontsize=20)
    plt.xlabel('horizontal displacement',fontsize=20)
    plt.ylabel('vertical displacement',fontsize=20)
    plt.plot(0,0,'.',color='k',label='centre')
    plt.plot(sphere_x,sphere_y,color='k',label = 'surface')
    plt.legend('upper right')

# Computes initial velocities
def params(g,R,H,d0):
    theta = np.arcsin(d0/R)
    v0 = np.sqrt(2*g*H)
    vx = v0 * np.sin(2*theta)
    vy = v0 * np.cos(2*theta)
    #if np.cos(2*theta) < 0:
     #   print('negative v0')
    if 2*g*H >= 1.01*(vx**2 + vy**2):
        print('conservation of energy not maintained:')
        print(vx**2 + vy**2 - 2*g*H)
    return vx,vy

# function
def func(t,d,h,vx,vy,g,R):
    return (d+vx*t)**2 + (h+vy*t-0.5*g*t**2)**2 -R**2

# function derivative
def derivFunc(t,d,h,vx,vy,g,R):
    return 2*(d+vx*t)*vx + 2*(h+vy*t-0.5*g*t**2)*(+vy-g*t)

# Calculates the root via Newton-Raphson
def newtonRaphson(t,d,h,vx,vy,g,R):
    h = func(t,d,h,vx,vy,g,R) / derivFunc(t,d,h,vx,vy,g,R)
    while abs(h) >= 0.001:
        h = func(t,d,h,vx,vy,g,R)/derivFunc(t,d,h,vx,vy,g,R)
        print(h)
        # x(i+1) = x(i) - f(x) / f'(x)
        t = t - h
    return t

# Plots parabola of ball
def plotPath(g,R,d0,H,variable,i):
    # Important points:
    initial = [d0,circlise([d0])[0]+H]
    impact = [initial[0], initial[1]-H]
    # Drop:
    drop = np.linspace(initial[1],impact[1],100)
    x_drop = [initial[0] for ele in drop]

    #Resolution of velocities
    vx, vy = params(g,R,H,d0)
    tmax = root(-0.5*g,vy,impact[1]) # root for s = ut - 0.5gt^2, a = -0.5g, b = vy, c = -s
    #print((2 * vy) / g - tmax) #check if tmax for zero displacement agrees
    t = tmax*np.linspace(0,1,500)#[:,None] #create timespan

    #calculation for displacement per timestep
    x = vx*t
    y = vy*t - 0.5*g*t**2
    delta_d = impact[0]+x
    delta_h = impact[1]+y

    #cut list directly before ball enters sphere
    delta_d,delta_h = phasecheck(delta_d,delta_h)

    #Graphing
    if variable.lower() == 'h':
        plt.plot(delta_d,delta_h,color='C{0}'.format(i+1),label = 'h = {0}'.format(H)) #plot ball path
        plt.plot(initial[0],initial[1],'x',color='C{0}'.format(i+1)) # Initial point
        plt.legend()
    else:
        plt.plot(delta_d,delta_h,color='C{0}'.format(i+1),label = 'd = {0}'.format(d0)) #plot ball path
        plt.plot(initial[0],initial[1],'x',color='C{0}'.format(i+1)) # Initial point
        plt.plot(x_drop,drop,'--',color='C{0}'.format(i+1)) # drop
        plt.legend()

# Newton Raphson attempt to plot trajectory
# Getting stuck in N-R, infinite loop - could be incorrect function, derivative, or function does not have roots.
def modelBounces(g,R,d0,H,variable,i):
    # Important points:
    initial = [d0,circlise([d0])[0]+H]
    impact = [initial[0], initial[1]-H]
    # Drop:
    drop = np.linspace(initial[1],impact[1],100)
    x_drop = [initial[0] for ele in drop]

    # Breaks between here and 'print statement
    #Resolution of velocities
    vx, vy = params(g,R,H,d0)
    t0 = 2*vx/g # In theory, this is near the crossing point as we have net zero vert. disp.
    print('gets to here')
    tmax = newtonRaphson(t0,impact[0],impact[1],vx,vy,g,R)
    t = tmax*np.linspace(0,1,500)#[:,None] #create timespan
    print('gets to here')
    #calculation for displacement per timestep
    x = vx*t
    y = vy*t - 0.5*g*t**2
    delta_d = impact[0]+x
    delta_h = impact[1]+y

    # Graphing
    if variable.lower() == 'h':
        plt.plot(delta_d,delta_h,color='C{0}'.format(i+1),label = 'h = {0}'.format(H)) #plot ball path
        plt.plot(initial[0],initial[1],'x',color='C{0}'.format(i+1)) # Initial point
        plt.legend()
    else:
        plt.plot(delta_d,delta_h,color='C{0}'.format(i+1),label = 'd = {0}'.format(d0)) #plot ball path
        plt.plot(initial[0],initial[1],'x',color='C{0}'.format(i+1)) # Initial point
        plt.plot(x_drop,drop,'--',color='C{0}'.format(i+1)) # drop
        plt.legend()
    onsphere = False
    if delta_d[-1]**2 + delta_h[-1]**2 <= R**2:
        onsphere = True
    return delta_d, delta_h, onsphere

#Initial parameters
g = 9.81
R = 1
d0 = 0.2
H = 0.6
variable = 'd'
if variable.lower() == 'h':
    initialPlot(R,'height')    
    Hrange = np.linspace(0,0.6,7)
    for i,H in enumerate(Hrange):
        plotPath(g,R,d0,round(H,2),variable,i)
elif variable.lower() == 'd':
    initialPlot(R,'distance')
    #drange = np.linspace(0,0.25,6)
    drange = (0,0.1,0.15,0.2,0.25) #nicer diagram, omits 0.05 for clarity
    for i,d in enumerate(drange):
        plotPath(g,R,round(d,2),H,variable,i)

else:
    initialPlot(R,'N-R')    
    drange = [0.2]#np.linspace(0.1,0.2,2)
    for i,d in enumerate(drange):
        onsphere = True
        while onsphere:
            modelBounces(g,R,round(d,2),H,variable,i)
            #delta_d,delta_h,onsphere = modelBounces(g,R,round(d,2),H,variable,i)
            
plt.axis('square')
plt.savefig('variable{0}.pdf'.format(variable))
plt.show()  