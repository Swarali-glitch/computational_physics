import numpy as np
import matplotlib.pyplot as plt

m = 0.01
k = 0.0001
g = 10

steps = 10  # You can change the number of steps you want to divide interval into

v = np.zeros(steps, dtype = float)

def DE(t,v):
  return (g-(k*pow(v,2)/m))


def RungeKutta(t_0, v_0, t, h):
  n = int((t-t_0)/h)
  v = v_0
  for i in range(1,n+1):
    s1 = h * DE(t_0, v)
    s2 = h * DE(t_0+(h/2),v +(s1/2))
    s3 = h * DE(t_0+(h/2),v +(s2/2))
    s4 = h * DE(t_0+h, v+s3)

    v = v + (1/6)*(s1+ (2*s2) + (2*s3) + s4)
    t_0 = t_0 + h
  return v

t_0 = 0
v_0 = 0
t = np.linspace(0,10,num=steps)

for i in range(0,steps):
  v[i] = RungeKutta(t_0, v_0, t[i], 0.2)

print(t)
print(v)

fig,ax = plt.subplots()
ax.plot(t,v)
ax.set(xlabel='time (s)', ylabel='velocity (m/s)',
       title='velocity vs. time for particle falling under gravity acted by air drag')
ax.grid()
plt.show()
