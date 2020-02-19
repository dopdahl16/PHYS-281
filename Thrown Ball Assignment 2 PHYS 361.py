import math
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

def newline(p1, p2):
    ax = plt.gca()
    xmin, xmax = ax.get_xbound()

    if(p2[0] == p1[0]):
        xmin = xmax = p1[0]
        ymin, ymax = ax.get_ybound()
    else:
        ymax = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmax-p1[0])
        ymin = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmin-p1[0])

    l = mlines.Line2D([xmin,xmax], [ymin,ymax])
    ax.add_line(l)
    return l

m = 1
k = 1
h = 0.001
x0 = 0
y0 = 0
xd0 = 86.6
yd0 = -50
x = [x0]
y = [y0]
xd = [xd0]
yd = [yd0]

for i in range(1040):
    xd.append(xd[i] - (h*m*k)*math.sqrt(xd[i]**2 + yd[i]**2)*xd[i])
    yd.append(yd[i] - (h*m*k)*math.sqrt(xd[i]**2 + yd[i]**2)*yd[i] + 9.81*h)
    x.append(x[i] + h*xd[i])
    y.append(y[i] + h*yd[i])
        
for i in range(len(y)):
    y[i] = -1*y[i]

plt.scatter(x,y,s=1,marker='o')

p1 = [0,0]
p2 = [7,0]

newline(p1,p2)

plt.title("Thrown ball w/ air resistance")
plt.xlabel("Distance (X)")
plt.ylabel("Height (Y)")

plt.show()