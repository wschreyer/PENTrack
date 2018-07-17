#!/usr/bin/env python
# Graphs projection of spin vector on magnetic field as a function of time
# Requires neutronspin.out file

def main():
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    t = []
    x = []
    y = []
    z = []
    sx = []
    sy = []
    sz = []
    wx = []
    wy = []
    wz = []
    bx = []
    by = []
    bz = []
    proj = []

    try:
        with open ('000000000000neutronspin.out',"r") as f1:
            lines = f1.readlines()[8:]
            for num, line in enumerate(lines):
                text = line.split()
                t.append(float( text[2]) )
                x.append(float( text[3]) )
                y.append(float( text[4]) )
                z.append(float( text[5]) )
                sx.append(float( text[6]) )
                sy.append(float( text[7]) )
                sz.append(float( text[8]) )
                wx.append(float( text[9]) )
                wy.append(float( text[10]) )
                wz.append(float( text[11]) )
                bx.append(float( text[12]) )
                by.append(float( text[13]) )
                bz.append(float( text[14]) )
    except IOError:
        print("Error reading neutronspin.out")
        return

    # Calculate spin projection on magnetic fields
    for bxi, byi, bzi, sxi, syi, szi in zip (bx, by, bz, sx, sy, sz):
        if (norm(bxi, byi, bzi) == 0):
            proj.append(5)    # Some random number > 1
            continue
        proj.append( (sxi*bxi + syi*byi + szi*bzi)/norm(bxi, byi, bzi) )

    fig1 = plt.figure(1)
    plt.plot(x, proj)
    plt.grid(True)
    plt.title('x vs Spin Proj')
    plt.xlabel('x [m]')
    plt.ylabel('Projection of S on B')

    fig2 = plt.figure(2)
    plt.plot(t, proj)
    plt.grid(True)
    plt.title('t vs Spin Proj')
    plt.xlabel('t [s]')
    plt.ylabel('Projection of S on B')

    plt.show()

    return

def norm (x, y, z):
    #Finds norm of cartesian vector components
    import numpy as np
    return np.sqrt(x*x + y*y + z*z)

if ( __name__ == '__main__' ):
    main()
