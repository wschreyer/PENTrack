#!/usr/bin/env python
def main():
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    x = []
    y = []
    z = []
    bx = []
    by = []
    bz = []
    x0 = []
    y0 = []
    z0 = []
    bx0 = []
    by0 = []
    bz0 = []
    sx = []         # Spin vectors
    sy = []
    sz = []
    tSpin = []      # time interval (corresponding to spin vectors)

    # Read in files
    try:
        with open ('BFCut.out',"r") as f1:
            lines = f1.readlines()[1:]
            for num, line in enumerate(lines):
                text = line.split(' ')
                x.append(float( text[0]) )
                y.append(float( text[1]) )
                z.append(float( text[2]) )
                bx.append(float( text[3]) )
                by.append(float( text[7]) )
                bz.append(float( text[11]) )
    except IOError:
        print("Error reading BFCut.out")
        return

    try:
        with open ('comsolField.txt',"r") as f1:
            lines = f1.readlines()[8:]
            for num, line in enumerate(lines):
                text = line.split()
                x0.append(float( text[0]) )
                y0.append(float( text[1]) )
                z0.append(float( text[2]) )
                bx0.append(float( text[3]) )
                by0.append(float( text[4]) )
                bz0.append(float( text[5]) )
    except IOError:
        print("Error reading comsolField.txt")
        return

    try:
        with open ('000000000000neutronspin.out',"r") as f1:
            lines = f1.readlines()[1:]
            for num, line in enumerate(lines):
                tSpin.append(float( line.split(' ')[2]) )
                sx.append(float( line.split(' ')[6]) )
                sy.append(float( line.split(' ')[7]) )
                sz.append(float( line.split(' ')[8]) )
    except IOError:
        print("Error reading neutron spin file")
        return

    if (len(x) != len(bx)):
        print("Error: not an equal number of xyz coords and b field values")
        return

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111, projection='3d')
    plt.title('Neutron spin')
    ax1.scatter(sx, sy, sz)
    ax1.set_xlim3d(-1,1)
    ax1.set_ylim3d(-1,1)
    ax1.set_zlim3d(-1,1)
    ax1.set_xlabel('P(x)')
    ax1.set_ylabel('P(y)')
    ax1.set_zlabel('P(z)')
    ax1.view_init(30,220)

    fig2 = plt.figure(2)
    plt.quiver(z[::5], y[::5], bz[::5], by[::5], units='x')
    plt.title('Vector Field zy plane')

    fig3 = plt.figure(3)
    plt.quiver(z0, y0, bz0, by0, units='x')
    plt.title('Vector Field zy plane (Original COMSOL)')




    plt.show()

    return

if ( __name__ == '__main__' ):
    main()
