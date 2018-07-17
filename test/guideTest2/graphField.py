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
                if (len(text) == 3):
                    x.append(float( text[0]) )
                    y.append(float( text[1]) )
                    z.append(float( text[2]) )
                else:
                    bx.append(float( text[0]) )
                    by.append(float( text[4]) )
                    bz.append(float( text[8]) )
    except IOError:
        print("Error reading BFCut.out")
        return

    try:
        with open ('segSolBField.txt',"r") as f1:
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
        print("Error reading comsol text file")
        return

    if (len(x) != len(bx)):
        print("Error: not an equal number of xyz coords and b field values")
        return

    fig2 = plt.figure(2)
    plt.quiver(x[::5], z[::5], bx[::5], bz[::5], units='x')
    plt.title('Vector Field xz plane')

    fig3 = plt.figure(3)
    plt.quiver(x0, z0, bx0, bz0, units='x')
    plt.title('Vector Field xz plane (Original COMSOL)')




    plt.show()

    return

if ( __name__ == '__main__' ):
    main()
