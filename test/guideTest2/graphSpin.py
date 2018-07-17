#!/usr/bin/env python
# Makes some quick and dirty plots

def main():
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    nbins = 200        # number of histogram bins for polarization plots
    Bx = []         # Magnetic field seen by neutron
    By = []
    Bz = []
    tBField = []    # time interval (corresponding to B field vectors)
    Bnorm = []
    x = []          # Position of neutron
    y = []
    z = []
    sxStart = []    # starting and ending values for a given neutron
    syStart = []
    szStart = []
    bxStart = []
    byStart = []
    bzStart = []
    sxEnd = []
    syEnd = []
    szEnd = []
    bxEnd = []
    byEnd = []
    bzEnd = []

    neutronCounter = 0

    try:
        with open ('000000000000neutrontrack.out',"r") as f2:
            lines = f2.readlines()[1:]
            for num, line in enumerate(lines):
            # For each neutron, append B field start/end values
            # Kinda ugly method of doing this but whatever
                Bx.append(float( line.split(' ')[12]) )
                By.append(float( line.split(' ')[16]) )
                Bz.append(float( line.split(' ')[20]) )
                tBField.append(float( line.split(' ')[3]) )
                x.append(float( line.split(' ')[4]) )
                y.append(float( line.split(' ')[5]) )
                z.append(float( line.split(' ')[6]) )
                Bnorm.append( norm(Bx[num], By[num], Bz[num]))
                if neutronCounter < int( line.split(' ')[1]):
                    bxStart.append(Bx[num])
                    byStart.append(By[num])
                    bzStart.append(Bz[num])
                    if neutronCounter != 0:
                        bxEnd.append(Bx[num-1])
                        byEnd.append(By[num-1])
                        bzEnd.append(Bz[num-1])
                    neutronCounter += 1
            # END FOR

            bxEnd.append(Bx[-1])
            byEnd.append(By[-1])
            bzEnd.append(Bz[-1])

    except IOError:
        print("Error reading neutron track file")
        return

    try:
        with open ('000000000000neutronend.out',"r") as f3:
            lines = f3.readlines()[1:]
            for num, line in enumerate(lines):
                sxStart.append(float( line.split(' ')[10]) )
                syStart.append(float( line.split(' ')[11]) )
                szStart.append(float( line.split(' ')[12]) )
                sxEnd.append(float( line.split(' ')[26]) )
                syEnd.append(float( line.split(' ')[27]) )
                szEnd.append(float( line.split(' ')[28]) )

    except IOError:
        print("Error reading neutron end file")
        return

    # Shows depolarization of neutrons
    # Find projection of S vector (unit vector) onto B vector (not unit vector)
    startProj = []      # Projection of S_start onto B_start for multiple neutrons
    endProj = []        # Projection of S_end onto B_end for multiple neutrons
    depolCount = 0
    for bxS, byS, bzS, sxS, syS, szS in zip (bxStart, byStart, bzStart, sxStart, syStart, szStart):
        if (norm(bxS, byS, bzS) == 0):
            continue
        startProj.append( (sxS*bxS + syS*byS + szS*bzS)/norm(bxS, byS, bzS) )
    for bxE, byE, bzE, sxE, syE, szE in zip (bxEnd, byEnd, bzEnd, sxEnd, syEnd, szEnd):
        if (norm(bxE, byE, bzE) == 0):
            depolCount += 1
            continue
        endProj.append( (sxE*bxE + syE*byE + szE*bzE)/norm(bxE, byE, bzE) )

    print("No. of neutrons ending in zero B field:\t", depolCount)

    fig1, ax1 = plt.subplots()
    ax1.hist(startProj, bins = nbins)
    ax1.set_xlabel('Spin projection of S on B')
    ax1.set_ylabel('Number of neutrons')
    ax1.set_title('Polarization of neutrons at start')

    fig2, ax2 = plt.subplots()
    ax2.hist(endProj, bins = nbins)
    ax2.set_xlabel('Spin projection of S on B')
    ax2.set_ylabel('Number of neutrons')
    ax2.set_title('Polarization of neutrons at end')

    plt.show()

    return

def norm (x, y, z):
    #Finds norm of cartesian vector components
    import numpy as np
    return np.sqrt(x*x + y*y + z*z)



if ( __name__ == '__main__' ):
    main()
