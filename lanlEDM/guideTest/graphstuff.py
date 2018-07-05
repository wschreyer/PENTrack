#!/usr/bin/env python
# Makes some quick and dirty plots
# Accepts flag arguments so that user isn't spammed with plots
# Use [--help] or [-h] to see options

def main():
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import argparse
    import sys

    nbins = 200        # number of histogram bins for polarization plots
    sx = []         # Spin vectors
    sy = []
    sz = []
    tSpin = []      # time interval (corresponding to spin vectors)
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

    # Read in command line flags
    parser = argparse.ArgumentParser(description='Plot some graphs from .out files')
    parser.add_argument("-s", "--spin", help = "Prints graphs related to neutron spin", action="store_true")
    parser.add_argument("-f", "--fourier", help = "Prints fourier transform of spin components", action="store_true")
    parser.add_argument("-b", "--bfield", help = "Prints magnetic field seen by neutron", action="store_true")
    parser.add_argument("-t", "--traj", help = "Prints trajectory of neutron", action="store_true")
    parser.add_argument("-p", "--pol", help = "Prints initial and final polarization of the neutron(s)", action="store_true")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        print("No flags given. Use [-h] or [--help]")
        sys.exit(1)

    # Read in files
    try:
        with open ('000000000000neutronspin.out',"r") as f1:
            lines = f1.readlines()[1:]
            for num, line in enumerate(lines):
                tSpin.append(float( line.split(' ')[2]) )
                sx.append(float( line.split(' ')[3]) )
                sy.append(float( line.split(' ')[4]) )
                sz.append(float( line.split(' ')[5]) )
    except IOError:
        print("Error reading neutron spin file")
        return

    try:
        with open ('000000000000neutrontrack.out',"r") as f2:
            lines = f2.readlines()[1:]
            for num, line in enumerate(lines):
                Bx.append(float( line.split(' ')[12]) )
                By.append(float( line.split(' ')[16]) )
                Bz.append(float( line.split(' ')[20]) )
                tBField.append(float( line.split(' ')[3]) )
                x.append(float( line.split(' ')[4]) )
                y.append(float( line.split(' ')[5]) )
                z.append(float( line.split(' ')[6]) )
                Bnorm.append( norm(Bx[num], By[num], Bz[num]))

                 # For each neutron, append B field start/end values
                 # Kinda ugly method of doing this but whatever
                if neutronCounter < int( line.split(' ')[1]):
                    bxStart.append(Bx[num])
                    byStart.append(By[num])
                    bzStart.append(Bz[num])
                    if neutronCounter != 0:
                        bxEnd.append(Bx[num-1])
                        byEnd.append(By[num-1])
                        bzEnd.append(Bz[num-1])
                    neutronCounter += 1
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

    # Plot stuff
    if args.spin:       # Spin related graphs
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
        plt.plot(tSpin, sz)
        plt.grid(True)
        plt.title('t vs Sz')
        plt.xlabel('t [sec]')
        plt.ylabel('P(z)')

        fig3 = plt.figure(3)
        plt.plot(tSpin, sy)
        plt.grid(True)
        plt.title('t vs Sy')
        plt.xlabel('t [sec]')
        plt.ylabel('P(y)')

        fig4 = plt.figure(4)
        plt.plot(tSpin, sx)
        plt.grid(True)
        plt.title('t vs Sx')
        plt.xlabel('t [sec]')
        plt.ylabel('P(x)')

    if args.fourier:        # Fourier transform graphs
        # Calculate fourier transform of oscillation in Sy and Sz
        timestep = tSpin[1]
        fourSy = np.fft.fft(sy)
        fourSz = np.fft.fft(sz)
        freqs = np.fft.fftfreq( len(sy), d=timestep )

        fig5 = plt.figure(5)
        plt.plot(freqs, np.abs(fourSy) )
        plt.grid(True)
        plt.xlabel('[Hz]')
        plt.title('Fourier transform of Sy')

        fig6 = plt.figure(6)
        plt.plot(freqs, np.abs(fourSz) )
        plt.grid(True)
        plt.xlabel('[Hz]')
        plt.title('Fourier transform of Sz')

    if args.traj:            # prints trajectory of neutron
        fig7, ax7 = plt.subplots(3, sharex=True)
        ax7[0].plot(tBField,x)
        ax7[0].set_title('Neutron Position')
        ax7[0].set_ylabel('x [m]')
        ax7[1].plot(tBField,y)
        ax7[1].set_ylabel('y [m]')
        ax7[2].plot(tBField,z)
        ax7[2].set_ylabel('z [m]')
        ax7[2].set_xlabel('t [sec]')

        fig8 = plt.figure(8)
        ax8 = fig8.add_subplot(111, projection='3d')
        plt.title('Neutron Position')
        ax8.scatter(x, y, z)
        ax8.set_xlim3d(-2.5,2.5)
        ax8.set_ylim3d(-0.2,0.2)
        ax8.set_zlim3d(-0.2,0.2)
        ax8.set_xlabel('x [m]')
        ax8.set_ylabel('y [m]')
        ax8.set_zlabel('z [m]')
        ax8.view_init(30,220)

    if args.bfield:         # prints Bfield seen by the neutron(s)
        fig9 = plt.figure(9)
        plt.plot(tBField, Bnorm)
        plt.grid(True)
        plt.title('t vs |B|')
        plt.xlabel('t')
        plt.ylabel('|B| [Tesla]')

    if args.pol:            # Shows depolarization of neutrons
        # Find projection of S vector (unit vector) onto B vector (not unit vector)
        startProj = []      # Projection of S_start onto B_start for multiple neutrons
        endProj = []        # Projection of S_end onto B_end for multiple neutrons
        for bxS, byS, bzS, sxS, syS, szS in zip (bxStart, byStart, bzStart, sxStart, syStart, szStart):
            startProj.append( (sxS*bxS + syS*byS + szS*bzS)/norm(bxS, byS, bzS) )
        for bxE, byE, bzE, sxE, syE, szE in zip (bxEnd, byEnd, bzEnd, sxEnd, syEnd, szEnd):
            endProj.append( (sxE*bxE + syE*byE + szE*bzE)/norm(bxE, byE, bzE) )

        fig10, ax10 = plt.subplots()
        ax10.hist(startProj, bins = nbins)
        ax10.set_xlabel('Spin projection of S on B')
        ax10.set_ylabel('Number of neutrons')
        ax10.set_title('Polarization of neutrons at start')

        fig11, ax11 = plt.subplots()
        ax11.hist(endProj, bins = nbins)
        ax11.set_xlabel('Spin projection of S on B')
        ax11.set_ylabel('Number of neutrons')
        ax11.set_title('Polarization of neutrons at end')

    plt.show()

    return

def norm (x, y, z):
    #Finds norm of cartesian vector components
    import numpy as np
    return np.sqrt(x*x + y*y + z*z)



if ( __name__ == '__main__' ):
    main()
