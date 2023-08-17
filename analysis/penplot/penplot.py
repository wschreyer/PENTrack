import os as os
import glob as glob
import scipy as scipy
import h5py as h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
# from pyvistaqt import BackgroundPlotter
import pyvistaqt as pvqt
import vtkplotlib as vpl
import stl as stl
from stl.mesh import Mesh
import magpylib as magpy
import magtspect as magt
from colorsys import hls_to_rgb
import numpy as np 
from datetime import datetime
# import imageio.v3 as iio
# import pygifsicle as pygifsicle
from itertools import compress

def get_distinct_colors(n):
    colors = []
    for i in np.arange(360., 0., -360. / n):
        h = i / 360.
        l = (50 + np.random.rand() * 10) / 100.
        s = (90 + np.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    return colors


plt.ion()
plt.show()

lighting = True

keydic = {
        'neutrontrack':'nt', 'neutronhit':'nh', 'neutronend':'ne', 'neutronsnapshot':'ns', 'neutronspin':'nS',
        'protontrack':'pt', 'protonhit':'ph', 'protonend':'pe', 'protonsnapshot':'ps', 'protonspin':'pS',
        'electrontrack':'et', 'electronhit':'eh', 'electronend':'ee', 'electronsnapshot':'es', 'electronspin':'eS',
}

magpy.defaults.display.style.current.arrow.width = 4

# stl_file_paths = glob.glob(workingdir + "in/CAD/*.stl") + glob.glob(workingdir + "in/CAD/store/*.stl")

            
class data:
    cbsize = 30
    
    def __init__(self, dfile="000000000002", dfilext=".h5", keys=('nt'), loadstl=True, plotter=None, ptdir="/home/sly/Work/Physics/Neutrons/tSPECT/PENtrack/pentrack/"):
        self.ptdir = ptdir
        self.keys = keys
        self.dfile = dfile
        self.dfilext = dfilext
        self.df = {}
        self.stlfiles = []
        # self.plotter.stlactors = []
        # self.plotter.logsactor = {}
        # self.pactor = None
        self.loadata()
        # self.simcount = int(str(self.df['config']['GLOBAL'][3][1]).split("'")[1].replace("\\", ""))
        self.simcount = int(str(self.df['config']["GLOBAL"][3][1]).split(" ")[1].split("\\")[0])
        # self.simtime = int(str(self.df['config']['GLOBAL'][4][1]).split("'")[1].replace("\\", ""))
        self.simtime = int(str(self.df['config']["GLOBAL"][4][1]).split(" ")[1].split("\\")[0])
        if loadstl:
            self.loadstl()
        # self.axes = None
        self.plotter = plotter
                
    def loadata(self):
        if (("h5" in self.dfilext) or ("hdf5" in self.dfilext))  :
            fh5 = h5py.File(self.ptdir + "out/" + self.dfile + ".h5",'r')
            for key in fh5.keys():
                if (key not in 'config') and (key in keydic):
                    # print("loading " + key + " in dataframe df")
                    
                    self.df[keydic[key]] = pd.DataFrame(fh5[key][:])
                else:
                    self.df["config"] = fh5["config"]

        elif (("out" in self.dfilext) or ("txt" in self.dfilext)):
            for key in keydic:
                outpath = self.ptdir + "out/" + self.dfile + key + ".out"
                if os.path.exists(outpath):
                    outpathc = outpath
                    print("loading " + key + " in dataframe df")
                    self.df[keydic[key]] = pd.read_csv(outpathc, sep=" ")

            # # Todo :  iterator for large file?
            # spinlist = []
            # data_iter = pd.read_hdf(workingdir + "out/saved/" + config_name + ".h5", key='neutronspin', columns=['x', 'y', 'z', 'Sx', 'Sy', 'Sz'], iterator=True, chunksize=100, where=None, start=50000000)
            # it = 0
            # for chunk in data_iter:
            #    if (it > 10):
            #       break
    
            #    it = it+1
            #    #train cnn on chunk here
            #    spinlist.append(chunk[['x', 'y', 'z', 'Sx', 'Sy', 'Sz']])
            #    # print(chunk[['x', 'y', 'z', 'Sx', 'Sy', 'Sz']])
            #    # print(chunk)
            #    print(spinlist)

        elif ".root" in self.dfile:
            print("not implemented")
            return

    def loadstl(self, config=None):
        if config is None:
            config = self.df["config"]
            
        geolist = config["GEOMETRY"]
        for geo in geolist:
            strgeo = str(geo[1]).split("\\t")
            if "ignored" not in strgeo[1]:
                self.stlfiles.append(self.ptdir + "in/" + strgeo[1])

    def loadallstl(self):
        '''
        Todo
        '''
        # stl_file_paths = glob.glob(self.ptdir + "in/STLs/*.stl") + glob.glob(self.ptdir + "in/STLs/store/*.stl")
        for geo in geolist:
            strgeo = str(geo[1]).split("\\t")
            if "ignored" not in strgeo[1]:
                self.allstlfiles.append(self.ptdir + "in/" + strgeo[1])

                
    class SetVisibilityCallback:
        """Helper callback to keep a reference to the actor being modified."""

        def __init__(self, actors):
            self.actors = actors

        def __call__(self, state):
            if isinstance(self.actors, (list, tuple, set, dict)):
                for actor in self.actors:
                    # if isinstance(actor, (list, tuple, set, dict)):
                    #     for actor2 in actor:
                    #         actor2.SetVisibility(state)
                    # else:
                        actor.SetVisibility(state)
            else:
                self.actors.SetVisibility(state)
                            
            
                
    class ScreenshotCallback:
        """Helper callback to keep a reference to the actor being modified."""

        def __init__(self, plotter, dfile, spath):
            self.plotter = plotter
            self.dfile = dfile
            self.spath = spath
        def __call__(self, state):
            if state:
                current_datetime = datetime.now()
                formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")
                filename = self.spath + "/screenshot-%s-%s.png" % (self.dfile, formatted_datetime)
                self.plotter.screenshot(filename)
            
    class RecordAnimCallback:
        """Helper callback to keep a reference to the actor being modified."""

        def __init__(self, plotter, dfile, spath, fps):
            self.plotter = plotter
            self.dfile = dfile
            self.spath = spath
            self.gifname = None
            self.fps = fps
        def __call__(self, state):
            self.plotter.saveanim = state
            if state:
                print("recording animation")
                current_datetime = datetime.now()
                formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")
                self.gifname = self.spath + "/animation-%s-%s" % (self.dfile.split("/")[-1], formatted_datetime)
                # self.plotter.open_gif(self.gifname+'.gif', subrectangles=True, fps=self.fps)
                self.plotter.open_movie(self.gifname+'.mp4', framerate=self.fps, quality=4)
            else:
                print("saving animation (please wait)")
                # self.plotter.open_gif(self.spath + "/dummy.gif")
                self.plotter.open_movie(self.spath + "/dummy.mp4")
                # if self.gifname is not None:
                    # pygifsicle.optimize(self.gifname) # For overwriting the original one

        
    def initplotter(self):
        plotter = pvqt.BackgroundPlotter(window_size=(1920, 1080), multi_samples=4)
        plotter.show_bounds(color="black")
        plotter.set_background('aliceblue', top='white')
        plotter.enable_depth_peeling() # Slow down animation
        return plotter

    
    def opacity_stl(self, opacity):
        for actor in self.plotter.stlactors:
            actor.GetProperty().SetOpacity(opacity)

            
    def plotstl(self, opacity=1, stlfiles=None):
        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.plotter.stlactors = []
        stlfiles = self.stlfiles if stlfiles is None else stlfiles

        pos = 2*self.cbsize
        colors = get_distinct_colors(len(stlfiles))

        for stlfile, color in zip(stlfiles, colors):
            mesh = pv.read(stlfile).scale([1000.0, 1000.0, 1000.0], inplace=False)
            # mesh = mesh.rotate_y(-90, point=axes.origin, inplace=False)
            actor = self.plotter.add_mesh(mesh, opacity=opacity, color=color, lighting=lighting)
            
            self.plotter.stlactors.append(actor)
            
            self.plotter.add_checkbox_button_widget(self.SetVisibilityCallback(actor), value=True, position=(5.0, pos), size=self.cbsize, color_on=color)
            pos = pos + self.cbsize + (self.cbsize // 10)
        
        self.plotter.add_slider_widget(self.opacity_stl, (0.0, 1.0), value=opacity, pointa=(0.04, 0.6), pointb=(0.04, 0.98),  slider_width=0.01, tube_width=0.002, color="black", style='modern')

        return self.plotter

    
    def opacity_mag(self, opacity):
        for actor in self.plotter.magactors:
            actor.GetProperty().SetOpacity(opacity)

            
    def plotmag(self, opacity=1):
        
        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.plotter.magactors = []
      
        co = magt.coil(wx=20, wr=50).source.rotate_from_angax(90, axis='y', anchor=0)
        octu = magt.octupole(nzseg=1, nring=24, stype="cyl").source.rotate_from_angax(90, axis='y', anchor=0)
        comp = magt.octupole(nzseg=1, nring=5, stype="cyl", comp=True).source.rotate_from_angax(90, axis='y', anchor=0)
        sf1 = magt.spinflipper(isf=1, kseg=20).source#.rotate_from_angax(-90, axis='y', anchor=0)
        sf2 = magt.spinflipper(isf=2, kseg=20).source#.rotate_from_angax(-90, axis='y', anchor=0)

        size = 30
        pos = 2*size
        sources = [co, octu, comp, sf1, sf2]
        prevactors = self.plotter.actors.copy()
        colors = get_distinct_colors(len(sources))
        for source, color in zip(sources, colors):
            source.set_children_styles(color=color)
            source.set_children_styles(opacity=opacity)
             
            oldvactors = self.plotter.actors.copy()
            self.plotter = source.show(backend='pyvista', canvas=self.plotter, return_fig=True, style_magnetization_show=False)
            actor = [v for k,v in self.plotter.actors.items() if k not in oldvactors]

            self.plotter.add_checkbox_button_widget(self.SetVisibilityCallback(actor), value=True, position=(40, pos), size=self.cbsize, color_on=color)
            pos = pos + size + (size // 10)

        self.plotter.magactors = [v for k,v in self.plotter.actors.items() if k not in prevactors]
        self.plotter.add_slider_widget(self.opacity_mag, (0, 1), value=opacity, pointa=(0.1, 0.6), pointb=(0.1, 0.98),  slider_width=0.01, tube_width=0.002, color="black", style='modern')

        return self.plotter


    def plotlogs(self, state="start", ptype="neutron", pselect=None, color=None):
        '''
        state = "start", "end", or "hit"
        '''
        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        
        if not hasattr(self.plotter, "logsactor"):
            self.plotter.logsactor = {}
             
        st = '' if "hit" in state else state
        ltype = 'h' if "hit" in state else 'e'

        pos = 0
        if (ptype == 'n' or ptype == "neutron"):
            key = 'n'
        elif (ptype == 'p' or ptype == "proton"):
            key = 'p'
            pos = 2*self.cbsize
        elif (ptype == 'e' or ptype == "electron"):
            key = 'e'
            pos = 3*self.cbsize
            

        if state == "start":
            color = "black" if color is None else color
        elif state == "end":
            color = "red" if color is None  else color
            pos += self.cbsize
        else:
            color = "pink" if color is None  else color
            pos += 2*self.cbsize

        df_log = da.df[key+ltype]
        
        if pselect is not None:
            df_log = df_log[df_log['particle'].isin(pselect)]
            
        positions = df_log[['x'+st, 'y'+st, 'z'+st]].values
        cloud = pv.PolyData(1000 * positions)
        
        actor = self.plotter.add_mesh(1000 * positions, color=color, render_points_as_spheres=True, point_size=6, lighting=lighting)
        self.plotter.logsactor[key+ltype+state] = actor
        
        self.plotter.add_checkbox_button_widget(self.SetVisibilityCallback(actor), value=True, position=(300 + pos, 1), size=self.cbsize, color_on=color)

        
        return self.plotter


    def speedanime(self, speed):
        self.plotter.dt = self.plotter.idt * np.sign(int(speed)) * speed**2

        
    def linlin(self, x, xmin, xmax):
        conv = (((x - xmin) * (1 - 0)) / (xmax - xmin)) + 0
        conv = conv * (1 - (conv < 0))
        conv = conv * (1 - (conv > 1)) + (conv > 1)
        return conv

    def push(self, a, n):
        a = np.roll(a, 1, axis=0)
        a[0] = n
        return a

    
    def animate(self, ti=0, tf=None, dt=0.005, fps=10, spath="/home/sly/Work/Physics/Neutrons/tSPECT/PENtrack/plots/animations/penanalyse/", pselect=None, minp=0, maxp=None, minE=None, maxE=None, minH=None, maxH=None):
        """
        Main function for animation
        """

        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.plotter.saveanim = False
        self.plotter.disable_depth_peeling()
        
        # Load and filter data
        df_nt = self.df['nt'] #.sort_values(['particle', 't'])

        maxp = self.simcount if maxp is None else maxp

        # select particle from seelct list, sort entries.
        if pselect is not None:
            df_nt = df_nt[df_nt['particle'].isin(pselect)]
        else:
            df_nt = df_nt[(df_nt['particle'] >= minp) & (df_nt['particle'] <= maxp)]

        
        df_nt = df_nt.sort_values(['particle', 't'])

        # select entries with time between bounds
        tf = tf if (tf is not None) and (tf < self.simtime) else self.simtime
        ti = ti if ti < self.simtime else self.simtime
        
        df_nt = df_nt[df_nt['t'] <= tf]
        df_nt = df_nt[df_nt['t'] >= ti]

        # plist = df_nt['particle'].unique().astype(int)
                
        df_nt_first = df_nt.groupby(['particle']).first()
        df_nt_last = df_nt.groupby(['particle']).last()

        # Get min and max kineatic (E) and total energies (H)
        minE = df_nt.min()['E'] if minE is None else minE
        maxE = df_nt.max()['E'] if maxE is None else maxE

        minH = df_nt.min()['H'] if minH is None else minH
        maxH = df_nt.max()['H'] if maxH is None else maxH

        interpolatevals = ['x', 'y', 'z', 'polarisation'] # E  H
        interpolation_functions = []
        plist = []
        for particle, data in df_nt.groupby('particle'):
            positions = data[interpolatevals].values
            times = data['t'].values
            interfunc = scipy.interpolate.interp1d(times, positions, axis=0, bounds_error=False, fill_value=(positions[0], positions[-1]))
            if interfunc is not None:
                interpolation_functions.append(interfunc)
                plist.append(particle)
        plist = np.array(plist).astype(int)
              
        # Create cloud of points
        self.plotter.pcloud = pv.PolyData(1000 * df_nt_first[['x', 'y', 'z']].values)

        # self.plotter.add_mesh(self.plotter.pcloud.delaunay_2d(), color='red', line_width=2)
        
        self.plotter.plines = []
        self.plotter.plactors = []
        for inipos in df_nt_first[['x', 'y', 'z']].values:
            pline = pv.PolyData(np.ones((10, 3)) * inipos[None])
            plactor = self.plotter.add_points(pline, color=np.random.rand(3), style='points', point_size=6)
            plactor.SetVisibility(False)
            self.plotter.plactors.append(plactor)
            self.plotter.plines.append(pline)

        self.plotter.add_checkbox_button_widget(self.SetVisibilityCallback(self.plotter.plactors), value=False, position=(100, 60), size=self.cbsize, color_on="goldenrod")

        rgba = np.ones((len(plist), 4))
        rgba[:, 0] = 120
        rgba[:, 1] = 255//2 - 100*df_nt_first['polarisation'].values
        rgba[:, 2] = 255//2 + 100*df_nt_first['polarisation'].values
        rgba[:, 3] = 1 # opacity
        self.plotter.pcloud["rgba"] = rgba
        
        
        self.plotter.pactor = self.plotter.add_points(self.plotter.pcloud, render_points_as_spheres=True, point_size=8, scalars="rgba", rgba=True, lighting=lighting)
        # self.plotter.pactor = self.plotter.add_points(cloud, point_size=8, scalars="rgba", rgba=True, lighting=lighting)
        
        self.plotter.add_checkbox_button_widget(self.SetVisibilityCallback(self.plotter.pactor), value=True, position=(100, 1), size=self.cbsize, color_on="aquamarine")

        self.plotter.pcloud["pID"] = [f"{int(i)}" for i in plist]

        self.plotter.labelpactor = self.plotter.add_point_labels(self.plotter.pcloud, "pID", font_size=14, point_size=0.1)

        self.plotter.add_checkbox_button_widget(self.SetVisibilityCallback(self.plotter.labelpactor), value=False, position=(100, 30), size=self.cbsize, color_on="burlywood")

        self.plotter.labelpactor.SetVisibility(False)        

        # # self.plotter.rgba = 'H'
        self.plotter.rgba = 'polarisation'
        # self.plotter.rgba = ''

        # self.plotter.add_checkbox_button_widget(self.RGBACallback(self.plotter), value=False, position=(100, 80), size=self.cbsize, color_on="gray", color_off='chocolate')


        self.plotter.add_checkbox_button_widget(self.RecordAnimCallback(self.plotter, self.dfile, spath, fps=fps), value=False, position=(800, 1), size=self.cbsize, color_on="maroon")
        
        self.plotter.add_checkbox_button_widget(self.ScreenshotCallback(self.plotter, self.dfile, spath), value=False, position=(800, 40), size=self.cbsize, color_on="gray", color_off='lime')

        self.plotter.add_checkbox_button_widget(self.ScreenshotCallback(self.plotter, self.dfile, spath), value=False, position=(800, 40), size=self.cbsize, color_on="gray", color_off='lime')

        
        self.plotter.time_stamp = ti
        self.plotter.idt = dt
        self.plotter.titletxt = self.plotter.add_text(self.dfile, position='upper_right', font_size=18, color="black")
        self.plotter.timetxt = self.plotter.add_text("t=%g [s]"%self.plotter.time_stamp, position='lower_right', font_size=18, color="black")
        
        def update_animation():
            if ((self.plotter.dt < 0) and (self.plotter.time_stamp > ti)) or ((self.plotter.dt > 0) and (self.plotter.time_stamp < tf)):
                self.plotter.time_stamp += self.plotter.dt
                self.plotter.timetxt.SetText(1, "t={:8.3f} [s]".format(self.plotter.time_stamp))

                p2up = np.array(self.plotter.time_stamp < df_nt_last['t']).astype(bool)                    

                ####### based on preselection. Beta
                
                # interpolated_vals = []
                # # for interpolation_func, pline in zip(interpolation_functions, self.plotter.plines):
                # for interpolation_func in compress(interpolation_functions, p2up):
                #     valinter = interpolation_func(self.plotter.time_stamp)                #     interpolated_vals.append(valinter)
                #     # pline.points = self.push(pline.points, 1000 * valinter[:3])
           
                # interpolated_vals = np.array(interpolated_vals)
                # newpoints = 1000 * interpolated_vals[:, 0:3]
                # self.plotter.pcloud.points[p2up] = newpoints

                # if self.plotter.plactors[0].GetVisibility():
                #     for particle, pline in enumerate(list(self.plotter.plines[i-1] for i in plist[p2up])):
                #         self.push(pline.points, newpoints[particle-1])

                ###############
                
                interpolated_vals = []
                for interpolation_func in interpolation_functions:
                    valinter = interpolation_func(self.plotter.time_stamp)
                    interpolated_vals.append(valinter)
                    # pline.points = self.push(pline.points, 1000 * valinter[:3])
           
                interpolated_vals = np.array(interpolated_vals)
                newpoints = 1000 * interpolated_vals[:, 0:3]
                self.plotter.pcloud.points = newpoints

                if self.plotter.plactors[0].GetVisibility():
                    for particle, pline in enumerate(self.plotter.plines):
                        pline.points = self.push(pline.points, newpoints[particle])


                ###########
                
                rgba = self.plotter.pcloud["rgba"]
                if self.plotter.rgba == 'polarisation':
                    rgba[:,1] = 255//2 - 100*interpolated_vals[:, 3]
                    rgba[:,2] = 255//2 + 100*interpolated_vals[:, 3]
                    
                elif self.plotter.rgba == 'E':
                    rgba[:,1:3] = 255 * self.linlin(interpolated_vals[:, 4], minE, maxE)[:, None]
                    
                elif self.plotter.rgba == 'H':
                    rgba[:,0] = 255 * self.linlin(interpolated_vals[:, 5], minH, maxH)
                    rgba[:,1] = 255 - 255 * self.linlin(interpolated_vals[:, 5], minH, maxH)

                rgba[:, 3] = p2up
                self.plotter.pcloud["rgba"] = rgba

                # self.plotter.render()
                self.plotter.update()
                    
                if self.plotter.saveanim:
                    self.plotter.write_frame()
                    

        self.plotter.add_checkbox_button_widget(self.speedanime, value=True, position=(900, 1), size=self.cbsize, color_on="aquamarine")

        self.plotter.add_slider_widget(self.speedanime, (-5, 5), value=1, pointa=(0.5, 0.03), pointb=(0.8, 0.03), slider_width=0.02, tube_width=0.005, color="aquamarine", style='modern', title_opacity=0)

        self.plotter.add_callback(update_animation, int(1000/fps), None)
        
        return self.plotter



# name of datafile
# dfile = "000000000105"
# dfile = "000000000112"
# dfile = "000000000017"
dfile = "000000000034"
# dfile = "/saved/sf-conductor62"

# dfile = "000000000201"
# instantiate data object
da = data(dfile)

# plots stl surface and magnetic source
pl = da.plotmag(opacity=0.01)
pl = da.plotstl(opacity=0.01)

df_ne = da.df['ne']

# pselect = df_ne[df_ne['xend'] > df_ne['xstart'] + 0.5]['particle']
pselect = None
# pselect = np.arange(1, 20)
     
# plots neutrons start, end, and hits point.
pl = da.plotlogs(ptype="n", state="start", pselect=pselect, color="lightgreen")
pl = da.plotlogs(ptype="n", state="end", pselect=pselect, color="deeppink")
pl = da.plotlogs(ptype="n", state="hit", pselect=pselect, color="deepskyblue")

# # play animation
# pl = da.animate(ti=0, tf=None, dt=0.001, fps=10, pselect=pselect, minp=0, maxp=None, minE=0, maxE=5e-8, minH=0, maxH=3e-7)

pl = da.animate(ti=0, tf=10, dt=0.002, fps=20, pselect=pselect, minp=0, maxp=None, minE=0, maxE=5e-8, minH=0, maxH=3e-7)



# plt.figure()
# [plt.plot(da.df["nt"][da.df["nt"]["particle"] == p]['t'], da.df["nt"][da.df["nt"]["particle"] == p]['H']) for p in da.df['ne']["particle"]]




# [plt.plot(da.df['nt'][da.df['nt']['particle'] == p]['x'], da.df['nt'][da.df['nt']['particle'] == p]['polarisation']) for p in np.arange(10)]





# dfs = da.df['nt'][da.df['nt']['particle'] == 31]
# plt.figure()
# plt.subplot(2, 2, 1)
# plt.plot(dfs['t'], dfs['H'])
# plt.xlabel('t')
# plt.ylabel('H')


# plt.subplot(2, 2, 2)
# plt.plot(dfs['t'], dfs['E'])
# plt.xlabel('t')
# plt.ylabel('E')


# plt.subplot(2, 2, 3)
# plt.plot(dfs['t'], dfs['dBxdx'] + dfs['dBydy'] + dfs['dBzdz'])
# plt.xlabel('t')
# plt.ylabel(r'$\sum_i dB_i/dx_i$')


# plt.subplot(2, 2, 4)
# plt.plot(dfs['t'], dfs['dBxdy'] - dfs['dBydx'])
# plt.plot(dfs['t'], dfs['dBxdz'] - dfs['dBzdx'])
# plt.plot(dfs['t'], dfs['dBydz'] - dfs['dBzdy'])
# plt.ylabel(r'$dB_i/dx_j - dB_j/dx_i$')






#####################


# plt.figure()
# plt.hist(da.df['nt']['polarisation'])


