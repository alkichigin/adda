import sys, os, multiprocessing, numpy as np
if sys.path[1] != os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper") : sys.path.insert(1,os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper")) #This is to import ADDA Wrapper from parent directory
import ADDAwrapper as aw

# ADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = aw.addaexec_find(mode="seq"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-1, #number of parallel processes is equal to the number of processor cores minus 1

    mp_files = [os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Au_JHW.csv")], #file with refractive index of the particle, each string contains: ev,mp_re,mp_im
    ev_range = np.arange(0.5,2.5+.01,0.05), #[eV]. Used in "spectrum_" functions. (ev_min,ev_max): range from ev_min[eV] to ev_max[eV]
    ev = 3.55, #[eV]. Used in "scan_" and "extrapolation_" functions.
    
    #Used in "scan_" functions. Beam propagation must be orthogonal to the grid.
    #So "prop" must be "0 0 whatever" and rotations with "orient" must be made by 90 degrees.
    scan_x_range = (-110,120), #[nm], (x_left, x_right)
    scan_y_range = (-130,130), #[nm], (y_bottom, y_top)
    scan_step = 1 #dipoles per each step, must be an integer >= 1
    #The beam must always blast exactly in the middle between the dipoles,
    #so start and stop coordinates will be adjusted, covering more area than you entered. Obligatory to use "no_vol_cor" with scan.
)

import math
Dx = (209*math.sqrt(3)/2)
#print(Dx)

# ADDA command line arguments
adda_cmdlineargs = dict(
    # Particle
    shape = f"prism 3 {10/Dx}",
    size = Dx, #[nm]
    grid = 128, #dipoles per x-axis size of the particle
    surf = "5 2.15 0",
    
    # Beam
    beam = "electron 200", #Energy[keV]
    beam_center = f"{Dx/2} 0 0", # x[nm] y[nm] z[nm] - beam center coordinates 
    prop = "0 0 -1", #beam propagation direction vector
    
    # Precision and performance
    eps = 2, #Residual norm
    
    # Additional options
    sym = "enf", #Do not simulate second polarization
    scat_matr = "none", #Do not calculate the Mueller matrix
    no_vol_cor = "", #Disable volume correction
    iter = "qmr2", #Iterative solver
    pol = "igt_so", #Polarizability prescription
    int = "igt_so", #Interaction term
    Csca = "", #Calculate Csca with the Romberg integral. Needed to properly calculate Cathodoluminesce
    alldir_inp = os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/Csca_upper_integration.txt")
)

### Executing commands
if __name__ == '__main__': 
    
    import csv
    import numpy as np
    import matplotlib.pyplot as plt
    
    # aw_parameters["dirname"] = os.path.abspath(__file__ + "/../")
    # grid = adda_cmdlineargs["grid"]
    # aw.geom_gen(aw_parameters,adda_cmdlineargs)
    # shape1 = adda_cmdlineargs["shape"].split()[0]
    # data = np.genfromtxt(f"geom/{shape1}.geom", delimiter=' ', dtype="int")
    # xs1, ys1, zs1 = data[:,0], data[:,1], data[:,2]
    # ms1 = np.full(shape=len(xs1),fill_value=2,dtype="int")
    # # xs1, ys1, zs1 = xs1, zs1, -ys1 #Mx(90)
    # # xs1, ys1, zs1 = ys1, -xs1, zs1 #Mz(90)
    # xs1 = xs1 - min(xs1)
    # ys1 = ys1 - min(ys1)
    # zs1 = zs1 - min(zs1)
    # with open("geom/shape1.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerows(zip(xs1,ys1,zs1))
    
    # times = 2
    # laynum = int(30*adda_cmdlineargs["grid"]/adda_cmdlineargs["size"])
    # zs1 = zs1 + laynum
    # xs0, ys0, zs0 = [], [], []
    # for i in range(laynum):
    #     for j in range(0,times*(max(ys1)+1)):
    #         for k in range(0,times*(max(xs1)+1)):
    #             xs0.append(k),ys0.append(j),zs0.append(i)
    # xs0 = np.asarray(xs0) - int((times-1)*max(xs1+1)/2)
    # ys0 = np.asarray(ys0) - int((times-1)*max(ys1+1)/2)
    # ms0 = np.full(shape=len(xs0),fill_value=1,dtype="int")
    
    # xs, ys, zs = np.concatenate((xs0, xs1), axis=None), np.concatenate((ys0, ys1), axis=None), np.concatenate((zs0, zs1), axis=None)
    # zs = zs - min(zs)
    # ms = np.concatenate((ms0, ms1), axis=None)
    # with open("geom/sum.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerow(["Nmat=2"])
    #     writer.writerows(zip(xs,ys,zs,ms))


    # fig1,ax1 = aw.geometry("geom/sum.geom", ["deepskyblue","gold"])
    # ax1.view_init(30,120)
    # plt.axis('off')
    # fig1.savefig("sum.svg", bbox_inches='tight')
    # fig1.savefig("sum.png", bbox_inches='tight', dpi=1200)
    
    # fig1,ax1 = aw.geometry("geom/prism.geom", ["gold"])
    # ax1.view_init(30,120)
    # plt.axis('off')
    # fig1.savefig("prism.svg", bbox_inches='tight')
    # fig1.savefig("prism.png", bbox_inches='tight',dpi=1200)

    
    
    # lw = 2
    # alp = 0.8
    # fig,ax = aw.plot_create()
    
    # adda_cmdlineargs["beam_center"] = f"{Dx/2} {209/2} 0"
    # aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_side")
    # aw.spectrum_execute(aw_parameters,adda_cmdlineargs)
    # aw.spectrum_collect("Peels",aw_parameters)
    # data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv", delimiter=',')[1:]
    # ax.plot(data[:,0], data[:,1], label="Новая теория", color="royalblue", linewidth=lw, alpha=alp)
    

    # ax.grid(which="both", axis="x", linestyle="None")
    # ax.set_xlim([min(data[:,0]),max(data[:,0])])
    # ax.set_ylabel("Probability, eV$^{-1}$")
    # #ax.grid(which="both", axis="x", linestyle="dotted", zorder=0)
    # #ax.grid(False)
    # ax.legend()
    # aw.plot_setaspect(ax)
    # fig.savefig("result.svg", bbox_inches='tight')
    
    aw_parameters["dirname"] = os.path.abspath(__file__ + "/../" + "spectrumline")
    aw_parameters["spectrumline_startpoint"] = (Dx/2,0) # (x,y) [nm]
    aw_parameters["spectrumline_endpoint"] = (Dx/2,209/2) # (x,y) [nm]
    aw_parameters["spectrumline_points"] = 10 # how many points, including startpoint and endpoint
    aw.spectrumline_execute(aw_parameters,adda_cmdlineargs) # Execute simulation
    aw.spectrumline_collect("Peels",aw_parameters) # Collect EELS spectrum
    aw.spectrumline_plot("Peels",aw_parameters) # Plot EELS spectrum
    aw.spectrumline_collect("Pcl",aw_parameters) # Collect CL spectrum
    aw.spectrumline_plot("Pcl",aw_parameters) # Plot CL spectrum
    
    # aw_parameters["parallel_procs"] = 20
    
    # aw_parameters["dirname"] = "/mnt/scratch/ws/aakichigin/202207291312ws1/scan1.75"
    # aw_parameters["ev"] = 1.75
    # aw.scan_execute(aw_parameters,adda_cmdlineargs)
    # aw.scan_collect("Peels",aw_parameters)
    # aw.scan_plot("Peels",aw_parameters, details=False)
    
    # aw_parameters["dirname"] = "/mnt/scratch/ws/aakichigin/202207291312ws1/scan1.4"
    # aw_parameters["ev"] = 1.4
    # aw.scan_execute(aw_parameters,adda_cmdlineargs)
    # aw.scan_collect("Peels",aw_parameters)
    # aw.scan_plot("Peels",aw_parameters, details=False)
    
    # aw_parameters["dirname"] = "/mnt/scratch/ws/aakichigin/202207291312ws1/scan1.25"
    # aw_parameters["ev"] = 1.25
    # aw.scan_execute(aw_parameters,adda_cmdlineargs)
    # aw.scan_collect("Peels",aw_parameters)
    # aw.scan_plot("Peels",aw_parameters, details=False)
    
    # aw_parameters["dirname"] = "/mnt/scratch/ws/aakichigin/202207291312ws1/scan0.85"
    # aw_parameters["ev"] = 0.85
    # aw.scan_execute(aw_parameters,adda_cmdlineargs)
    # aw.scan_collect("Peels",aw_parameters)
    # aw.scan_plot("Peels",aw_parameters, details=False)
    
    


