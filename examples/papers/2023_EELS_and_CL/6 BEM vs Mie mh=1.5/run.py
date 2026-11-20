import sys, os, multiprocessing, numpy as np
if sys.path[1] != os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper") : sys.path.insert(1,os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper")) #This is to import ADDA Wrapper from parent directory
import ADDAwrapper as aw

# ADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = os.path.abspath(__file__ + "/../../../../../src/seq/adda"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-4, #number of parallel processes is equal to the number of processor cores minus 1

    mp_files = [os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper/refractive_index/" + "Ag_Palik_Garcia.csv")], #file with refractive index of the particle, each string contains: ev,mp_re,mp_im
    ev_range = np.arange(1,5+.01,0.05).round(2), #[eV]. Used in "spectrum_" functions. (ev_min,ev_max): range from ev_min[eV] to ev_max[eV]
    ev = 3.0, #[eV]. Used in "scan_" and "extrapolation_" functions.
    
    spectrumline_startpoint = (0,50), # (x,y) [nm]
    spectrumline_endpoint = (50,50), # (x,y) [nm]
    spectrumline_points = 10, #how many points, including startpoint and endpoint
    
    #Used in "scan_" functions. Beam propagation must be orthogonal to the grid.
    #So "prop" must be "0 0 whatever" and rotations with "orient" must be made by 90 degrees.
    scan_x_range = (0,30), #[nm], (x_left, x_right)
    scan_y_range = (0,30), #[nm], (y_bottom, y_top)
    scan_step = 1 #dipoles per each step, must be an integer >= 1
    #The beam must always blast exactly in the middle between the dipoles,
    #so start and stop coordinates will be adjusted, covering more area than you entered. Obligatory to use "no_vol_cor" with scan.
)

# ADDA command line arguments
adda_cmdlineargs = dict(
    # Particle
    shape = "sphere",
    size = 150, #[nm]
    grid = 128, #dipoles per x-axis size of the particle
    mhost = "1.5 0", #refractive index of the host medium
    
    # Beam
    beam = "electron 100", #Energy[keV]
    beam_center = "100 0 0", # x[nm] y[nm] z[nm] - beam center coordinates 
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
    alldir_inp = os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper/Csca_default.txt")
)

### Executing commands
if __name__ == '__main__': 
    
    # Execute spectrum with extrapolation at each energy ev
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_with_extrapolation_mh1.5")
    aw.spectrum_with_extrapolation_execute(aw_parameters,adda_cmdlineargs)
    # Collect and plot fit EELS probabilities spectrum (with errorbar)
    aw.spectrum_with_extrapolation_collect("Peels",aw_parameters)
    aw.spectrum_with_extrapolation_collect("Pcl",aw_parameters)
    
    # adda_cmdlineargs["grid"] = 192
    # aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "extrapolation_3.0")
    # aw.extrapolation_execute(aw_parameters,adda_cmdlineargs)
    # aw.extrapolation_collect("Peels",aw_parameters)
    # aw.extrapolation_plot("Peels",aw_parameters)

    import numpy as np

    lw = 2
    alp = 0.8
    
    fig,ax = aw.plot_create()
    
    data = np.genfromtxt("EELS_mie_mh1.5.csv", delimiter=',')[1:]
    ax.plot(data[:,0], data[:,1], label="Lorenz-Mie\n(scaled)", color="black", linewidth=lw, linestyle="--", alpha=alp)
    
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels_fit.csv", delimiter=',')[1:]
    ax.plot(data[:,0], data[:,1], label="DDA", color=aw.color_for_plot("Peels"), linewidth=lw, alpha=alp)
    
    data = np.genfromtxt("EELS_bem_mh1.5.csv", delimiter=',')[1:]
    ax.plot(data[:,0], data[:,1], label="BEM", color="green", linestyle="none", marker="8", fillstyle="none", linewidth=lw, alpha=alp)
    
    ax.grid(which="both", axis="x", linestyle="None")
    ax.set_xlim([min(data[:,0]),max(data[:,0])])
    ax.set_ylabel("EELS probability, eV$^{-1}$")
    #ax.grid(False)
    if min(data[:,1]) > 0:
        ax.set_ylim(bottom = 0)
    ax.legend(loc="upper left")
    aw.plot_setaspect(ax)
    fig.savefig("EELS.svg", bbox_inches='tight')
    
    
    
    fig,ax = aw.plot_create()
    
    data = np.genfromtxt("CL_mie_mh1.5.csv", delimiter=',')[1:]
    ax.plot(data[:,0], data[:,1], label="Lorenz-Mie (scaled)", color="black", linewidth=lw, linestyle="--", alpha=alp)
    
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Pcl_fit.csv", delimiter=',')[1:]
    ax.plot(data[:,0], data[:,1], label="DDA", color=aw.color_for_plot("Pcl"), linewidth=lw, alpha=alp)
    
    data = np.genfromtxt("CL_bem_mh1.5.csv", delimiter=',')[1:]
    ax.plot(data[:,0], data[:,1], label="BEM", color="green", linestyle="none", marker="8", fillstyle="none", linewidth=lw, alpha=alp)
    ax.plot(data[:,0], data[:,1]*1.5, label="Corrected BEM", color="green", linestyle="none", marker="x", fillstyle="none", linewidth=lw, alpha=alp)
    
    ax.grid(which="both", axis="x", linestyle="None")
    ax.set_xlim([min(data[:,0]),max(data[:,0])])
    ax.set_ylabel("CL probability, eV$^{-1}$")
    #ax.grid(False)
    if min(data[:,1]) > 0:
        ax.set_ylim(bottom = 0)
    ax.legend(loc="lower right")
    aw.plot_setaspect(ax)
    fig.savefig("CL.svg", bbox_inches='tight')

