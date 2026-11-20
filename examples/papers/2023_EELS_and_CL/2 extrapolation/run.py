import sys, os, multiprocessing
if sys.path[1] != os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper") : sys.path.insert(1,os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper")) #This is to import ADDA Wrapper from parent directory
import ADDAwrapper as aw

# ADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = os.path.abspath(__file__ + "/../../../../../src/seq/adda"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-1, #number of parallel processes is equal to the number of processor cores minus 1

    mp_files = [os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper/refractive_index/" + "Ag_Palik_Garcia.csv")], #file with refractive index of the particle, each string contains: ev,mp_re,mp_im
    ev_range = (1,5), #[eV]. Used in "spectrum_" functions. (ev_min,ev_max): range from ev_min[eV] to ev_max[eV]
    ev = 3.6, #[eV]. Used in "scan_" and "extrapolation_" functions.
    
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
    mhost = "1 0", #refractive index of the host medium
    
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
    
    # Execute extrapolation for single energy ev
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "extrapolation_3.6")
    aw.extrapolation_execute(aw_parameters,adda_cmdlineargs)
    # Collect and plot EELS probabilities for different y~1/grid with extrapolated value in y=0 (with errorbar)
    aw.extrapolation_collect("Peels",aw_parameters)
    #aw.extrapolation_plot("Peels",aw_parameters)
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv",delimiter=',')[1:]
    fig,ax = aw.plot_create()
    ax.plot((-1,1), (0.003454231,0.003454231), label="EELS Lorenz-Mie", color="black", linestyle="--", linewidth=3, zorder=3)
    ax.plot(data[:,1], data[:,2], label=r"P$_{\rmEELS}$ (simulated)", color=aw.color_for_plot("Peels"), marker="o", markersize=12, linestyle="none", zorder=2)
    ys_fitted = np.linspace(data[:,1][0],0,100)
    results_fit = np.genfromtxt(f"{aw_parameters['dirname']}/Peels_fit.csv",delimiter=',')[1:]
    a = results_fit[:,0]
    error = results_fit[:,1]
    points_fitted = a[0] + a[1]*ys_fitted + a[2]*ys_fitted**2
    ax.plot(ys_fitted, points_fitted, label=r"P$_{\rmEELS}$ (fit)", color="black", linewidth=3, zorder=3)
    ax.errorbar(0, a[0], yerr=error[0], color="black", linestyle="", marker="s", capsize=3, barsabove=True, label = "Confidence interval")
    ax.set_xlabel("y = kd|m|")
    ax.set_xlim([-0.005,max(data[:,1])*1.05])
    #ax.grid(False)
    ax.legend()
    ax.set_ylabel("Probability, eV$^{-1}$")
    aw.plot_setaspect(ax)
    plt.savefig("Figure.svg", bbox_inches='tight')
    


