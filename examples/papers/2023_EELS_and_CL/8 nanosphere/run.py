import sys, os, multiprocessing, numpy as np
if sys.path[1] != os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper") : sys.path.insert(1,os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper")) #This is to import ADDA Wrapper from parent directory
import ADDAwrapper as aw

# ADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = os.path.abspath(__file__ + "/../../../../../src/seq/adda"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-4, #number of parallel processes is equal to the number of processor cores minus 1

    mp_files = [os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper/refractive_index/" + "Ag_JHW.csv")], #file with refractive index of the particle, each string contains: ev,mp_re,mp_im
    ev_range = np.arange(2,5+.01,0.05).round(2), #[eV]. Used in "spectrum_" functions. (ev_min,ev_max): range from ev_min[eV] to ev_max[eV]
    ev = 3.55, #[eV]. Used in "scan_" and "extrapolation_" functions.
    
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

import math

# ADDA command line arguments
adda_cmdlineargs = dict(
    # Particle
    shape = "sphere",
    size = 18.3, #[nm]
    grid = 128, #dipoles per x-axis size of the particle
    mhost = f"{math.sqrt(3.2)} 0", #refractive index of the host medium
    
    # Beam
    beam = "electron 120", #Energy[keV]
    beam_center = "12.4 0 0", # x[nm] y[nm] z[nm] - beam center coordinates 
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
    
    # Execute spectrum
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum")
    aw.spectrum_execute(aw_parameters,adda_cmdlineargs)
    # Collect and plot EELS probabilities
    aw.spectrum_collect("Peels",aw_parameters)
    
    from scipy.ndimage import gaussian_filter1d
    import numpy as np
    
    def fwhm2sigma(fwhm):
        return fwhm / np.sqrt(8 * np.log(2))
    
    fwhm_ev = 0.15
    fwhm = fwhm_ev*20 #20 points per eV
    sigma = fwhm2sigma(fwhm)
    
    lw = 2
    alp = 0.9
    
    fig,ax = aw.plot_create()
    
    data = np.genfromtxt("Raza.csv", delimiter=',')[1:]
    xs = data[:,0]
    ys = data[:,1]
    ys = ys/np.max(ys)
    ax.plot(xs, ys, label="Experiment", marker = ".", linestyle="none", color="black", linewidth=lw, alpha=alp)
    
    data = np.genfromtxt(aw_parameters['dirname'] + "/Peels.csv", delimiter=',')[1:]
    x = data[:,0]
    y = data[:,1]
    y = y/np.max(y)
    ax.plot(x, y, label="DDA", color=aw.color_for_plot("Peels"), linewidth=lw, alpha=alp)
    y_fwhm = gaussian_filter1d(y, sigma)
    y_fwhm = y_fwhm/np.max(y_fwhm)
    ax.plot(x, y_fwhm, label=r"DDA ($\ast$PSF)", color="purple", linewidth=lw, alpha=alp)
    #ax.fill_between(x, 0, y_fwhm, color=aw.color_for_plot("Peels"), alpha=0.2)
    
    ax.set_ylabel(r"EELS probability, arb. units")
    ax.grid(which="both", axis="x", linestyle="None")
    ax.set_xlim([min(xs),max(xs)])
    ax.legend(loc="lower right", ncol=3, bbox_to_anchor=(0., 1.02, 1., .102), borderaxespad=0.)
    aw.plot_setaspect(ax)
    fig.savefig("Peels_vs_Raza.svg", bbox_inches='tight')

