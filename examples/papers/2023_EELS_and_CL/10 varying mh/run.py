import sys, os, multiprocessing
if sys.path[1] != os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper") : sys.path.insert(1,os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper")) #This is to import ADDA Wrapper from parent directory
import ADDAwrapper as aw

# ADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = os.path.abspath(__file__ + "/../../../../../src/seq/adda"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-1, #number of parallel processes is equal to the number of processor cores minus 1

    mp_file = os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper/refractive_index/" + "Ag_Palik_Garcia.csv"), #file with refractive index of the particle, each string contains: ev,mp_re,mp_im
    ev_range = (1,5), #[eV]. Used in "spectrum_" functions. (ev_min,ev_max): range from ev_min[eV] to ev_max[eV]
    ev = 2., #[eV]. Used in "scan_" and "extrapolation_" functions.
    
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
    grid = 32, #dipoles per x-axis size of the particle
    
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
    #int = "igt_so", #Interaction term
    Csca = "", #Calculate Csca with the Romberg integral. Needed to properly calculate Cathodoluminesce
    alldir_inp = os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper/Csca_default.txt")
)

### Executing commands
if __name__ == '__main__': 

    import numpy as np
    import math
    lw = 2
    alp = 0.8
    r = 510.99895
    vc = math.sqrt(1-math.pow(r/(100+r),2))
    fig,ax = aw.plot_create()
    
    ax.axvline(x=1/vc, color="black", linestyle="dotted")

    a = np.linspace(1.5, 2.5, 201)
    b = [(str(round(i,3)) + " 0") for i in a]
    # b = ["1 " + (str(round(i,2))) for i in a]
    # #print(b)
    
    adda_cmdlineargs["lambda"] = aw.ev_to_nm(2)
    adda_cmdlineargs["m"] = "0.130999989 3.880000094"
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "varyany_mh_ev" + str(aw_parameters["ev"]))
    aw.varyany_execute(aw_parameters,adda_cmdlineargs,"mhost",b)
    aw.varyany_collect("Peels",aw_parameters)
    aw.varyany_collect("Pext",aw_parameters)
    
    data1 = np.genfromtxt(f"{aw_parameters['dirname']}/Pext.csv", delimiter=',', dtype="S")[1:]
    data1 = np.char.decode(data1)
    xs = [float(i.split(" ")[0]) for i in data1[:,0]]
    data2 = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv", delimiter=',')[1:]
    ax.plot(xs, data2[:,1].astype("float"), label="P" + r"$_{\rm EELS}$", color="royalblue", linewidth=2*lw, alpha=alp)
    ax.plot(xs, data1[:,1].astype("float"), label="P" + r"$_{\rm ext}$", color="black", linewidth=lw, alpha=alp)
    
    xs = np.asarray(xs)
    
    x1 = 1.82
    y1 = 2.5989760563E-03
    x2 = 1.83
    y2 = 2.7554533887E-03
    k = (y2-y1)/(x2-x1)
    b = y1 - k*x1
    y0 = k*(1/vc) + b
    
    xs = xs[xs > 1.82]
    #print(xs)
    ys = y0/(vc*xs)**2
    ax.plot(xs, ys, label="P" + r"$_{\rm EELS}$ scaled", color="red", linewidth=lw, alpha=alp)
    
    
    ax.set_xlabel(r"Refractive index of the host medium $\it{m}_{\rmh}$")
    ax.set_ylabel("Probability, eV$^{-1}$") 
    ax.set_xlim([1.5,2.5])
    ax.legend()
    aw.plot_setaspect(ax)
    ax.grid(which="both", axis="x", linestyle="None")
    fig.savefig("Figure.svg", bbox_inches='tight')
    fig.savefig(f"{str(aw_parameters['ev'])}.svg", bbox_inches='tight')



