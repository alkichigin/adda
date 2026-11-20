import sys, os, multiprocessing, numpy as np
if sys.path[1] != os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper") : sys.path.insert(1,os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper")) #This is to import ADDA Wrapper from parent directory
import ADDAwrapper as aw

# ADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = os.path.abspath(__file__ + "/../../../../../src/seq/adda"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-4, #number of parallel processes is equal to the number of processor cores minus 1

    mp_files = [os.path.abspath(__file__ + "/../../../../../misc/ADDAwrapper/refractive_index/" + "Ag_Palik_Garcia.csv")], #file with refractive index of the particle, each string contains: ev,mp_re,mp_im
    ev_range = np.arange(1,5+.01,0.05), #[eV]. Used in "spectrum_" functions. (ev_min,ev_max): range from ev_min[eV] to ev_max[eV]
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
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh1")
    aw.spectrum_execute(aw_parameters,adda_cmdlineargs)
    aw.spectrum_collect("Peels",aw_parameters)
    aw.spectrum_collect("Pcl",aw_parameters)
    
    adda_cmdlineargs["mhost"] = "1.5 0"
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh1.5")
    aw.spectrum_execute(aw_parameters,adda_cmdlineargs)
    aw.spectrum_collect("Peels",aw_parameters)
    aw.spectrum_collect("Pcl",aw_parameters)
    
    adda_cmdlineargs["mhost"] = "2 0"
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh2")
    aw.spectrum_execute(aw_parameters,adda_cmdlineargs)
    aw.spectrum_collect("Peels",aw_parameters)
    aw.spectrum_collect("Pcl",aw_parameters)
    aw.spectrum_collect("Pcltot",aw_parameters)
    
    # adda_cmdlineargs["mhost"] = "3 0"
    # aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh3")
    # #aw.spectrum_execute(aw_parameters,adda_cmdlineargs)
    # aw.spectrum_collect("Peels",aw_parameters)
    # aw.spectrum_collect("Pcl",aw_parameters)
    # aw.spectrum_collect("Pcltot",aw_parameters)
    

    import numpy as np

    lw = 2
    alp = 1
    
    fig1,ax1 = aw.plot_create()
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh1")
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv", delimiter=',')[1:]
    ax1.plot(data[:,0], data[:,1], label=""+r"${\it m}$"+r"$_{\rm h}$=1", color=aw.color_for_plot("Peels"), linewidth=lw, alpha=alp)
    
    ax1.plot(np.NaN, np.NaN, color='none', label=' ')
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh1.5")
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv", delimiter=',')[1:]
    ax1.plot(data[:,0], data[:,1], label=""+r"${\it m}$"+r"$_{\rm h}$=1.5", color="purple", linestyle="-", linewidth=lw, alpha=alp)
    
    
    
    
    
    
    ##### SCALING
    
    import math
    r = 510.99895
    K = float(adda_cmdlineargs["beam"].split()[1])
    
    aw_parameters["mp_files"] = [os.path.abspath(__file__ + "/../" + "Ag_Palik_Garcia_div1.5.csv")]
    mh = 1.5
    v_c = math.sqrt(1-(r/(K+r))**2)
    KK = r/math.sqrt(1-mh*mh*v_c*v_c) - r
    
    adda_cmdlineargs_copy = dict(adda_cmdlineargs)
    adda_cmdlineargs_copy["mhost"] = "1 0"
    adda_cmdlineargs_copy["size"] = 225
    adda_cmdlineargs_copy["beam"] = f"electron {KK}"
    adda_cmdlineargs_copy["beam_center"] = "150 0 0"
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh1.5_scaled")
    aw.spectrum_execute(aw_parameters,adda_cmdlineargs_copy)
    aw.spectrum_collect("Peels",aw_parameters)
    aw.spectrum_collect("Pcl",aw_parameters)
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh1.5_scaled")
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv", delimiter=',')[1:]
    ax1.plot(data[:,0], data[:,1]/mh, label=""+r"${\it m}$"+r"$_{\rm h}$=1.5 (scaled)", color="purple", linestyle="none", marker="o", fillstyle="none", markersize=5, alpha=alp)
    
    
    
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh2")
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv", delimiter=',')[1:]
    ax1.plot(data[:,0], data[:,1], label=""+r"${\it m}$"+r"$_{\rm h}$=2", color="red", linestyle="-", linewidth=lw, alpha=alp)
    
    aw_parameters["mp_files"] = [os.path.abspath(__file__ + "/../" + "Ag_Palik_Garcia_div2.csv")]
    mh = 2
    v_c = math.sqrt(1-(r/(K+r))**2)
    beta_eps = v_c*mh
    
    adda_cmdlineargs_copy = dict(adda_cmdlineargs)
    adda_cmdlineargs_copy["mhost"] = "1 0"
    adda_cmdlineargs_copy["size"] = 300
    adda_cmdlineargs_copy["beam"] = "electron 999999"
    adda_cmdlineargs_copy["beam_center"] = "200 0 0"
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh2_scaled")
    aw.spectrum_execute(aw_parameters,adda_cmdlineargs_copy)
    aw.spectrum_collect("Peels",aw_parameters)
    aw.spectrum_collect("Pcl",aw_parameters)
    aw.spectrum_collect("Pcltot",aw_parameters)
    
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv", delimiter=',')[1:]
    ax1.plot(data[:,0], data[:,1]/(mh*beta_eps**2), label=""+r"${\it m}$"+r"$_{\rm h}$=2 (scaled)", color="red", linestyle="none", marker="o", fillstyle="none", markersize=5, alpha=alp)
    
    
    
    # aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh3")
    # data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv", delimiter=',')[1:]
    # ax1.plot(data[:,0], data[:,1], label=""+r"${\it m}$"+r"$_{\rm h}$=3", color="orange", linestyle="-", linewidth=lw, alpha=alp)
    
    # aw_parameters["mp_files"] = [os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Ag_Palik_Garcia_div3.csv")]
    # mh = 3
    # v_c = math.sqrt(1-(r/(K+r))**2)
    # beta_eps = v_c*mh
    
    # adda_cmdlineargs_copy = dict(adda_cmdlineargs)
    # adda_cmdlineargs_copy["mhost"] = "1 0"
    # adda_cmdlineargs_copy["size"] = 450
    # adda_cmdlineargs_copy["beam"] = "electron 999999"
    # adda_cmdlineargs_copy["beam_center"] = "300 0 0"
    
    # aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh3_scaled")
    # #aw.spectrum_execute(aw_parameters,adda_cmdlineargs_copy)
    # aw.spectrum_collect("Peels",aw_parameters)
    # aw.spectrum_collect("Pcl",aw_parameters)
    # aw.spectrum_collect("Pcltot",aw_parameters)
    
    # data = np.genfromtxt(f"{aw_parameters['dirname']}/Peels.csv", delimiter=',')[1:]
    # ax1.plot(data[:,0], data[:,1]/(mh*beta_eps**2), label=""+r"${\it m}$"+r"$_{\rm h}$=3 (scaled)", color="orange", linestyle="none", marker="o", fillstyle="none", markersize=5, alpha=alp)
    
    
    ax1.grid(which="both", axis="x", linestyle="None")
    ax1.set_xlim([min(data[:,0]),max(data[:,0])])
    ax1.set_ylabel("EELS probability, eV$^{-1}$")
    #ax.grid(False)
    if min(data[:,1]) > 0:
        ax1.set_ylim(bottom = 0)
    ax1.legend(loc="lower right", ncol=3, bbox_to_anchor=(0., 1.02, 1., .102), borderaxespad=0.)
    aw.plot_setaspect(ax1)
    fig1.savefig("Plot_EELS.svg", bbox_inches='tight')
    
    
    
    fig2,ax2 = aw.plot_create()
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh1")
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Pcl.csv", delimiter=',')[1:]
    ax2.plot(data[:,0], data[:,1], label=""+r"${\it m}$"+r"$_{\rm h}$=1", color=aw.color_for_plot("Pcl"), linewidth=lw, alpha=alp)
    
    ax2.plot(np.NaN, np.NaN, color='none', label=' ')
    
    dirname = os.path.abspath(__file__ + "/../" + "spectrum_mh1.5")
    data = np.genfromtxt(f"{dirname}/Pcl.csv", delimiter=',')[1:]
    ax2.plot(data[:,0], data[:,1], label=""+r"${\it m}$"+r"$_{\rm h}$=1.5", color="purple", linestyle="-", linewidth=lw, alpha=alp)
    
    
    mh=1.5
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh1.5_scaled")
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Pcl.csv", delimiter=',')[1:]
    ax2.plot(data[:,0], data[:,1]/mh, label=""+r"${\it m}$"+r"$_{\rm h}$=1.5 (scaled)", color="purple", linestyle="none", marker="o", fillstyle="none", markersize=5, alpha=alp)
    
    dirname = os.path.abspath(__file__ + "/../" + "spectrum_mh2")
    data = np.genfromtxt(f"{dirname}/Pcl.csv", delimiter=',')[1:]
    ax2.plot(data[:,0], data[:,1], label=""+r"${\it m}$"+r"$_{\rm h}$=2", color="red", linestyle="-", linewidth=lw, alpha=alp)
    
    mh=2
    v_c = math.sqrt(1-(r/(K+r))**2)
    beta_eps = v_c*mh
    
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh2_scaled")
    data = np.genfromtxt(f"{aw_parameters['dirname']}/Pcl.csv", delimiter=',')[1:]
    ax2.plot(data[:,0], data[:,1]/(mh*beta_eps**2), label=""+r"${\it m}$"+r"$_{\rm h}$=2 (scaled)", color="red", linestyle="none", marker="o", fillstyle="none", markersize=5, alpha=alp)
    
    
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_mh3")
    # data = np.genfromtxt(f"{dirname}/Pcl.csv", delimiter=',')[1:]
    # ax2.plot(data[:,0], data[:,1], label=""+r"${\it m}$"+r"$_{\rm h}$=3", color="orange", linestyle="-", linewidth=lw, alpha=alp)
    
    # mh=3
    # v_c = math.sqrt(1-(r/(K+r))**2)
    # beta_eps = v_c*mh
    
    # aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum_mh3_scaled")
    # data = np.genfromtxt(f"{aw_parameters['dirname']}/Pcl.csv", delimiter=',')[1:]
    # ax2.plot(data[:,0], data[:,1]/(beta_eps**2), label=""+r"${\it m}$"+r"$_{\rm h}$=3 (scaled)", color="orange", linestyle="none", marker="o", fillstyle="none", markersize=5, alpha=alp)
    
    
    
    ax2.grid(which="both", axis="x", linestyle="None")
    ax2.set_xlim([min(data[:,0]),max(data[:,0])])
    ax2.set_ylabel("CL probability, eV$^{-1}$")
    #ax2.grid(False)
    if min(data[:,1]) > 0:
        ax2.set_ylim(bottom = 0)
    ax2.legend(loc="lower right", ncol=3, bbox_to_anchor=(0., 1.02, 1., .102), borderaxespad=0.)
    aw.plot_setaspect(ax2)
    fig2.savefig("Plot_CL.svg", bbox_inches='tight')
