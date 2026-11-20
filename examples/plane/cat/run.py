import sys, os, multiprocessing, numpy as np
if sys.path[1] != os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper") : sys.path.insert(1,os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper")) #This is to import ADDA Wrapper from parent directory
import ADDAwrapper as aw

# ADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = aw.addaexec_find(mode="seq"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-1, #number of parallel processes is equal to the number of processor cores minus 1

    mp_files = [os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Al_nm_McPeak.csv")], #files with refractive index of the particle, each string contains: ev,mp_re,mp_im
    
    nm_range = np.arange(100,800+1,10),
    
    spectrumline_startpoint = (0,30), # (x,y) [nm]
    spectrumline_endpoint = (0,30), # (x,y) [nm]
    spectrumline_points = 10, #how many points, including startpoint and endpoint
    
    #Used in "scan_" functions. Beam propagation must be orthogonal to the grid.
    #So "prop" must be "0 0 whatever" and rotations with "orient" must be made by 90 degrees.
    scan_x_range = (-30,30), #[nm], (x_left, x_right)
    scan_y_range = (-30,30), #[nm], (y_bottom, y_top)
    scan_step = 1 #dipoles per each step, must be an integer >= 1
    #The beam must always blast exactly in the middle between the dipoles,
    #so start and stop coordinates will be adjusted, covering more area than you entered. Obligatory to use "no_vol_cor" with scan.
)

# ADDA command line arguments
adda_cmdlineargs = dict(
    # Particle
    shape = f"read cat.geom",
    size = 199, #[nm]
    grid = 199, #dipoles per x-axis size of the particle
    
    # Beam
    beam = "plane", #Energy[keV]
    beam_center = "0 0 0", # x[nm] y[nm] z[nm] - beam center coordinates 
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
    alldir_inp = os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/Csca_default.txt"),
    store_int_field = ""
)

### Executing commands
if __name__ == '__main__': 
    
    import csv,os
    import numpy as np
    
    # os.makedirs("geom",exist_ok=True)
    # aw_parameters["dirname"] = os.path.abspath(__file__ + "/../")
    

    
    # aw_parameters['dirname'] = os.path.abspath(__file__ + "/../")
    # aw.geom_gen(aw_parameters,adda_cmdlineargs)
    # data = np.genfromtxt("cat.geom", delimiter=' ', dtype="int")
    # xs1, ys1, zs1 = data[:,0], data[:,1], data[:,2]
    # xs1, ys1, zs1 = xs1, -zs1, ys1 #Mx(90)
    # xs1, ys1, zs1 = ys1, -xs1, zs1 #Mz(90)
    # xs1 = xs1 - min(xs1)
    # ys1 = ys1 - min(ys1)
    # zs1 = zs1 - min(zs1)
    # with open("cat_rotated.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerows(zip(xs1,ys1,zs1))
    
    fig1,ax1 = aw.geometry("cat_rotated.geom", ["silver"])
    ax1.view_init(30,-30)
    fig1.savefig("cat_axis.png", bbox_inches='tight', dpi=1200)
    # ax1.axis('off')
    # fig1.savefig("cat.pdf", bbox_inches='tight')
    
    
    
    adda_cmdlineargs["shape"] = "read cat_rotated.geom"
    adda_cmdlineargs["size"] = 38
    adda_cmdlineargs["grid"] = 38
    aw_parameters['dirname'] = os.path.abspath(__file__ + "/../" + "spectrum")
    
    #aw.spectrum_execute(aw_parameters,adda_cmdlineargs)
    aw.spectrum_collect("Cext",aw_parameters)
    aw.spectrum_plot("Cext",aw_parameters)
    aw.spectrum_collect("Csca",aw_parameters)
    aw.spectrum_plot("Csca",aw_parameters)

    
    from matplotlib.colors import LogNorm
    lam = "400.0"
    
    aw.field_projection(aw_parameters,"y","z",lam,"Y")
    aw.scan_plot(f"field_projection_{lam}",aw_parameters,norm=LogNorm(vmin=0.01, vmax=1))
    
    aw.field_projection(aw_parameters,"x","z",lam,"Y")
    aw.scan_plot(f"field_projection_{lam}",aw_parameters,norm=LogNorm(vmin=0.01, vmax=1))
    
    aw.field_projection(aw_parameters,"y","x",lam,"Y")
    aw.scan_plot(f"field_projection_{lam}",aw_parameters,norm=LogNorm(vmin=0.01, vmax=1))



