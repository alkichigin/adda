import sys, os, multiprocessing
if sys.path[1] != os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper") : sys.path.insert(1,os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper")) #This is to import ADDA Wrapper from parent directory
import ADDAwrapper as aw

mp_files = [os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Ag_nm_McPeak.csv")]

# ADDAWrapper parameters
aw_parameters = dict(
    adda_exec = aw.addaexec_find(mode="seq"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-1, #number of parallel processes is equal to the number of processor cores minus 1

    mp_files = mp_files, #files with refractive index of the particle, each string contains: ev,mp_re,mp_im
    
    nm_range = range(300,1000+1,10),
    
    # spectrumline_startpoint = (0,200), # (x,y) [nm]
    # spectrumline_endpoint = (25,220), # (x,y) [nm]
    # spectrumline_points = 10, #how many points, including startpoint and endpoint
    
    #Used in "scan_" functions. Beam propagation must be orthogonal to the grid.
    #So "prop" must be "0 0 whatever" and rotations with "orient" must be made by 90 degrees.
    
    #The beam must always blast exactly in the middle between the dipoles,
    #so start and stop coordinates will be adjusted, covering more area than you entered. Obligatory to use "no_vol_cor" with scan.
)

# ADDA command line arguments
adda_cmdlineargs = dict(
    # Particle
    shape = f"chebyshev 1 15",
    size = 100, #[nm]
    grid = 21, #dipoles per x-axis size of the particle
    #surf = "17.5 1.5 0",
    
    # Beam
    beam = "plane", #plane wave for today
    #beam_center = "32.5 70 0", # x[nm] y[nm] z[nm] - beam center coordinates 
    prop = "0 0 -1", #beam propagation direction vector
    
    # Precision and performance
    eps = 2, #Residual norm
    
    # Additional options
    sym = "enf", #Do not simulate second polarization
    scat_matr = "none", #Do not calculate the Mueller matrix
    no_vol_cor = "", #Disable volume correction
    iter = "qmr", #Iterative solver
    pol = "igt_so", #Polarizability prescription
    int = "igt_so", #Interaction term
    # Csca = "", #Calculate Csca with the Romberg integral. Needed to properly calculate Cathodoluminesce
    # alldir_inp = os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/Csca_default.txt"),
)

### Executing commands
if __name__ == '__main__': 
    
    import csv,os,math
    import numpy as np
    
    os.makedirs("geom",exist_ok=True)
    os.makedirs("visual",exist_ok=True)

    
    grid = 61 # Size of the grid, chosen to be odd for symmetry
    # lower grid - roughly discretized particle, but faster simulations
    # higher grid - prettier discretization, but slowing down simulations
    
    ''' Making 3D box of (grid x grid x grid) cubical dipoles '''
    xs0, ys0, zs0 = [], [], [] #making 3D grid
    for i in range(grid):
        for j in range(grid):
            for k in range(grid):
                xs0.append(k),ys0.append(j),zs0.append(i)
    xs0 = np.asarray(xs0)
    ys0 = np.asarray(ys0)
    zs0 = np.asarray(zs0)
    xs0 = xs0 - int(max(xs0+1)/2)
    ys0 = ys0 - int(max(ys0+1)/2)
    zs0 = zs0 - int(max(zs0+1)/2)
    data = np.array([xs0,ys0,zs0]).T
    with open("geom/box.geom", 'w') as file:
        writer = csv.writer(file, delimiter=' ', lineterminator='\n')
        writer.writerows(data)
    
    # fig1,ax1 = aw.geometry("geom/box.geom",["silver"])
    # ax1.view_init(30,-60)
    # fig1.savefig("visual/box.png", bbox_inches='tight', dpi=600)
    # # ax1.axis('off')
    # # fig1.savefig("visual/box_axoff.png", bbox_inches='tight', dpi=600)

    
    
    '''Making a spiked star-shaped structure using modified Chebyshev particle'''
    r0 = grid/3.5 # r0 should be calculated precisely to fit the grid, but we are using some approximate value
    e = 1 # proportional to the amplitude of spikes
    n = 6 # proportional to the number of vertical spikes
    m = n # proportional to the number of horizontal spikes
    spsh = 19/3 # "spike sharpening", should be an odd number, makes spikes sharper
    data = np.genfromtxt("geom/box.geom", delimiter=' ', dtype="int")
    theta = np.arctan(np.sqrt(data[:,0]**2 + data[:,1]**2)/data[:,2])
    phi = np.arctan(data[:,1]/data[:,0]) 
    data = data[np.sqrt(data[:,0]**2 + data[:,1]**2 + data[:,2]**2) < r0*(1+e*(np.cos(n*theta)*np.cos(m*phi))**spsh)] # Chebyshev shape: r=r0(1+εcos(nθ)cos(mφ))
    
    with open("geom/chebyshev.geom", 'w') as file:
        writer = csv.writer(file, delimiter=' ', lineterminator='\n')
        writer.writerows(data)
    
    # fig1,ax1 = aw.geometry("geom/chebyshev.geom",["silver"])
    # ax1.view_init(30,-60)
    # fig1.savefig("visual/chebyshev.png", bbox_inches='tight', dpi=600)
    # # ax1.axis('off')
    # # fig1.savefig("visual/chebyshev_axoff.png", bbox_inches='tight', dpi=600)
    
    
    
    '''Making a sphere to place inside chebyshev particle'''
    R = grid/4 #sphere radius
    data = np.genfromtxt("geom/box.geom", delimiter=' ', dtype="int")
    data = data[data[:,0]**2 + data[:,1]**2 + data[:,2]**2 < R**2] #sphere parametric equation
    with open("geom/sphere.geom", 'w') as file:
        writer = csv.writer(file, delimiter=' ', lineterminator='\n')
        writer.writerows(data)
        
    # fig1,ax1 = aw.geometry("geom/sphere.geom",["gold"])
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/sphere.png", bbox_inches='tight', dpi=600)
    # # ax1.axis('off')
    # # fig1.savefig("visual/sphere_axoff.png", bbox_inches='tight', dpi=600)
    
    
    '''Excluding the same sphere from the chebyshev particle where they intersect'''
    data = np.genfromtxt("geom/chebyshev.geom", delimiter=" ", dtype=int)
    data = data[data[:,0]**2 + data[:,1]**2 + data[:,2]**2 >= R**2] #sphere parametric equation
    with open("geom/chebyshev_exc.geom", 'w') as file:
        writer = csv.writer(file, delimiter=' ', lineterminator='\n')
        writer.writerows(data)
    
    # fig1,ax1 = aw.geometry("geom/chebyshev_exc.geom",["silver"])
    # ax1.view_init(30,-60)
    # fig1.savefig("visual/chebyshev_exc.png", bbox_inches='tight', dpi=600)
    # # ax1.axis('off')
    # # fig1.savefig("visual/chebyshev_exc_axoff.png", bbox_inches='tight', dpi=600)
    
    
    
    '''Combining the chebyshev_exc and the sphere geometries as two different materials'''
    data = np.genfromtxt("geom/chebyshev_exc.geom", delimiter=" ", dtype=int)
    xs1, ys1, zs1 = data[:,0], data[:,1], data[:,2]
    ms1 = np.full(shape=len(xs1),fill_value=1,dtype="int") # material #1
    
    data = np.genfromtxt("geom/sphere.geom", delimiter=" ", dtype=int)
    xs2, ys2, zs2 = data[:,0], data[:,1], data[:,2]
    ms2 = np.full(shape=len(xs2),fill_value=2,dtype="int") # material #2
    
    xs = np.concatenate((xs1, xs2)) 
    ys = np.concatenate((ys1, ys2)) 
    zs = np.concatenate((zs1, zs2)) 
    ms = np.concatenate((ms1, ms2)) 
    
    with open("geom/nanoparticle.geom", 'w') as file:
        writer = csv.writer(file, delimiter=' ', lineterminator='\n')
        writer.writerow(["Nmat=2"]) # The number of materials used
        writer.writerows(zip(xs,ys,zs,ms))
    
    fig1,ax1 = aw.geometry("geom/nanoparticle.geom",["silver", "gold"],half=True)
    ax1.view_init(30,-60)
    fig1.savefig("visual/nanoparticle.png", bbox_inches='tight', dpi=600)
    ax1.axis('off')
    fig1.savefig("visual/nanoparticle_axoff.png", bbox_inches='tight', dpi=600)
    
    
    
    '''Simulating the absorption spectrum'''
    aw_parameters["mp_files"] = [os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Ag_nm_McPeak.csv"), # material #1
                                  os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Au_nm_McPeak.csv")] # material #2
    adda_cmdlineargs["shape"] = "read geom/nanoparticle.geom"
    grid_real = max(max(xs)-min(xs)+1, max(ys)-min(ys)+1, max(zs)-min(zs)+1) #since we did not calculate the r0 parameter exactly, we need to re-estimate the real grid
    
    adda_cmdlineargs["grid"] = grid_real 
    adda_cmdlineargs["size"] = 100 #nm
    
    aw_parameters["dirname"] = os.path.abspath(__file__ + "/../" + "spectrum")
    aw.spectrum_execute(aw_parameters,adda_cmdlineargs)
    aw.spectrum_collect("Cabs",aw_parameters)
    aw.spectrum_plot("Cabs",aw_parameters)








