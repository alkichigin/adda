import sys, os, multiprocessing
if sys.path[1] != os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper") : sys.path.insert(1,os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper")) #This is to import ADDA Wrapper from parent directory
import ADDAwrapper as aw

mp_fileGr = []
mp_fileGr.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "graphene_nm_Song.csv"))

mp_fileAu = []
mp_fileAu.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Au_nm_McPeak.csv"))

mp_fileAuSur = []
mp_fileAuSur.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Au_nm_McPeak.csv"))
mp_fileAuSur.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "const_nm_1.5.csv"))

mp_fileGrSur = []
mp_fileGrSur.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "graphene_nm_Song.csv"))
mp_fileGrSur.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "const_nm_1.5.csv"))

mp_fileGrAu = []
mp_fileGrAu.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "graphene_nm_Song.csv"))
mp_fileGrAu.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Au_nm_McPeak.csv"))

mp_fileGrAuSur = []
mp_fileGrAuSur.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "graphene_nm_Song.csv"))
mp_fileGrAuSur.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "Au_nm_McPeak.csv"))
mp_fileGrAuSur.append(os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/refractive_index/" + "const_nm_1.5.csv"))

# ADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = aw.addaexec_find(mode="seq"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-4, #number of parallel processes is equal to the number of processor cores minus 1

    mp_file = mp_fileGr, #files with refractive index of the particle, each string contains: ev,mp_re,mp_im
    
    nm_range = range(300,1500+1,10),
    
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

d = 1 #[nm] - dipole size
l1 = 50
l2 = 100
l3 = 140

h1 = 15
h2 = 25
h3 = 65

# ADDA command line arguments
adda_cmdlineargs = dict(
    # Particle
    shape = f"cylinder .2",
    size = 50, #[nm]
    grid = 8, #dipoles per x-axis size of the particle
    mhost = "1.5 0",
    
    # Beam
    beam = "plane", #Energy[keV]
    #beam_center = "25 0 0", # x[nm] y[nm] z[nm] - beam center coordinates 
    prop = "0 0 1", #beam propagation direction vector
    
    # Precision and performance
    eps = 2, #Residual norm
    
    # Additional options
    sym = "enf", #Do not simulate second polarization
    scat_matr = "none", #Do not calculate the Mueller matrix
    no_vol_cor = "", #Disable volume correction
    iter = "qmr2", #Iterative solver
    pol = "igt_so", #Polarizability prescription
    int = "igt_so", #Interaction term
    # Csca = "", #Calculate Csca with the Romberg integral. Needed to properly calculate Cathodoluminesce
    # alldir_inp = os.path.abspath(__file__ + "/../../../../misc/ADDAwrapper/Csca_integration.txt")
)

### Executing commands
if __name__ == '__main__': 
    
    import csv,os
    import numpy as np
    
    os.makedirs("geom",exist_ok=True)
    aw_parameters["dirname"] = os.path.abspath(__file__ + "/../")
    
    # with open("geom/Sur_for_Au.txt", "w") as file:
    #     file.write("0 0\n")
    #     file.write(f"{h3/2} 0\n")
    #     file.write(f"{h3/2} {l3}\n")
    #     file.write(f"0 {l3}\n")
    #     file.write(f"0 {l3 - (l3-l1)/2}\n")
    #     file.write(f"{h1/2} {l3 - (l3-l1)/2}\n")
    #     file.write(f"{h1/2} {(l3-l1)/2}\n")
    #     file.write(f"0 {(l3-l1)/2}\n")
    # adda_cmdlineargs["shape"] = "axisymmetric geom/Sur_for_Au.txt"
    # adda_cmdlineargs["grid"] = round(h3/d)
    # grid = adda_cmdlineargs["grid"]
    # aw.geom_gen(aw_parameters,adda_cmdlineargs)
    # shape1 = adda_cmdlineargs["shape"].split()[0]
    # Sur_for_Au = np.genfromtxt(f"geom/{shape1}.geom", delimiter=' ', dtype="int")
    # Sur_for_Au[:,0], Sur_for_Au[:,1], Sur_for_Au[:,2] = Sur_for_Au[:,0], Sur_for_Au[:,2], -Sur_for_Au[:,1] #Mx(90)
    # Sur_for_Au[:,0], Sur_for_Au[:,1], Sur_for_Au[:,2] = Sur_for_Au[:,0]-min(Sur_for_Au[:,0]), Sur_for_Au[:,1]-min(Sur_for_Au[:,1]), Sur_for_Au[:,2]-min(Sur_for_Au[:,2])
    # with open("geom/Sur_for_Au.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerows(zip(Sur_for_Au[:,0],Sur_for_Au[:,1],Sur_for_Au[:,2]))
        
    # with open("geom/Sur_for_GrAu.txt", "w") as file:
    #     file.write("0 0\n")
    #     file.write(f"{h3/2} 0\n")
    #     file.write(f"{h3/2} {l3}\n")
    #     file.write(f"0 {l3}\n")
    #     file.write(f"0 {l3 - (l3-l1)/2}\n")
    #     file.write(f"{h1/2} {l3 - (l3-l1)/2}\n")
    #     file.write(f"{h1/2} {l3 - (l3-l2)/2}\n")
    #     file.write(f"{h2/2} {l3 - (l3-l2)/2}\n")
    #     file.write(f"{h2/2} {(l3-l2)/2}\n")
    #     file.write(f"{h1/2} {(l3-l2)/2}\n")
    #     file.write(f"{h1/2} {(l3-l1)/2}\n")
    #     file.write(f"0 {(l3-l1)/2}\n")
    # adda_cmdlineargs["shape"] = "axisymmetric geom/Sur_for_GrAu.txt"
    # adda_cmdlineargs["grid"] = round(h3/d)
    # grid = adda_cmdlineargs["grid"]
    # aw.geom_gen(aw_parameters,adda_cmdlineargs)
    # shape1 = adda_cmdlineargs["shape"].split()[0]
    # Sur_for_GrAu = np.genfromtxt(f"geom/{shape1}.geom", delimiter=' ', dtype="int")
    # Sur_for_GrAu[:,0], Sur_for_GrAu[:,1], Sur_for_GrAu[:,2] = Sur_for_GrAu[:,0], Sur_for_GrAu[:,2], -Sur_for_GrAu[:,1] #Mx(90)
    # Sur_for_GrAu[:,0], Sur_for_GrAu[:,1], Sur_for_GrAu[:,2] = Sur_for_GrAu[:,0]-min(Sur_for_GrAu[:,0]), Sur_for_GrAu[:,1]-min(Sur_for_GrAu[:,1]), Sur_for_GrAu[:,2]-min(Sur_for_GrAu[:,2])
    # with open("geom/Sur_for_GrAu.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerows(zip(Sur_for_GrAu[:,0],Sur_for_GrAu[:,1],Sur_for_GrAu[:,2]))
        

    # with open("geom/Sur_for_Gr.txt", "w") as file:
    #     file.write("0 0\n")
    #     file.write(f"{h3/2} 0\n")
    #     file.write(f"{h3/2} {l3}\n")
    #     file.write(f"{h1/2} {l3}\n")
    #     file.write(f"{h1/2} {l3 - (l3-l2)/2}\n")
    #     file.write(f"{h2/2} {l3 - (l3-l2)/2}\n")
    #     file.write(f"{h2/2} {(l3-l2)/2}\n")
    #     file.write(f"{h1/2} {(l3-l2)/2}\n")
    #     file.write(f"{h1/2} {l3}\n")
    #     file.write(f"0 {l3}\n")
    # adda_cmdlineargs["shape"] = "axisymmetric geom/Sur_for_Gr.txt"
    # adda_cmdlineargs["grid"] = round(h3/d)
    # grid = adda_cmdlineargs["grid"]
    # aw.geom_gen(aw_parameters,adda_cmdlineargs)
    # shape1 = adda_cmdlineargs["shape"].split()[0]
    # Sur_for_Gr = np.genfromtxt(f"geom/{shape1}.geom", delimiter=' ', dtype="int")
    # Sur_for_Gr[:,0], Sur_for_Gr[:,1], Sur_for_Gr[:,2] = Sur_for_Gr[:,0], Sur_for_Gr[:,2], -Sur_for_Gr[:,1] #Mx(90)
    # Sur_for_Gr[:,0], Sur_for_Gr[:,1], Sur_for_Gr[:,2] = Sur_for_Gr[:,0]-min(Sur_for_Gr[:,0]), Sur_for_Gr[:,1]-min(Sur_for_Gr[:,1]), Sur_for_Gr[:,2]-min(Sur_for_Gr[:,2])
    # with open("geom/Sur_for_Gr.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerows(zip(Sur_for_Gr[:,0],Sur_for_Gr[:,1],Sur_for_Gr[:,2]))

    
    # with open("geom/Gr.txt", "w") as file:
    #     file.write(f"{h1/2} 0\n")
    #     file.write(f"{h2/2} 0\n")
    #     file.write(f"{h2/2} {l2}\n")
    #     file.write(f"{h1/2} {l2}\n")
    # adda_cmdlineargs["shape"] = "axisymmetric geom/Gr.txt"
    # adda_cmdlineargs["grid"] = round(h2/d)
    # grid = adda_cmdlineargs["grid"]
    # aw.geom_gen(aw_parameters,adda_cmdlineargs)
    # shape1 = adda_cmdlineargs["shape"].split()[0]
    # dataGr = np.genfromtxt(f"geom/{shape1}.geom", delimiter=' ', dtype="int")
    # dataGr[:,0], dataGr[:,1], dataGr[:,2] = dataGr[:,0], dataGr[:,2], -dataGr[:,1] #Mx(90)
    # dataGr[:,0], dataGr[:,1], dataGr[:,2] = dataGr[:,0]-min(dataGr[:,0]), dataGr[:,1]-min(dataGr[:,1]), dataGr[:,2]-min(dataGr[:,2])
    # dataGr[:,0] = dataGr[:,0] + int((max(Sur_for_GrAu[:,0])-max(dataGr[:,0])+1)/2)
    # dataGr[:,1] = dataGr[:,1] + int((max(Sur_for_GrAu[:,1])-max(dataGr[:,1])+1)/2)
    # dataGr[:,2] = dataGr[:,2] + int((max(Sur_for_GrAu[:,2])-max(dataGr[:,2])+1)/2)
    # msGr = np.full(shape=len(dataGr[:,0]),fill_value=1,dtype="int")
    # with open("geom/Gr.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerows(zip(dataGr[:,0],dataGr[:,1],dataGr[:,2]))
    
    
    # with open("geom/Au.txt", "w") as file:
    #     file.write("0 0\n")
    #     file.write(f"{h1/2} 0\n")
    #     file.write(f"{h1/2} {l1}\n")
    #     file.write(f"0 {l1}\n")
    # adda_cmdlineargs["shape"] = "axisymmetric geom/Au.txt"
    # adda_cmdlineargs["grid"] = round(h1/d)
    # grid = adda_cmdlineargs["grid"]
    # aw.geom_gen(aw_parameters,adda_cmdlineargs)
    # shape1 = adda_cmdlineargs["shape"].split()[0]
    # dataAu = np.genfromtxt(f"geom/{shape1}.geom", delimiter=' ', dtype="int")
    # dataAu[:,0], dataAu[:,1], dataAu[:,2] = dataAu[:,0], dataAu[:,2], -dataAu[:,1] #Mx(90)
    # dataAu[:,0], dataAu[:,1], dataAu[:,2] = dataAu[:,0]-min(dataAu[:,0]), dataAu[:,1]-min(dataAu[:,1]), dataAu[:,2]-min(dataAu[:,2])
    # dataAu[:,0] = dataAu[:,0] + int((max(Sur_for_GrAu[:,0])-max(dataAu[:,0])+1)/2)
    # dataAu[:,1] = dataAu[:,1] + int((max(Sur_for_GrAu[:,1])-max(dataAu[:,1])+1)/2)
    # dataAu[:,2] = dataAu[:,2] + int((max(Sur_for_GrAu[:,2])-max(dataAu[:,2])+1)/2)
    # msAu = np.full(shape=len(dataAu[:,0]),fill_value=2,dtype="int")
    # with open("geom/Au.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerows(zip(dataAu[:,0],dataAu[:,1],dataAu[:,2]))
        
    # data, ms = np.concatenate((dataGr,dataAu)), np.concatenate((msGr,msAu), axis=None)
    # with open("geom/GrAu.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerow(["Nmat=2"])
    #     writer.writerows(zip(data[:,0],data[:,1],data[:,2],ms))
        
    # msSur_for_Gr = np.full(shape=len(Sur_for_Gr[:,0]),fill_value=2,dtype="int")
    # data, ms = np.concatenate((dataGr,Sur_for_Gr)), np.concatenate((msGr,msSur_for_Gr), axis=None)
    # with open("geom/GrSur.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerow(["Nmat=2"])
    #     writer.writerows(zip(data[:,0],data[:,1],data[:,2],ms))
        
    # msSur_for_GrAu = np.full(shape=len(Sur_for_GrAu[:,0]),fill_value=3,dtype="int")
    # data, ms = np.concatenate((dataGr,dataAu,Sur_for_GrAu)), np.concatenate((msGr,msAu,msSur_for_GrAu), axis=None)
    # with open("geom/GrAuSur.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerow(["Nmat=3"])
    #     writer.writerows(zip(data[:,0],data[:,1],data[:,2],ms))
        
    # msSur_for_Au = np.full(shape=len(Sur_for_Au[:,0]),fill_value=2,dtype="int")
    # msAu = np.full(shape=len(dataAu[:,0]),fill_value=1,dtype="int")
    # data, ms = np.concatenate((dataAu,Sur_for_Au)), np.concatenate((msAu,msSur_for_Au), axis=None)
    # with open("geom/AuSur.geom", 'w') as file:
    #     writer = csv.writer(file, delimiter=' ', lineterminator='\n')
    #     writer.writerow(["Nmat=2"])
    #     writer.writerows(zip(data[:,0],data[:,1],data[:,2],ms))



    os.makedirs("visual",exist_ok=True)
    
    # fig1,ax1 = aw.geometry("geom/Gr.geom",["dimgrey"], half=True)
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/Gr.png", bbox_inches='tight', dpi=600)
    # ax1.axis('off')
    # fig1.savefig("visual/Gr_axoff.png", bbox_inches='tight', dpi=600)
    
    # fig1,ax1 = aw.geometry("geom/GrAu.geom",["dimgrey","gold"], half=True)
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/GrAu.png", bbox_inches='tight', dpi=600)
    # ax1.axis('off')
    # fig1.savefig("visual/GrAu_axoff.png", bbox_inches='tight', dpi=600)
    
    # fig1,ax1 = aw.geometry("geom/Sur_for_GrAu.geom",["dimgrey"], half=True)
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/Sur_for_GrAu.png", bbox_inches='tight', dpi=600)
    # ax1.axis('off')
    # fig1.savefig("visual/Sur_for_GrAu_axoff.png", bbox_inches='tight', dpi=600)
    
    # fig1,ax1 = aw.geometry("geom/Sur_for_Gr.geom",["dimgrey"], half=True)
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/Sur_for_Gr.png", bbox_inches='tight', dpi=600)
    # ax1.axis('off')
    # fig1.savefig("visual/Sur_for_Gr_axoff.png", bbox_inches='tight', dpi=600)
    
    # fig1,ax1 = aw.geometry("geom/GrSur.geom",["dimgrey","silver"], half=True)
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/GrSur.png", bbox_inches='tight', dpi=600)
    # ax1.axis('off')
    # fig1.savefig("visual/GrSur_axoff.png", bbox_inches='tight', dpi=600)
    
    # fig1,ax1 = aw.geometry("geom/GrAuSur.geom",["dimgrey","gold","silver"], half=True)
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/GrAuSur.png", bbox_inches='tight', dpi=600)
    # ax1.axis('off')
    # fig1.savefig("visual/GrAuSur_axoff.png", bbox_inches='tight', dpi=600)
    
    # fig1,ax1 = aw.geometry("geom/AuSur.geom",["gold","silver"], half=True)
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/AuSur.png", bbox_inches='tight', dpi=600)
    # ax1.axis('off')
    # fig1.savefig("visual/AuSur_axoff.png", bbox_inches='tight', dpi=600)
    
    # fig1,ax1 = aw.geometry("geom/Sur_for_Au.geom",["silver"], half=True)
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/Sur_for_Au.png", bbox_inches='tight', dpi=600)
    # ax1.axis('off')
    # fig1.savefig("visual/Sur_for_Au_axoff.png", bbox_inches='tight', dpi=600)
    
    # fig1,ax1 = aw.geometry("geom/Au.geom",["gold","silver"], half=True)
    # ax1.view_init(30,-30)
    # fig1.savefig("visual/Au.png", bbox_inches='tight', dpi=600)
    # ax1.axis('off')
    # fig1.savefig("visual/Au_axoff.png", bbox_inches='tight', dpi=600)
    
    
    # aw_parameters["mp_files"] = mp_fileAu
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_Au")
    # adda_cmdlineargs["shape"] = "read geom/Au.geom"
    # adda_cmdlineargs["size"] = l1
    # adda_cmdlineargs["grid"] = round(adda_cmdlineargs["size"]/d)
    # #aw.spectrum_nm_execute(aw_parameters,adda_cmdlineargs,dirname)
    # aw.spectrum_collect("Cext",dirname)
    # aw.spectrum_nm_plot("Cext",dirname)
    
    # aw_parameters["mp_files"] = mp_fileAuSur
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_AuSur")
    # adda_cmdlineargs["shape"] = "read geom/AuSur.geom"
    # adda_cmdlineargs["size"] = l3
    # adda_cmdlineargs["grid"] = round(adda_cmdlineargs["size"]/d)
    # adda_cmdlineargs["store_int_field"] = ""
    # #aw.spectrum_nm_execute(aw_parameters,adda_cmdlineargs,dirname)
    # aw.spectrum_collect("Cext",dirname)
    # aw.spectrum_nm_plot("Cext",dirname)
    # aw.delkeys_silent(adda_cmdlineargs, ["store_int_field"])
    
    
    # aw_parameters["mp_file"] = mp_fileGr
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_Gr")
    # adda_cmdlineargs["shape"] = "read geom/Gr.geom"
    # adda_cmdlineargs["size"] = l2
    # adda_cmdlineargs["grid"] = round(adda_cmdlineargs["size"]/d)
    # #aw.spectrum_nm_execute(aw_parameters,adda_cmdlineargs,dirname)
    # aw.spectrum_collect("Cext",dirname)
    # aw.spectrum_nm_plot("Cext",dirname)
    
    # aw_parameters["mp_file"] = mp_fileGrSur
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_GrSur")
    # adda_cmdlineargs["shape"] = "read geom/GrSur.geom"
    # adda_cmdlineargs["size"] = l3
    # adda_cmdlineargs["grid"] = round(adda_cmdlineargs["size"]/d)
    # adda_cmdlineargs["store_int_field"] = ""
    # #aw.spectrum_nm_execute(aw_parameters,adda_cmdlineargs,dirname)
    # aw.spectrum_collect("Cext",dirname)
    # aw.spectrum_nm_plot("Cext",dirname)
    # aw.delkeys_silent(adda_cmdlineargs, ["store_int_field"])
    

    # aw_parameters["mp_file"] = mp_fileGrAu
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_GrAu")
    # adda_cmdlineargs["shape"] = "read geom/GrAu.geom"
    # adda_cmdlineargs["size"] = l2
    # adda_cmdlineargs["grid"] = round(adda_cmdlineargs["size"]/d)
    # #aw.spectrum_nm_execute(aw_parameters,adda_cmdlineargs,dirname)
    # aw.spectrum_collect("Cext",dirname)
    # aw.spectrum_nm_plot("Cext",dirname)
    
    # aw_parameters["mp_file"] = mp_fileGrAuSur
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_GrAuSur")
    # adda_cmdlineargs["shape"] = "read geom/GrAuSur.geom"
    # adda_cmdlineargs["size"] = l3
    # adda_cmdlineargs["grid"] = round(adda_cmdlineargs["size"]/d)
    # adda_cmdlineargs["store_int_field"] = ""
    # #aw.spectrum_nm_execute(aw_parameters,adda_cmdlineargs,dirname)
    # aw.spectrum_collect("Cext",dirname)
    # aw.spectrum_nm_plot("Cext",dirname)
    # aw.delkeys_silent(adda_cmdlineargs, ["store_int_field"])
    
    
    # fig,ax = aw.plot_create()
    
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_GrSur")
    # data = np.genfromtxt(f"{dirname}/Cext.csv", delimiter=',')[1:]
    # ax.plot(data[:,0], data[:,1], label="Graphene", color="black", linewidth=3)
    
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_GrAuSur")
    # data = np.genfromtxt(f"{dirname}/Cext.csv", delimiter=',')[1:]
    # ax.plot(data[:,0], data[:,1], label="Graphene/Au", color="red", linewidth=3)
    
    # #ax.grid(which="both", axis="x", linestyle="None")
    # ax.set_xlim([min(data[:,0]),max(data[:,0])])
    # ax.set_ylim([1*1e3,5*1e3])
    # ax.set_xlabel("Wavelength, nm")
    # ax.set_ylabel(r"C$_{\rm ext}$, nm$^2$")
    # #ax.grid(which="both", axis="x", linestyle="dotted", zorder=0)
    # #ax.grid(False)
    # ax.legend()
    # #ax.legend(loc="upper center")
    # aw.plot_setaspect(ax)
    # fig.savefig("spectra.png", bbox_inches='tight', dpi=600)
    
    
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_GrAuSur")
    # data = np.genfromtxt(f"{dirname}/ADDA_output/1100.0/IntField-Y", delimiter=' ', dtype="float", skip_header=1)
    # data = np.array(data[:,0:4])
    # data = data[data[:,0] == -37.5]
    # with open(f"{dirname}/Field1100GrAu.csv", 'w') as file:
    #     writer = csv.writer(file, delimiter=',', lineterminator='\n')
    #     writer.writerow(["y","z","|E|^2"])
    #     writer.writerows(zip(data[:,1],data[:,2],data[:,3]))
    # aw.scan_plot("Field1100GrAu",dirname)
    
    # dirname = os.path.abspath(__file__ + "/../" + "spectrum_GrSur")
    # data = np.genfromtxt(f"{dirname}/ADDA_output/1100.0/IntField-Y", delimiter=' ', dtype="float", skip_header=1)
    # data = np.array(data[:,0:4])
    # data = data[data[:,0] == -37.5]
    # with open(f"{dirname}/Field1100Gr.csv", 'w') as file:
    #     writer = csv.writer(file, delimiter=',', lineterminator='\n')
    #     writer.writerow(["y","z","|E|^2"])
    #     writer.writerows(zip(data[:,1],data[:,2],data[:,3]))
    # aw.scan_plot("Field1100Gr",dirname)
    
    dirname = os.path.abspath(__file__ + "/../" + "spectrum_AuSur")
    data = np.genfromtxt(f"{dirname}/ADDA_output/1100.0/IntField-Y", delimiter=' ', dtype="float", skip_header=1)
    data = np.array(data[:,0:4])
    data = data[data[:,0] == -37.5]
    with open(f"{dirname}/Field1100Au.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["y","z","|E|^2"])
        writer.writerows(zip(data[:,1],data[:,2],data[:,3]))
    aw.scan_plot("Field1100Au",dirname)
    
    dirname = os.path.abspath(__file__ + "/../" + "spectrum_AuSur")
    data = np.genfromtxt(f"{dirname}/ADDA_output/830.0/IntField-Y", delimiter=' ', dtype="float", skip_header=1)
    data = np.array(data[:,0:4])
    data = data[data[:,0] == -37.5]
    with open(f"{dirname}/Field830Au.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["y","z","|E|^2"])
        writer.writerows(zip(data[:,1],data[:,2],data[:,3]))
    aw.scan_plot("Field830Au",dirname)

