import os, shutil, re, csv, time, multiprocessing, tqdm, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
 
def addaexec_find(mode="seq"):
    if mode!="seq" and mode!="mpi" and mode!="ocl":
        print(f"ERROR: unkwnown mode {mode}")
        return
    if mode == "seq": 
        name="adda"
    else: 
        name=f"adda_{mode}" 
    if os.system(name  + " > " + os.devnull + " 2>&1") == 0 : return name
    cmdline = os.path.abspath(__file__ + f"/../../../src/{mode}/{name}")
    if os.system(cmdline  + " -V > " + os.devnull + " 2>&1") == 0 : return cmdline
    cmdline = os.path.abspath(__file__ + f"/../../../win64/{name}")
    if os.name == "nt" and os.system(cmdline  + " -V > " + os.devnull + " 2>&1") == 0 : return cmdline
    print(f"ERROR: No working {name} binary found")

def label_for_plot(match):
    if match[0] == "P":
        if match[1:] == "ext": return r"$\it{P}$" + r"$_{\rm ext}$" + ", eV$^{-1}$"
        return r"$\it{" + match[0] + r"}$" + r"$_{\rm " + match[1:].upper() + "}$" + ", eV$^{-1}$"
    elif match[0] == "C":
        return r"$\it{" + match[0] + r"}$" + r"$_{\rm " + match[1:] + "}$" + ", nm$^2$"
    elif match[0] == "Q":
        return r"$\it{" + match[0] + r"}$" + r"$_{\rm " + match[1:] + "}$"
    else:
        return match

def label_for_plot_arbunits(match):
    if match[0] == "P":
        return match[0] + r"$_{\rm " + match[1:].upper() + "}$" + ", arb. units"
    elif match[0] == "C":
        return match + ", arb. units"
    elif match[0] == "Q":
        return match
    else:
        return match

def color_for_plot(match):
    if match == "Peels":
        return "royalblue"
    elif match == "Pcl":
        return "violet"
    else:
        return "black"

def print_log(string, dirname=False, silent=False):
    if silent == False:
        print(string)
    if dirname != False:
        with open(dirname + "/log.txt", 'a') as file:
            file.write(string + "\n")

def delkeys_silent(dictionary, keys):
    for key in keys:
        if key in dictionary:
            del dictionary[key]

def mp_range_read(mp_file,ev_min,ev_max):
    mdata = np.genfromtxt(mp_file,delimiter=',')
    mdata = np.atleast_2d(np.asarray(mdata, dtype=float))
    mdata = mdata[np.argsort(mdata[:,0])]
    mdata_slice = mdata[(mdata[:,0] <= ev_min),0]
    ev_min_nearest = max(mdata_slice) if len(mdata_slice) != 0 else mdata[0,0]
    mdata_slice = mdata[(mdata[:,0] >= ev_max),0]
    ev_max_nearest = min(mdata_slice) if len(mdata_slice) != 0 else mdata[-1,0]
    mdata = mdata[(mdata[:,0] >= ev_min_nearest) & (mdata[:,0] <= ev_max_nearest),:]
    return mdata

def mp_single_read(mp_file_or_data, ev):
    if isinstance(mp_file_or_data, str):
        mdata = np.genfromtxt(mp_file_or_data, delimiter=',')
    else:
        mdata = mp_file_or_data
    mdata = np.atleast_2d(np.asarray(mdata, dtype=float))
    mdata = mdata[np.argsort(mdata[:,0])]
    xs = mdata[:,0]
    n_re = mdata[:,1]
    n_im = mdata[:,2]
    if len(xs) == 1:
        return np.array([ev, n_re[0], n_im[0]], dtype=float)
    n_re_interp = np.interp(ev, xs, n_re)
    n_im_interp = np.interp(ev, xs, n_im)
    return np.array([ev, n_re_interp, n_im_interp], dtype=float)

def sweep_values_read(aw_parameters):
    if 'ev_range' in aw_parameters.keys():
        return [float(x) for x in aw_parameters['ev_range']]
    return [float(x) for x in aw_parameters['nm_range']]

def sweep_value_to_dirname(value):
    return pretty_value(value)

def any_value_to_dirname(value):
    if isinstance(value, (int, float, np.integer, np.floating)):
        return sweep_value_to_dirname(value)
    return str(value)

def pretty_value(value):
    if isinstance(value, (int, float, np.integer, np.floating)):
        return format(float(value), '.12g')
    return str(value)

def cmdline_construct(aw_parameters,adda_cmdlineargs):
    cmdline = aw_parameters["adda_exec"]
    for arg in adda_cmdlineargs:
        cmdline += f" -{arg} {adda_cmdlineargs[arg]}"
    return cmdline

def ev_to_nm(ev):
    return 1239.8419842361123824 / ev

def nm_to_ev(nm):
    return 1239.8419842361123824 / nm

def parse_value(file,match):
    with open(file, "r") as file:
            for line in file:
                if (match+"\t") in line:
                    value = float(re.findall("[ \t][-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?",line)[0])
                    break
    return value

def midpoints(x):
    sl = ()
    for i in range(x.ndim):
        x = (x[sl + np.index_exp[:-1]] + x[sl + np.index_exp[1:]]) / 2.0
        sl += np.index_exp[:]
    return x

def geom_gen(aw_parameters,adda_cmdlineargs,dirname="geom"):
        cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
        #cmdline = aw_parameters["adda_exec"] 
        os.makedirs(dirname, exist_ok=True)
        cmdline += f" -dir {dirname}"
        cmdline += " -prognosis"
        cmdline += " -save_geom"
        cmdline += " > " + os.devnull
        os.system(cmdline)
        #print(cmdline)

def dipole_middles(point1, point2, d, oddX, oddY):
    (x_left,y_bottom) = point1
    (x_right,y_top) = point2
    left = math.floor(x_left/d) - oddX
    right = math.ceil(x_right/d) + oddX
    bottom = math.floor(y_bottom/d) - oddY
    top = math.ceil(y_top/d) + oddY
    # print(left,right)
    # print(bottom,top)
    # print(oddX,oddY)
    # print()
    x0s = np.linspace(left*d, right*d, abs(round(right - left))+1)
    x0s = np.around(x0s,2)
    y0s = np.linspace(bottom*d, top*d, abs(round(top - bottom))+1)
    y0s = np.around(y0s,2)
    return x0s, y0s

def scan_grid(aw_parameters,adda_cmdlineargs):
    # if "coords" in aw_parameters:
    #     print("\n! Warning ! Coords already calculated & prop already rotated\n")
    x_deg = aw_parameters["x_deg"]*(math.pi/180)
    z_deg = aw_parameters["z_deg"]*(math.pi/180)
    shift = aw_parameters["shift"]
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    x_left,x_right = aw_parameters["scan_x_range"]
    y_bottom,y_top = aw_parameters["scan_y_range"]
    step = aw_parameters["scan_step"]
    
    tempdir = "_temp_rotating"
    geom_gen(aw_parameters,adda_cmdlineargs,tempdir)
    #print(os.listdir('temp_rotating'))
    for file in os.listdir(tempdir):
        if file.endswith(".geom"):
            data = geometry_read(f"_temp_rotating/{file}")
            #print(data)
            gx, gy, gz = len(np.unique(data[:,0])), len(np.unique(data[:,1])), len(np.unique(data[:,2]))
            #print(gx,gy,gz)
    shutil.rmtree(tempdir, ignore_errors=True)
    
    #adjusting area so the points are exactly in the middle between the dipoles
    d = size/grid #nm
    oddX = 0.5*(gx % 2)
    oddY = 0.5*(gy % 2)
    #print(oddX,oddY)
    x0s, y0s = dipole_middles((x_left,y_bottom),(x_right,y_top),d,oddX,oddY)
    x0s, y0s = x0s[0::step], y0s[0::step]
    # print(x0s)
    # print(y0s)
    # print()
    
    # print(aw_parameters["x_deg"])
    # print(aw_parameters["z_deg"])
    # print()
    
    #if rotated by (90 +/- 45) degrees, the odd parameter may change
    if 45 < (abs(aw_parameters["x_deg"]) % 180) < 135: 
        #print("KAWABUNGA")
        shiftY = abs((gy % 2) - (gz % 2))*d/2
        #print((gy % 2),(gz % 2),shiftY)
        y0s = y0s+shiftY
    
    if 45 < (abs(aw_parameters["z_deg"]) % 180) < 135:
        #print("KAWABUNGA")
        shiftX = abs((gx % 2) - (gy % 2))*d/2
        #print((gx % 2),(gy % 2),shiftX)
        x0s = x0s+shiftX
    
    # print(x0s)
    # print(y0s)
    
    #x0s, y0s = x0s[0::step]+d/2, y0s[0::step]+d/2
    
    px0, py0, pz0 = [float(x) for x in aw_parameters["prop"].split(" ")]
    #px0, py0, pz0 = [float(x) for x in adda_cmdlineargs["prop"].split(" ")]
    cos_th, sin_th = math.cos(x_deg), math.sin(x_deg)
    cos_phi, sin_phi = math.cos(z_deg), math.sin(z_deg)
    
    #Rotating along x-axis
    px1 = px0
    py1 = cos_th*py0 - sin_th*pz0
    pz1 = sin_th*py0 + cos_th*pz0
    
    #Rotating along z-axis
    px2 = cos_phi*px1 - sin_phi*py1
    py2 = sin_phi*px1 + cos_phi*py1
    pz2 = pz1
    
    xs, ys, zs = [], [], []
    gxs, gys = [], []
    for x0 in x0s:
        for y0 in y0s:            
            #Rotating along x-axis
            x1 = x0
            y1 = cos_th*y0 - sin_th*shift
            z1 = sin_th*y0 + cos_th*shift
            
            #Rotating along z-axis
            x2 = cos_phi*x1 - sin_phi*y1# + shift*px2
            y2 = sin_phi*x1 + cos_phi*y1# + shift*py2
            z2 = z1# + shift*pz2
            
            # #Projecting on z=0 plane
            # hsub = 18
            # if abs(pz2) > 1e-10:
            #     if abs(px2) > 1e-10 : x2 = -(px2/pz2)*(z2+hsub) + x2
            #     if abs(py2) > 1e-10 : y2 = -(py2/pz2)*(z2+hsub) + y2
            #     z2 = -hsub
            
            gxs.append(x0)
            gys.append(y0)
            xs.append(x2)
            ys.append(y2)
            zs.append(z2)
    arr1 = [0,-sin_th*(shift),cos_th*(shift)]
    aw_parameters["arrow"] = [cos_phi*arr1[0] - sin_phi*arr1[1],sin_phi*arr1[0] + cos_phi*arr1[1], arr1[2]]
    adda_cmdlineargs["prop"] = f"{round(px2,3)} {round(py2,3)} {round(pz2,3)}"
    aw_parameters["coords"] = np.array([xs,ys,zs,gxs,gys]).T

def geometry_read(geom_path):
    lines = []
    with open(geom_path, 'r') as fp:
        # read an store all lines into list
        lines = fp.readlines()
    
    with open("temp.geom", 'w') as fp:
        for number, line in enumerate(lines):
            if line[0] != "N":
                fp.write(line)
    data = np.genfromtxt("temp.geom", delimiter=' ', dtype="int")
    os.remove("temp.geom")
    #print(data)
    return data

def geometry(geom_path,colorlist=None,dpl=1,half=False,alpha=0.9, edgecolors="black", exterior_only=True): 
    data = geometry_read(geom_path)
    if half == True:
        data = data[data[:,2]<.5*(max(data[:,2]) + min(data[:,2]))] #half
    xs, ys, zs = data[:,0], data[:,1], data[:,2]
    xs, ys, zs = xs-min(xs), ys-min(ys), zs-min(zs)
    x, y, z = dpl*np.indices((max(xs)-min(xs)+2,max(ys)-min(ys)+2,max(zs)-min(zs)+2)).astype(float)
    x = x - max(xs*dpl)/2 - max(xs%2)*0.5*dpl
    y = y - max(ys*dpl)/2 - max(ys%2)*0.5*dpl
    z = z - max(zs*dpl)/2 - max(zs%2)*0.5*dpl
    
    xc = midpoints(x)
    voxels = np.zeros(xc.shape, dtype=bool)
    colors = np.empty(xc.shape, dtype=object)
    if data.shape[1] == 4:
        ms = data[:,3] - 1
    else:
        if colorlist == None:
            ms = np.ones(xs.shape, dtype="int")
        else:
            ms = np.zeros(xs.shape, dtype="int")
    if colorlist == None:
        colorlist = ["deepskyblue", "silver", "gold", "yellowgreen", "tomato", "darkviolet", "peru", "lime"]
    for i in range(len(xs)):
        voxels[xs[i],ys[i],zs[i]] = True
        colors[xs[i],ys[i],zs[i]] = colorlist[ms[i]]
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    if exterior_only:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        nx, ny, nz = voxels.shape
        occupied = set(zip(xs, ys, zs))

        x0 = dpl*np.arange(nx + 1) - max(xs*dpl)/2 - max(xs%2)*0.5*dpl
        y0 = dpl*np.arange(ny + 1) - max(ys*dpl)/2 - max(ys%2)*0.5*dpl
        z0 = dpl*np.arange(nz + 1) - max(zs*dpl)/2 - max(zs%2)*0.5*dpl

        face_defs = [
            ((-1, 0, 0), lambda i, j, k: [(x0[i], y0[j], z0[k]), (x0[i], y0[j+1], z0[k]), (x0[i], y0[j+1], z0[k+1]), (x0[i], y0[j], z0[k+1])]),
            ((1, 0, 0), lambda i, j, k: [(x0[i+1], y0[j], z0[k]), (x0[i+1], y0[j], z0[k+1]), (x0[i+1], y0[j+1], z0[k+1]), (x0[i+1], y0[j+1], z0[k])]),
            ((0, -1, 0), lambda i, j, k: [(x0[i], y0[j], z0[k]), (x0[i], y0[j], z0[k+1]), (x0[i+1], y0[j], z0[k+1]), (x0[i+1], y0[j], z0[k])]),
            ((0, 1, 0), lambda i, j, k: [(x0[i], y0[j+1], z0[k]), (x0[i+1], y0[j+1], z0[k]), (x0[i+1], y0[j+1], z0[k+1]), (x0[i], y0[j+1], z0[k+1])]),
            ((0, 0, -1), lambda i, j, k: [(x0[i], y0[j], z0[k]), (x0[i+1], y0[j], z0[k]), (x0[i+1], y0[j+1], z0[k]), (x0[i], y0[j+1], z0[k])]),
            ((0, 0, 1), lambda i, j, k: [(x0[i], y0[j], z0[k+1]), (x0[i+1], y0[j], z0[k+1]), (x0[i+1], y0[j+1], z0[k+1]), (x0[i], y0[j+1], z0[k+1])]),
        ]

        faces = []
        facecolors = []
        for i, j, k in occupied:
            for (di, dj, dk), make_face in face_defs:
                if (i + di, j + dj, k + dk) not in occupied:
                    faces.append(make_face(i, j, k))
                    facecolors.append(colors[i, j, k])

        collection = Poly3DCollection(
            faces, facecolors=facecolors, edgecolors=edgecolors,
            linewidth=0.1, alpha=alpha, zorder=1
        )
        ax.add_collection3d(collection)
        ax.set_xlim(x0[0], x0[-1])
        ax.set_ylim(y0[0], y0[-1])
        ax.set_zlim(z0[0], z0[-1])
    else:
        ax.voxels(x, y, z, voxels, facecolors=colors, edgecolors=edgecolors, linewidth=0.1, alpha=alpha, zorder=1)
    
    # ax.set_xlim([2*min(x), 2*max(x)])
    # ax.set_ylim([0, max(y)])
    # ax.set_zlim([0, max(z)])
    ax.set(xlabel='x', ylabel='y', zlabel='z')
    plot_setrcparams()
    #ax.set_box_aspect([max(xs)-min(xs)+1,max(ys)-min(ys)+1,max(zs)-min(zs)+1])
    ax.set_box_aspect([ax.get_xlim()[1]-ax.get_xlim()[0]+1,ax.get_ylim()[1]-ax.get_ylim()[0]+1,ax.get_zlim()[1]-ax.get_zlim()[0]+1])
    ax.view_init(30, -120)
    return fig,ax


def scan_geometry(aw_parameters, geom_path, colorlist=None, dpl=1, half=False, alpha=0.9, edgecolors="black",
                  scan_color="dodgerblue", scan_alpha=0.35, scan_size=18, show_arrow=False,
                  arrow_color="royalblue", arrow_scale=None, depthshade=False, exterior_only=True):
    """Visualize a geometry together with the scan grid generated by scan_grid().

    Parameters
    ----------
    aw_parameters : dict
        Must contain aw_parameters["coords"], usually produced by scan_grid().
    geom_path : str
        Path to .geom file that is displayed through geometry().
    dpl : float
        Physical dipole size in nm per dipole index. Use size/grid for the current geometry.
    show_arrow : bool
        If True and aw_parameters contains "arrow", draws the scan-plane normal / beam axis.

    Returns
    -------
    fig, ax : matplotlib figure and 3D axis
    """
    if "coords" not in aw_parameters:
        raise KeyError('aw_parameters["coords"] is missing. Call scan_grid() first.')

    fig, ax = geometry(
        geom_path, colorlist=colorlist, dpl=dpl, half=half,
        alpha=alpha, edgecolors=edgecolors, exterior_only=exterior_only,
    )

    coords = np.asarray(aw_parameters["coords"], dtype=float)
    if coords.ndim != 2 or coords.shape[1] < 3:
        raise ValueError('aw_parameters["coords"] must be an (N, >=3) array')

    xyz = coords[:, 0:3]
    ax.scatter(xyz[:,0], xyz[:,1], xyz[:,2], s=scan_size, c=scan_color, alpha=scan_alpha,
               depthshade=depthshade, linewidths=0, zorder=5)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    zmin, zmax = ax.get_zlim()

    xmin = min(xmin, float(np.min(xyz[:,0])))
    xmax = max(xmax, float(np.max(xyz[:,0])))
    ymin = min(ymin, float(np.min(xyz[:,1])))
    ymax = max(ymax, float(np.max(xyz[:,1])))
    zmin = min(zmin, float(np.min(xyz[:,2])))
    zmax = max(zmax, float(np.max(xyz[:,2])))

    if show_arrow and "arrow" in aw_parameters:
        arrow = np.asarray(aw_parameters["arrow"], dtype=float)
        nrm = np.linalg.norm(arrow)
        if nrm > 0:
            origin = np.mean(xyz, axis=0)
            span = max(xmax-xmin, ymax-ymin, zmax-zmin)
            length = 0.2*span if arrow_scale is None else arrow_scale
            direction = arrow / nrm * length
            ax.quiver(origin[0], origin[1], origin[2], direction[0], direction[1], direction[2],
                      color=arrow_color, linewidth=2.0, arrow_length_ratio=0.15, zorder=6)
            xmin = min(xmin, origin[0], origin[0] + direction[0])
            xmax = max(xmax, origin[0], origin[0] + direction[0])
            ymin = min(ymin, origin[1], origin[1] + direction[1])
            ymax = max(ymax, origin[1], origin[1] + direction[1])
            zmin = min(zmin, origin[2], origin[2] + direction[2])
            zmax = max(zmax, origin[2], origin[2] + direction[2])

    span = max(xmax-xmin, ymax-ymin, zmax-zmin)
    if span <= 0:
        span = 1.0
    xmid, ymid, zmid = (xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2
    halfspan = 0.55*span
    ax.set_xlim([xmid-halfspan, xmid+halfspan])
    ax.set_ylim([ymid-halfspan, ymid+halfspan])
    ax.set_zlim([zmid-halfspan, zmid+halfspan])
    ax.set_box_aspect([1,1,1])

    return fig, ax

def field_projection(aw_parameters, ax1, ax2, energy, file):
    axes = dict(x = 0, y = 1, z = 2)
    data = np.genfromtxt(f"{aw_parameters['dirname']}/ADDA_output/{energy}/{file}", delimiter=' ', skip_header=1)
    data = np.array(data[:,0:4])
        
    out = []
    for i in np.unique(data[:,axes[ax1]]):
        for j in np.unique(data[:,axes[ax2]]):
            out.append([i,j,0])
    out = np.asarray(out)
    
    for ind,s in enumerate(out):
        data1 = data[
            (data[:,axes[ax1]] == s[0]) & (data[:,axes[ax2]] == s[1])
            ]
        count = 0
        for i in data1:
            count = count + 1
            out[ind][2] += i[3]
        if count > 0: out[ind][2] = out[ind][2]/count
        #print(count)
    
    
    #out[:,2] = out[:,2]/max(out[:,2])
    
    with open(f"{aw_parameters['dirname']}/field_projection_{energy}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow([ax1,ax2,"|E|^2"])
        writer.writerows(out)

def plot_create():
    fig = plt.figure(constrained_layout=True)
    plot_setrcparams()
    ax = fig.add_subplot(1, 1, 1)
    #ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0,0))
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel("Energy, eV")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which="both", axis="x", linestyle="dotted", zorder=0)
    ax.tick_params(bottom=True, top=True, left=True, right=True, which = "both")
    ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    return fig, ax

def plot_setaspect(ax,ratio=.5625):
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ytop-ybottom))*ratio)

def plot_setrcparams():
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 16
    
    # SMALL_SIZE = 16
    # MEDIUM_SIZE = 18
    # BIGGER_SIZE = 16
    
    # SMALL_SIZE = 24
    # MEDIUM_SIZE = 32
    # BIGGER_SIZE = 32
    
    plt.rc('font', **{'family': 'serif', 'serif': 'Arial'})
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    
    #plt.rc('legend', fontsize=32)    # legend fontsize

def exec_cmdlines(cmdlines,parallel_procs):
    pool = multiprocessing.Pool(parallel_procs)
    result_list_tqdm = []
    for result in tqdm.tqdm(pool.imap(os.system, cmdlines, 1), total=len(cmdlines)):
        result_list_tqdm.append(result)
    return result_list_tqdm

def varyany_execute(aw_parameters,adda_cmdlineargs,var,var_range):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    dirname = aw_parameters["dirname"]
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Vary any parameter: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    if var in adda_cmdlineargs:
        del adda_cmdlineargs[var]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"Varying {var} from {pretty_value(var_range[0])} to {pretty_value(var_range[-1])}",dirname)
    cmdlines = []
    for i in var_range:
        cmdline_i = cmdline
        output_dir = os.path.abspath(dirname + f"/ADDA_output/{any_value_to_dirname(i)}")
        cmdline_i += f' -dir "{output_dir}"'
        cmdline_i += f" -{var} %s" % str(i)
        cmdline_i += " > " + os.devnull
        cmdlines.append(cmdline_i)
    #print(cmdlines)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)
    
def varyany_collect(match,aw_parameters,silent=False):
    dirname = aw_parameters["dirname"]
    xs = [d.name for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()]
    try:
        xs = sorted(xs, key=float)
    except ValueError:
        xs = sorted(xs)
    #print(xs)
    values = []
    for x in xs:
        values.append(parse_value(f"{dirname}/ADDA_output/{x}/CrossSec-Y",match))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["var",match])
        writer.writerows(zip(xs,values))
        print_log(f"Saved {dirname}/{match}.csv", silent=silent)
        
def varyany_plot(match,aw_parameters):
    dirname = aw_parameters["dirname"]
    data = np.genfromtxt(f"{dirname}/{match}.csv", delimiter=',')[1:]
    fig,ax = plot_create()
    ax.plot(data[:,0], data[:,1], label=label_for_plot(match), color=color_for_plot(match), marker=".", markersize=10)
    ax.set_xlim([min(data[:,0]),max(data[:,0])])
    ax.set_ylabel(label_for_plot(match))
    #ax.legend()
    plot_setaspect(ax)
    ax.set_xlabel("")
    fig.savefig(f"{dirname}/{match}.png", bbox_inches='tight')
    print_log(f"Saved {dirname}/{match}.pdf")

def spectrum_execute(aw_parameters,adda_cmdlineargs):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    dirname = aw_parameters["dirname"]
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Spectrum: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    delkeys_silent(adda_cmdlineargs, ["lambda"])
    if 'ms_files' in aw_parameters.keys():
        surfheight = adda_cmdlineargs["surf"].split(" ")[0]
        delkeys_silent(adda_cmdlineargs, ["surf"])
        msdata = np.genfromtxt(aw_parameters["ms_files"][0],delimiter=',')
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    mp_files = aw_parameters["mp_files"]
    if 'ev_range' in aw_parameters.keys():
        evs = sweep_values_read(aw_parameters)
    else:
        evs = sweep_values_read(aw_parameters)
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    mps = dict()
    for ev in evs:
        mps[ev] = ""
        for i in mdata:
            nk = mp_single_read(mdata[i], ev)
            nk = str(nk[1]) + " " + str(nk[2]) + " "
            mps[ev] += str(nk)
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    if 'ev_range' in aw_parameters.keys():
        print_log(f"Varying energy from {pretty_value(evs[0])} to {pretty_value(evs[-1])} eV",dirname)
    else:
        print_log(f"Varying wavelength from {pretty_value(evs[0])} to {pretty_value(evs[-1])} nm",dirname)
    cmdlines = []
    for ev in evs:
        cmdline_i = cmdline
        output_dir = os.path.abspath(dirname + f"/ADDA_output/{sweep_value_to_dirname(ev)}")
        cmdline_i += f' -dir "{output_dir}"'
        if 'ev_range' in aw_parameters.keys():
            cmdline_i += " -lambda %s" % ev_to_nm(ev)
        else:
            cmdline_i += " -lambda %s" % ev
        cmdline_i += " -m " + f" {mps[ev]} "
        if 'ms_files' in aw_parameters.keys():
            nk = mp_single_read(msdata, ev)
            nk = str(nk[1]) + " " + str(nk[2])
            cmdline_i += " -surf " + f" {surfheight} {nk}"
        cmdline_i += " > " + os.devnull
        cmdlines.append(cmdline_i)
    #print(cmdlines)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)
    
def spectrum_collect(match,aw_parameters,pol="Y",silent=False):
    dirname = aw_parameters["dirname"]
    ev_dirs = sorted([d.name for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()], key=float)
    evs = [float(ev_dir) for ev_dir in ev_dirs]
    values = []
    for ev_dir in ev_dirs:
        values.append(parse_value(f"{dirname}/ADDA_output/{ev_dir}/CrossSec-{pol}",match))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        if 'ev_range' in aw_parameters.keys():
            writer.writerow(["eV",match])
        else:
            writer.writerow(["nm",match])
        writer.writerows(zip(evs,values))
        print_log(f"Saved {dirname}/{match}.csv", silent=silent)

def spectrum_plot(match,aw_parameters,format="svg"):
    dirname = aw_parameters["dirname"]
    data = np.genfromtxt(f"{dirname}/{match}.csv", delimiter=',')[1:]
    fig,ax = plot_create()
    if 'ev_range' in aw_parameters.keys():
        ax.set_xlabel("Energy, eV")
    else:
        ax.set_xlabel("Wavelength, nm")
    #print(data)
    ax.plot(data[:,0], data[:,1], label=label_for_plot(match), color=color_for_plot(match), linewidth=3)
    ax.set_xlim([min(data[:,0]),max(data[:,0])])
    ax.set_ylabel(label_for_plot(match))
    #ax.legend()
    #plot_setaspect(ax)
    #ax.set_yscale('log')
    plt.savefig(f"{dirname}/{match}.{format}", bbox_inches='tight',dpi=600)
    print_log(f"Saved {dirname}/{match}.{format}")
    return fig,ax
    
def spectrumline_execute(aw_parameters,adda_cmdlineargs):
    dirname = aw_parameters["dirname"]
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    print()
    print_log("--- Spectrum for a set of points on a line: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    mp_files = aw_parameters["mp_files"]
    if 'ev_range' in aw_parameters.keys():
        evs = sweep_values_read(aw_parameters)
    else:
        evs = sweep_values_read(aw_parameters)
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    
    mps = dict()
    for ev in evs:
        mps[ev] = ""
        for i in mdata:
            nk = mp_single_read(mdata[i], ev)
            nk = str(nk[1]) + " " + str(nk[2]) + " "
            mps[ev] += str(nk)
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    point1 = aw_parameters["spectrumline_startpoint"]
    point2 = aw_parameters["spectrumline_endpoint"]
    howmanypoints = aw_parameters["spectrumline_points"]
    x_step, y_step = (point2[0]-point1[0])/(howmanypoints-1), (point2[1]-point1[1])/(howmanypoints-1)
    #adjusting area so the points are exactly in the middle between the dipoles
    d = size/grid #nm
    odd = 0.5*(grid % 2)
    x0s, y0s = dipole_middles(point1,point2,d,odd,odd)
    # print(x0s)
    # print(y0s)
    points = []
    for i in range(0,howmanypoints):
        x0_i = point1[0] + x_step*i
        y0_i = point1[1] + y_step*i
        x_i = min(x0s, key=lambda val: abs(val - x0_i))
        y_i = min(y0s, key=lambda val: abs(val - y0_i))
        #print(i, x0_i, x_i, y0_i, y_i)
        points.append((x_i,y_i))
    points = np.array(points)
    #print(points)
    
    # fig = plt.figure(constrained_layout=True)
    # ax = fig.add_subplot(1, 1, 1)
    # ax.plot(points[:,0], points[:,1], linewidth=3, marker="o")
    # #return
    
    beam_list = adda_cmdlineargs["beam_center"].split(" ")
    delkeys_silent(adda_cmdlineargs, ["lambda","m","beam_center"])
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    print_log(f"dipole size = {d} nm",dirname)
    print_log(f"Varying position from {tuple(points[0,:])} nm to {tuple(points[-1,:])} nm ({howmanypoints} points)",dirname)
    if 'ev_range' in aw_parameters.keys():
        print_log(f"Varying energy from {pretty_value(evs[0])} to {pretty_value(evs[-1])} eV",dirname)
    else:
        print_log(f"Varying wavelength from {pretty_value(evs[0])} to {pretty_value(evs[-1])} nm",dirname)
    
    cmdlines = []
    counter = 0
    for p in points:
        counter += 1
        point_dir = f"{dirname}/" + "{:03d}".format(counter) + f"_{p[0]}_{p[1]}"
        #print(point_dir)
        os.mkdir(point_dir)
        os.mkdir(point_dir+"/ADDA_output")
        beam_list[0], beam_list[1] = str(p[0]), str(p[1])
        beam_center = (" ").join(beam_list)
        #print(beam)
        for ev in evs:
            cmdline_i = cmdline
            cmdline_i += f" -beam_center {beam_center}"
            dir_fixed = os.path.abspath(point_dir + f"/ADDA_output/{sweep_value_to_dirname(ev)}")
            cmdline_i += f' -dir "{dir_fixed}"'
            if 'ev_range' in aw_parameters.keys():
                cmdline_i += " -lambda %s" % ev_to_nm(ev)
            else:
                cmdline_i += " -lambda %s" % ev
            cmdline_i += " -m" + f" {mps[ev]} "
            cmdline_i += " > " + os.devnull
            cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)
    
def spectrumline_collect(match, aw_parameters):
    dirname = aw_parameters["dirname"]
    dirs = sorted([d.name for d in os.scandir(dirname) if d.is_dir()])
    #print(dirs)
    points = []
    ys = []
    for dir_i in dirs:
        aw_parameters["dirname"] = f"{dirname}/{dir_i}"
        spectrum_collect(match,aw_parameters, silent=True)
        ys_i = np.genfromtxt(f"{dirname}/{dir_i}/{match}.csv", delimiter=',')[1:,1]
        ys.append(ys_i)
        points.append(dir_i.split("_"))
    aw_parameters["dirname"] = dirname
    #points = np.array(points)
    #print(points)
    xs = np.genfromtxt(f"{dirname}/{dirs[0]}/{match}.csv", delimiter=',')[1:,0]
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerows(zip(["Point no.", "x [nm]", "y [nm]"],*points))
        writer.writerow("-"*(len(points)+1))
        if 'ev_range' in aw_parameters.keys():
            valuenames = ["eV"] + [f'{match}']*len(points)
        else:
            valuenames = ["nm"] + [f'{match}']*len(points)
        
        writer.writerow(valuenames)
        writer.writerows(zip(xs,*ys))
        print(f"Saved to {dirname}/{match}.csv")

def spectrumline_plot(match, aw_parameters, format="svg"):
    dirname = aw_parameters["dirname"]
    #alldata = np.genfromtxt(f"{dirname}/{match}.csv", delimiter=',',dtype=None, encoding=None)
    data = np.genfromtxt(f"{dirname}/{match}.csv", delimiter=',')
    fig1,ax1 = plot_create()
    ax1.set_ylabel(label_for_plot_arbunits(match))
    #fig2,ax2 = plot_create()
    if 'ev_range' in aw_parameters.keys():
        ax1.set_xlabel("Energy, eV")
        #ax2.set_xlabel("Energy, eV")
    else:
        ax1.set_xlabel("Wavelength, nm")
        #ax2.set_xlabel("Wavelength, nm")
    #ax1.set_ylabel(r"$\it{P}$_{\rm EELS}$" + ", a. u.")
    ax1.set_ylabel(label_for_plot(match) + ", a. u.")
    #ax2.set_ylabel(label_for_plot(match) + ", a. u.")
    
    xs = data[5:,0]
    ys2 = np.zeros(np.shape(data.T[0][5:]))
    for point in data.T[1:]:
        num = int(point[0])
        ys = point[5:]/max(point[5:]) + .02*num
        ys2 += point[5:]
        #print(num)
        #print(ys)
        ax1.plot(xs, ys, color=color_for_plot(match), linewidth=1.5, zorder=(1-0.001*num))
        ax1.fill_between(xs, min(ys), ys, facecolor="white", alpha=.3, zorder=(1-0.001*num))
    ys2 /= num
    ax1.set_xlim([min(xs),max(xs)])
    #ax2.set_xlim([min(xs),max(xs)])
    #ax2.plot(xs, ys2, color=color_for_plot(match), linewidth=3)
    #ax1.set_yscale('log')
    plot_setaspect(ax1)
    #plot_setaspect(ax2)
    plt.savefig(f"{dirname}/{match}.{format}", bbox_inches='tight',dpi=600)
    print_log(f"Saved {dirname}/{match}.{format}")
    # fig2.savefig(f"{dirname}/{match}_averaged.svg", bbox_inches='tight', dpi=600)
    # print_log(f"Saved {dirname}/{match}_averaged.svg")
    return fig1, ax1, #fig2, ax2

def extrapolation_execute(aw_parameters,adda_cmdlineargs):
    dirname = aw_parameters["dirname"]
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Extrapolation: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    mp_files = aw_parameters["mp_files"]
    ev = aw_parameters["ev"]
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    mps = dict()
    mps[ev] = ""
    m_abs = 0
    for i in mdata:
        nk = mp_single_read(mdata[i], ev)
        if m_abs < math.sqrt(float(nk[1])**2 + float(nk[2])**2):
            m_abs = math.sqrt(float(nk[1])**2 + float(nk[2])**2)
        nk = str(nk[1]) + " " + str(nk[2]) + " "
        mps[ev] += str(nk) + " "
        
    lam = ev_to_nm(ev)
    adda_cmdlineargs["lambda"] = lam
    adda_cmdlineargs["m"] = f"{mps[ev]}"
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    y_min = (2*math.pi/lam)*(size/grid)*m_abs #y = k*d*|m|
    y_max = 4*y_min
    ys = np.exp(np.linspace(math.log(y_min), math.log(y_max), 9))
    grids = np.rint((2*math.pi/lam)*(size/ys)*m_abs).astype(np.int32)
    grids = np.array(list(dict.fromkeys(grids.tolist())), dtype=np.int32)
    #grids = np.asarray([128, 112, 96, 80, 64, 56, 48, 40, 32])
    del adda_cmdlineargs["grid"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    print_log(f"ev = {pretty_value(ev)}",dirname)
    print_log(f"m_p = {mps[ev]}",dirname)
    print_log(f"Varying grids: {grids}",dirname)
    cmdlines = []
    for i in grids:
        cmdline_i = cmdline
        dir_fixed = os.path.abspath(dirname + f"/ADDA_output/{i}")
        cmdline_i += f' -dir "{dir_fixed}"'
        cmdline_i += f" -grid {i}"
        cmdline_i += " > " + os.devnull
        cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)

def extrapolation_collect(match, aw_parameters):
    dirname = aw_parameters["dirname"]
    with open(f"{dirname}/ADDA_cmdlineargs.csv") as file:
        reader = csv.reader(file)
        adda_cmdlineargs = dict(reader)
    m_abs = 0
    mdata = re.split(" +",adda_cmdlineargs["m"])
    for i in range(int(len(mdata)/2)):
        if m_abs < math.sqrt(float(mdata[2*i])**2 + float(mdata[2*i+1])**2):
            m_abs = math.sqrt(float(mdata[2*i])**2 + float(mdata[2*i+1])**2)
    grids = np.array(sorted([int(d.name) for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()]))
    values = []
    for grid_i in grids:
        values.append(parse_value(f"{dirname}/ADDA_output/{grid_i}/CrossSec-Y",match))
    ys = (2*math.pi/float(adda_cmdlineargs["lambda"]))*(float(adda_cmdlineargs["size"])/grids)*m_abs #y = k*d*|m|
    weights = ys**-3
    fit,cov = np.polyfit(ys, values, 2, w=weights, cov=True)
    a = np.flip(fit)
    error = 2*np.sqrt(np.flip(np.diag(cov)))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["grids","ys","values"])
        writer.writerows(zip(grids,ys,values))
        print(f"Saved to {dirname}/{match}.csv")
    with open(f"{dirname}/{match}_fit.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["a[i]","error[i]"])
        writer.writerows(zip(a,error))
        print(f"Saved to {dirname}/{match}_fit.csv")
    
def extrapolation_plot(match, aw_parameters):
    dirname = aw_parameters["dirname"]
    data = np.genfromtxt(f"{dirname}/{match}.csv",delimiter=',')[1:]
    fig,ax = plot_create()
    ax.plot(data[:,1], data[:,2], label=label_for_plot(match)+" (simulated)", color=color_for_plot(match), marker="o", markersize=12, linestyle="none", zorder=2)
    ys_fitted = np.linspace(data[:,1][0],0,100)
    results_fit = np.genfromtxt(f"{dirname}/{match}_fit.csv",delimiter=',')[1:]
    a = results_fit[:,0]
    error = results_fit[:,1]
    points_fitted = a[0] + a[1]*ys_fitted + a[2]*ys_fitted**2
    ax.plot(ys_fitted, points_fitted, label=label_for_plot(match)+" (fit)", color="black", linewidth=3, zorder=3)
    ax.errorbar(0, a[0], yerr=error[0], color="black", linestyle="", marker="s", capsize=3, barsabove=True, label = "Confidence interval")
    ax.set_xlabel("y = kd|m|")
    ax.legend()
    plot_setaspect(ax)
    plt.savefig(f"{dirname}/{match}.svg", bbox_inches='tight')
    print_log(f"Saved {dirname}/{match}.svg")

def spectrum_with_extrapolation_execute(aw_parameters,adda_cmdlineargs):
    dirname = aw_parameters["dirname"]
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname,exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Spectrum with extrapolation: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    mp_files = aw_parameters["mp_files"]
    evs = sweep_values_read(aw_parameters)
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    mps = dict()
    m_abs = 0
    for ev in evs:
        mps[ev] = ""
        for i in mdata:
            nk = mp_single_read(mdata[i], ev)
            mre = float(nk[1])
            mim = float(nk[2])
            if m_abs < math.sqrt(mre**2 + mim**2):
                m_abs = math.sqrt(mre**2 + mim**2)
                maxmre = mre
                maxmim = mim
            nk = str(nk[1]) + " " + str(nk[2]) + " "
            mps[ev] += str(nk)
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    ev = evs[0]
    lam = ev_to_nm(ev)
    y_min = (2*math.pi/lam)*(size/grid)*m_abs #y = k*d*|m|
    y_max = 4*y_min
    ys = np.exp(np.linspace(math.log(y_min), math.log(y_max), 9))
    grids = np.rint((2*math.pi/lam)*(size/ys)*m_abs).astype(np.int32) #using same grids for all ev - they would not change
    grids = np.array(list(dict.fromkeys(grids.tolist())), dtype=np.int32)
    #grids = np.flip([220,186,156,132,110,92,78,66,56])
    adda_cmdlineargs["lambda"] = lam
    adda_cmdlineargs["m"] = f"{maxmre} {maxmim}"
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    del adda_cmdlineargs["lambda"]
    del adda_cmdlineargs["m"]
    del adda_cmdlineargs["grid"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    print_log(f"Varying energy from {pretty_value(evs[0])} to {pretty_value(evs[-1])} eV",dirname)
    print_log(f"Varying grids: {grids}",dirname)
    cmdlines = []
    for grid_i in grids:
        os.mkdir(f"{dirname}/ADDA_output/{grid_i}")
        for ev in evs:
            cmdline_i = cmdline
            dir_fixed = os.path.abspath(dirname + f"/ADDA_output/{grid_i}/{sweep_value_to_dirname(ev)}")
            cmdline_i += f' -dir "{dir_fixed}"'
            cmdline_i += f" -grid {grid_i}"
            cmdline_i += " -lambda %s" % ev_to_nm(ev)
            cmdline_i += " -m " + f" {mps[ev]} "
            cmdline_i += " > " + os.devnull
            cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)

def spectrum_with_extrapolation_collect(match, aw_parameters):
    dirname = aw_parameters["dirname"]
    with open(f"{dirname}/ADDA_cmdlineargs.csv") as file:
        reader = csv.reader(file)
        adda_cmdlineargs = dict(reader)
    mre = float(adda_cmdlineargs["m"].split(" ")[0])
    mim = float(adda_cmdlineargs["m"].split(" ")[1])
    m_abs = math.sqrt(mre**2 + mim**2)
    grids = np.array(sorted([int(d.name) for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()]))
    ev_dirs = sorted([d.name for d in os.scandir(f"{dirname}/ADDA_output/{grids[0]}") if d.is_dir()], key=float)
    evs = [float(ev_dir) for ev_dir in ev_dirs]
    ys = (2*math.pi/float(adda_cmdlineargs["lambda"]))*(float(adda_cmdlineargs["size"])/grids)*m_abs #y = k*d*|m|
    weights = ys**-3
    fit_values = []
    fit_errors = []
    for ev_dir in ev_dirs:
        values = []
        for grid_j in grids:
            values.append(parse_value(f"{dirname}/ADDA_output/{grid_j}/{ev_dir}/CrossSec-Y",match))
        fit,cov = np.polyfit(ys, values, 2, w=weights, cov=True)
        a = np.flip(fit)
        error = 2*np.sqrt(np.flip(np.diag(cov)))
        fit_values.append(a[0])
        fit_errors.append(error[0])
    with open(f"{dirname}/{match}_fit.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["ev",match,"error"])
        writer.writerows(zip(evs,fit_values,fit_errors))
        print(f"Saved to {dirname}/{match}_fit.csv")

def spectrum_with_extrapolation_plot(match,aw_parameters,format="svg"):
    dirname = aw_parameters["dirname"]
    data = np.genfromtxt(f"{dirname}/{match}_fit.csv",delimiter=',')[1:]
    fig,ax = plot_create()
    ax.plot(data[:,0], data[:,1], label=label_for_plot(match), color=color_for_plot(match), linewidth=3)
    ax.fill_between(data[:,0], data[:,1]-data[:,2], data[:,1]+data[:,2], label="Confidence interval", color="blue", alpha=0.2)
    ax.set_xlim([min(data[:,0]),max(data[:,0])])
    ax.legend()
    plot_setaspect(ax)
    plt.savefig(f"{dirname}/{match}.{format}", bbox_inches='tight',dpi=600)
    print_log(f"Saved {dirname}/{match}.{format}")

def scan_execute(aw_parameters,adda_cmdlineargs):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    dirname = aw_parameters["dirname"]
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Scan: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    mp_files = aw_parameters["mp_files"]
    if "ev" in aw_parameters.keys():
        ev = aw_parameters["ev"]
        adda_cmdlineargs["lambda"] = ev_to_nm(ev)
    else:
        ev = aw_parameters["nm"]
        adda_cmdlineargs["lambda"] = ev
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    mps = dict()
    mps[ev] = ""
    for i in mdata:
        nk = mp_single_read(mdata[i], ev)
        nk = str(nk[1]) + " " + str(nk[2]) + " "
        mps[ev] += str(nk) + " "
    adda_cmdlineargs["m"] = f"{mps[ev]}"
    if 'ms_files' in aw_parameters.keys():
        surfheight = adda_cmdlineargs["surf"].split(" ")[0]
        delkeys_silent(adda_cmdlineargs, ["surf"])
        msdata = np.genfromtxt(aw_parameters["ms_files"][0],delimiter=',')
        nk = mp_single_read(msdata, ev)
        nk = str(nk[1]) + " " + str(nk[2])
        adda_cmdlineargs["surf"] = f"{surfheight} {nk}"
    
    # size = adda_cmdlineargs["size"]
    # grid = adda_cmdlineargs["grid"]
    x_left,x_right = aw_parameters["scan_x_range"]
    y_bottom,y_top = aw_parameters["scan_y_range"]
    # step = aw_parameters["scan_step"]
    # #adjusting area so the points are exactly in the middle between the dipoles
    # d = size/grid #nm
    # print_log(f"dipole size = {d} nm",dirname)
    # odd = 0.5*(grid % 2)
    # x0s, y0s = dipole_middles((x_left,y_bottom),(x_right,y_top),d,odd)
    # x0s, y0s = x0s[0::step], y0s[0::step]
    
    del adda_cmdlineargs["beam_center"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    if "ev" in aw_parameters.keys():
        print_log(f"ev = {pretty_value(ev)}",dirname)
    else:
        print_log(f"nm = {pretty_value(ev)}",dirname)
    
    print_log(f"m_p = {mps[ev]}",dirname)
    print_log(f"Varying (x_left,x_right) = ({x_left},{x_right}) nm",dirname)
    print_log(f"Varying (y_bottom,y_top) = ({y_bottom},{y_top}) nm",dirname)
    cmdlines = []
    
    for i in aw_parameters["coords"]:
        cmdline_i = cmdline
        dir_fixed = os.path.abspath(dirname + f"/ADDA_output/{round(i[3],2)}_{round(i[4],2)}")
        cmdline_i += f' -dir "{dir_fixed}"'
        beam_center = f"{round(i[0],2)} {round(i[1],2)} {round(i[2],2)}"
        cmdline_i += f" -beam_center {beam_center}"
        cmdline_i += " > " + os.devnull
        cmdlines.append(cmdline_i)
    
    # for x0_i in x0s:
    #     for y0_i in y0s:
    #         cmdline_i = cmdline
    #         dir_fixed = os.path.abspath(dirname + f"/ADDA_output/{x0_i}_{y0_i}")
    #         cmdline_i += f' -dir "{dir_fixed}"'
    #         beam_list[0], beam_list[1] = str(x0_i), str(y0_i)
    #         beam_center = (" ").join(beam_list)
    #         cmdline_i += f" -beam_center {beam_center}"
    #         cmdline_i += " > " + os.devnull
    #         cmdlines.append(cmdline_i)
    
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)

def scan_collect(match, aw_parameters):
    aw_parameters = dict(aw_parameters)
    dirname = aw_parameters["dirname"]
    dirs = sorted([d.name for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()])
    xs = []
    ys = []
    values = []
    for dir_i in dirs:
        xy = dir_i.split("_")
        xs.append(float(xy[0]))
        ys.append(float(xy[1]))
        values.append(parse_value(f"{dirname}/ADDA_output/{dir_i}/CrossSec-Y",match))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["x","y",match])
        writer.writerows(zip(xs,ys,values))
        print(f"Saved to {dirname}/{match}.csv")

def scan_plot(match, aw_parameters, details=True, norm=None, format="svg", cmap="turbo"):
    aw_parameters = dict(aw_parameters)
    dirname = aw_parameters["dirname"]
    with open(f"{dirname}/ADDA_cmdlineargs.csv") as file:
        reader = csv.reader(file)
        adda_cmdlineargs = dict(reader)
    size = float(adda_cmdlineargs["size"])
    grid = float(adda_cmdlineargs["grid"])
    data = np.genfromtxt(f"{dirname}/{match}.csv",delimiter=',',dtype=None, encoding=None)
    axNames = data[0].astype("str")
    data = data[1:].astype("float")
    #print(axNames)
    xs = data[:,0]
    ys = data[:,1]
    zs = data[:,2]
    Nx = len(np.unique(xs))
    Ny = len(np.unique(ys))
    ind = np.lexsort((xs,ys))
    (x,y,z) = (xs[ind].reshape((Ny, Nx)),ys[ind].reshape((Ny, Nx)),zs[ind].reshape((Ny, Nx)))
    # print(x)
    # print(y)
    # print(z)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plot_setrcparams()
    ax.set_aspect('equal')
    d = size/grid
    im = ax.imshow(z, extent=(min(xs)-d/2, max(xs)+d/2, min(ys)-d/2, max(ys)+d/2), origin="lower", cmap=cmap, norm=norm)
    #plt.scatter(x, y, c=z, marker="s") # scatter is the most stable function for visualization, so use this for debugging
    ax.set_xlabel(f"{axNames[0]}, nm")
    ax.set_ylabel(f"{axNames[1]}, nm")
    
    if details == True:
        cbar = fig.colorbar(im)
        cbar.set_label(label_for_plot(axNames[2]))
        #cbar.formatter.set_powerlimits((0, 0))
    else:
        plt.axis('off')

    plt.savefig(f"{dirname}/{match}.{format}", bbox_inches='tight',dpi=600)
    print_log(f"Saved {dirname}/{match}.{format}")
