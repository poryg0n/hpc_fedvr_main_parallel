import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")   # <-- non-interactive backend
import matplotlib.pyplot as plt


def find_runs(base_dir, keyword=None):
    runs = []

    for name in os.listdir(base_dir):
        path = os.path.join(base_dir, name)

        if os.path.isdir(path):
            if keyword is None or keyword in name:
                runs.append(path)

    return runs

def load_data(filepath):
    data = np.loadtxt(filepath)
    x  = data[:,0]
    y1 = data[:,1]


    return x, y1

def read_params(filepath):
    params = {}
    with open(filepath) as f:
        for line in f:
            if "=" in line:
                key, val = line.split("=")
                params[key.strip()] = val.strip()
    return params



def get_plot_config(prefix):

    if prefix == "density":
        return {
            "xlabel": r"$x\,\rm (a.u.)$",
            "ylabel": r"$\mathcal{D}(x,t)$",
            "yscale": "log"
        }

    elif prefix == "pemd":
        return {
            "xlabel": r"$k\,(a.u.)$",
            "ylabel": r"$P(k)$",
            "yscale": "log"
        }

    elif prefix == "hhg":
        return {
            "xlabel": r"$\omega$",
            "ylabel": r"$Q^{(c)}(\omega)$",
            "yscale": "log"
        }

    else:
        return {
            "xlabel": "x",
            "ylabel": "y",
            "yscale": "linear"
        }


def format_param(val, mode="float"):
    try:
        val = float(val)
    except:
        return val

    if mode == "float":
        return f"{val:.3f}"
    elif mode == "sci":
        mant, exp = f"{val:.2e}".split("e")
        return rf"{mant} \times 10^{{{int(exp)}}}"
    else:
        return str(val)


# ---- execution part ----

output_dir = "data/exp/" 
runs = find_runs(output_dir, keyword=keyword)

if len(sys.argv) > 1:
    prefix = sys.argv[1]
else:
    prefix = "density"

if len(sys.argv) > 2:
    scale = sys.argv[2]
else:
    scale = "linear"

if len(sys.argv) > 3:
    keyword = sys.argv[3]
else:
    keyword = None



for run in runs:
#   print("Found directory:", run)

    data_file = os.path.join(run, f"{prefix}.dat")
    param_struct_file = os.path.join(run, "struct.dat")
    param_dyn_file = os.path.join(run, "dynamic.dat")

    if os.path.exists(data_file):
        x, y1 = load_data(data_file)


        plt.figure(figsize=(10, 4))

#       plt.xticks([0, 1, 2, 3])
#       plt.yticks([1e-3, 1e-2, 1e-1, 1])
        plt.xlim(min(x), max(x))
        plt.ylim(1e-14, 10)


        # ---- params textbox ----
        if os.path.exists(param_struct_file):
            struct_params = read_params(param_struct_file)
            dyn_params = read_params(param_dyn_file)

            textstr = "\n".join([
                rf"$Np = {struct_params.get('nnbr','?')}$",
                rf"$Ns = {struct_params.get('snbr','?')}$",
                rf"$xmax = {format_param(struct_params.get('xmax'))}$",
                f"\n",
                rf"$F_0 = {format_param(dyn_params.get('f0'))}$",
                rf"$\omega_0 = {format_param(dyn_params.get('omega0'))}$",
                rf"$\Delta t = {format_param(dyn_params.get('dt'), mode='sci')}$",
                rf"$n = {dyn_params.get('ntau','?')}$",
#               rf"$n_{{max}} = {params.get('nmax','?')}$",
#               rf"$dk = {params.get('dk','?')}$"
            ])

        else:
            textstr = "No params"

        plt.gca().text(
            0.05, 0.95, textstr,
            transform=plt.gca().transAxes,
            fontsize=10,
            verticalalignment='top',
#           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
#           bbox=dict(edgecolor='none', facecolor='white', alpha=0.8)
            bbox=None
        )

        plt.plot(x, y1, label=run)
#       plt.plot(x, y2, label=run)

        config = get_plot_config(prefix)
        
        plt.xlabel(config["xlabel"], fontsize=10)
        plt.ylabel(config["ylabel"], fontsize=10)
        
        if scale == "log":
            plt.yscale("log")
        else:
            plt.yscale(config["yscale"])


#       plt.legend()
        plt.grid()


        # build filename
        run_name = os.path.basename(run)
        filename = f"{prefix}_{run_name}.png"
        outpath = os.path.join(output_dir, filename)

        plt.savefig(outpath, dpi=150)
        plt.close()
#       plt.show()

        print("Loaded:", data_file)
    else:
        print("Missing file in:", run)
        
