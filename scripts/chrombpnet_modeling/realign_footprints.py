import h5py
import hdf5plugin
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse


def plot_motif(h5_file, outprefix):
    outputlen = 1000
    with h5py.File(h5_file, "r") as f:
        # Access the data stored in the file
        for motif in f.keys():
            motif_footprint = f[motif]['i0'][:]
            fig, ax = plt.subplots()
            ax.plot(range(-100,100),motif_footprint[outputlen//2-100:outputlen//2+100])
            ax.set_xticks([-100, 0 , 100])
            ax.set_ylim(0,0.005)
            ax.set_yticks([0, 0.0025, 0.005])
            fig.savefig(outprefix + ".{}.footprint.png".format(motif), bbox_inches = 'tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--h5_file',required=True,help="h5 file from motif footprinting")
    parser.add_argument('-op','--outprefix',required=True,help="Output prefix for png files")
    args = parser.parse_args()
    
    plot_motif(h5_file = args.h5_file, outprefix = args.outprefix)

    
