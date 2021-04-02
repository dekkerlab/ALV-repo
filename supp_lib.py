########################################
# Some notes on the installation and stuff ...
########################################
#!conda install -y -c bioconda bedtools
#!conda install -y -c bioconda ucsc-bedgraphtobigwig -> that was bad ... just a binary from UCSC - was a much better solution

#########################################
# assorted functions to help clean up notebooks
import numpy as np
import bioframe as bf
import bioframe
import pandas as pd
from sklearn.cluster import KMeans
from collections import OrderedDict
import bbi
from copy import copy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker

# now let's show a stackup of flipped and ordered EV1-s ...
cmap1 = copy(matplotlib.cm.get_cmap("RdBu_r"))
cmap1.set_bad(color='lightgrey')
cmap2 = copy(matplotlib.cm.get_cmap("RdYlBu_r"))
cmap2.set_bad(color='lightgrey')
cmap3 = copy(matplotlib.cm.get_cmap("Reds"))
cmap3.set_bad(color='lightgrey')


def rstacks(fname,regs,flank=250_000,binsize=50000):
    """
    generate a stackup from bigwig, using bbi
    """
    center_bin_ins = flank/binsize-.5
    arr = bbi.stackup(
        fname,
        regs.chrom,
        (regs.start+regs.end)/2-flank,
        (regs.start+regs.end)/2+flank,
        bins=int(2*flank/binsize),
        missing=np.nan,
        oob=np.nan,
        summary='mean'
    )
    return arr

def legacy_clst(ks, idxs, datdict, verbose=False):
    _ks = copy(ks)
    if _ks:
        k, nclust = _ks.pop(0)
        if verbose:
            print(f"clustering over {k} ...")
        dat = datdict[k][idxs,:]
        final_sorted_idx = []
        labels, centroids, _ = kmeans_missing(dat, nclust)
        clusters = np.argsort( np.mean(centroids, axis=1) )
        for i in clusters:
            clust_i_idx = np.argwhere(labels == i).flatten()
            clust_i_idx = idxs[clust_i_idx]
            sorted_index = clst( _ks, clust_i_idx, datdict )
            final_sorted_idx.append(sorted_index)
        return np.concatenate(final_sorted_idx)
    else:
        # for now let's return unaltered indexes ...
        return idxs



def recursive_data_ordering(ks, idxs, data_dict, verbose=False):
    """
    A clustering/sorting subroutine designed to organize multiple
    "stackup" arrays simultaneously (recursively).
    
    Notes
    -----
    Sorting makes the most sense within clusters, i.e. it should appear
    as the last step of the processing directives.
    
    Think of how to implement it as a class ...
    and chain those clustering/sorting operations.
    
    Parameters
    ----------
    ks : [( str, "sort", int, bool)] or [( str, "clust", int, bool)]
        List of recursive clustering/sorting directives, in
    idxs : array
        Array of indexes of the input data
    data_dict : dict
        Dictionary with the data that is used for clustering and sorting
    verbose : bool
        Enable verbose output to track clustering progress.

    Returns
    ----------
    idxs : array
        Array of indexes reordered to reflect the clustering/sorting result
    """
    _ks = copy(ks)
    if _ks:
        k, action, ngroups, ascending = _ks.pop(0)
        # extract relevant data from data_dict
        _data = data_dict[k][idxs]
        final_sorted_idx = []
        group_sizes = []
        if action == "clust":
            if verbose:
                print(f"clustering {len(_data)} elements of {k} into {ngroups} groups ...")
            # presumably - _data is 1D in this case
            assert len(_data.shape) > 1
            labels, centroids, _ = kmeans_missing(_data, ngroups)
            clusters = np.argsort( np.mean(centroids, axis=1) )
            cluster = clusters if ascending else clusters[::-1]
            # go into each subgroup/cluster:
            for i in clusters:
                clust_i_idx = np.argwhere(labels == i).flatten()
                clust_i_idx = idxs[clust_i_idx]
                _, sorted_index = recursive_data_ordering( _ks, clust_i_idx, data_dict, verbose )
                group_sizes.append( len(sorted_index) )
                final_sorted_idx.append( sorted_index )
        elif action == "sort":
            if verbose:
                print(f"sorting {len(_data)} elements of {k} into {ngroups} groups ...")
            # presumably - _data is 1D in this case
            assert len(_data.shape) == 1
            _sorted_idx = np.argsort( _data )
            _sorted_idx = _sorted_idx if ascending else _sorted_idx[::-1]
            _sorted_idx = idxs[_sorted_idx]
            # split into ngroups if requested
            for _grp_sorted_idx in np.array_split(_sorted_idx, ngroups):
                _, sorted_index = recursive_data_ordering( _ks, _grp_sorted_idx, data_dict, verbose )
                group_sizes.append( len(sorted_index) )
                final_sorted_idx.append( sorted_index )
        else:
            raise ValueError(f"{action} is not supported, use clust/sort.")
        # returning either way:
        return (group_sizes, np.concatenate(final_sorted_idx))
    else:
        # for now let's return unaltered indexes ...
        return ([1,], idxs)


# stollen from stackoverflow - source if missing
def kmeans_missing(X, n_clusters, max_iter=10):
    """Perform K-Means clustering on data with missing values.

    Args:
      X: An [n_samples, n_features] array of data to cluster.
      n_clusters: Number of clusters to form.
      max_iter: Maximum number of EM iterations to perform.

    Returns:
      labels: An [n_samples] vector of integer labels.
      centroids: An [n_clusters, n_features] array of cluster centroids.
      X_hat: Copy of X with the missing values filled in.
    """
    # Initialize missing values to their column means
    # oh WTF ? should I initialize them to their rows means ?!?!?!
    missing = ~np.isfinite(X)
    mu = np.nanmean(X, axis=0, keepdims=True) # - should I change that ?!
    X_hat = np.where(missing, mu, X)

    # cluster only when # of smaples is large enough
    if X.shape[0] > n_clusters:

        for i in range(max_iter):
            if i > 0:
                # initialize KMeans with the previous set of centroids. this is much
                # faster and makes it easier to check convergence (since labels
                # won't be permuted on every iteration), but might be more prone to
                # getting stuck in local minima.
                cls = KMeans(n_clusters, init=prev_centroids)
            else:
                # do multiple random initializations in parallel
                cls = KMeans(n_clusters)

            # perform clustering on the filled-in data
            labels = cls.fit_predict(X_hat)
            centroids = cls.cluster_centers_

            # fill in the missing values based on their cluster centroids
            X_hat[missing] = centroids[labels][missing]

            # when the labels have stopped changing then we have converged
            if i > 0 and np.all(labels == prev_labels):
                break

            prev_labels = labels
            prev_centroids = cls.cluster_centers_

        return labels, centroids, X_hat
    else:
        # don't bother with clustering if there aren't too many entries
        return np.arange(X_hat.shape[0]), X_hat, X_hat


def plot_spacing_footprint_distros(features):
    """
    little function to help us plot a couple of distributions
    a distribution of spacings between genomic features and
    a distribution of their footprints ....
    
    supply several genomic features DFs as a dictionary
    
    first one is going to be the main one ! use to define bins
    """
    # prepare figure/axis for spacings and footprint
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11,4))
    spacings = {}
    footprints = {}
    for k in features:
        _df = features[k]
        # calculate spacings, assuming everyhting is sorted
        spacings[k] = np.concatenate(
            [ np.diff(_df["start"].values) for _, _df in _df.groupby("chrom") ]
        )
        # calculate footprints of the features, after merging ...
        footprints[k] = (_df["end"] - _df["start"]).values

    # let's use the first features set as the more important one - sort of ...
    k0,*_ = features.keys()
    ref_spacings = spacings[k0]
    ref_footprints = footprints[k0]
    bins_spacing = np.geomspace(ref_spacings.min()+1,ref_spacings.max(),num=100)
    bins_footprint = np.geomspace(ref_footprints.min(),ref_footprints.max(),num=100)
    
    for k in features:
        ax1.hist(spacings[k], bins=bins_spacing, log=False, alpha=0.5, label=f"{k} {len(features[k])}")
        ax2.hist(footprints[k], bins=bins_footprint, log=False, alpha=0.5, label=f"{k} {len(features[k])}")

    ax1.set_xscale("log")
    ax1.set_xlabel("spacing between adjacent features, bp")
    ax1.set_ylabel("# of spacings")
    ax1.legend(frameon=False)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("feature footprint, bp")
    ax2.set_ylabel("# features")
    ax2.legend(frameon=False)
    
    return


def normalize_insulation_stackups_INPLACE(
    stackups_set,
    ins_keys,
    subtract = "mean_top",
    subtract_size = 7
):
    """
    take a dict of stackups for a subset of insulation stackups (ins_keys)
    subtract the maximum value (or mean of several highest values)
    from each row of a given insulation stackup in order to make
    the dip of the insulation stackup be "measured" in boundary stength units
    i.e. log fold change between dip and "shoulders" ... sorta kinda
    """
    # subtract = None
    # subtract = "shoulders"
    # subtract = "median_top"
    for k in ins_keys:
        _ddd = stackups_set[k]
        if subtract == "shoulders":
            shoulders = np.hstack( [_ddd[:,:subtract_size], _ddd[:,-subtract_size:]] )
            shoulders = np.nanmean(shoulders, axis=1, keepdims=True)
            _ddd -= shoulders
        elif subtract == "mean_top":
            _top = np.sort( _ddd )[:,-subtract_size:]
            mean_top = np.nanmean(_top, axis=1, keepdims=True)
            _ddd -= mean_top
        elif subtract == "median_top":
            _top = np.sort( _ddd )[:,-subtract_size:]
            median_top = np.nanmean(_top, axis=1, keepdims=True)
            _ddd -= median_top
        else:
            print("Operation not supported !!!")
        stackups_set[k] = _ddd
    # the end
    return

def plot_stackups(
                  extra_plots,
                  hmss,
                  titles,
                  limss,
                  cmps,
                  norms=None,
                  binsizes=None,
                  fname="xxx.png",
                  fillmissing=False,
                  interpolation="nearest"
                 ):
    """
    plot a buch of stackups ...
    """
    if extra_plots is not None:
        num_stackups = len(hmss) + len(extra_plots)
    else:
        extra_plots = []
        num_stackups = len(hmss)
    # let's figure out - how tall is this stackup
    stackup_height = len(hmss[0])*12/10_000
    figure_height = stackup_height + 2.5
    fig = plt.figure(
        figsize=(3.5*num_stackups, figure_height),
        facecolor="white",
        constrained_layout=True
    )
    gs = fig.add_gridspec(
        3,
        num_stackups,
        width_ratios=[1]*num_stackups,
        height_ratios=[
            0.95*2.5/figure_height,
            (figure_height-2.5)/figure_height,
            0.05*2.5/figure_height
        ]
    )

    ax_profile = {}
    ax_stackup = {}
    ax_cbar = {}
    for idx in range(len(extra_plots), num_stackups):
        ax_profile[idx] = fig.add_subplot(gs[0,idx])
        ax_stackup[idx] = fig.add_subplot(gs[1,idx])
        ax_cbar[idx] = fig.add_subplot(gs[2,idx])
    for idx in range(len(extra_plots)):
        ax_stackup[idx] = fig.add_subplot(gs[1,idx])
        
    for idx in range(len(extra_plots)):
        _y,_width,_color = extra_plots[idx]
        #         ax_stackup[idx].barh(_y,_width,1,color=_color,edgecolor="dimgray")
        ax_stackup[idx].step(_width,_y,color="dimgray")
        ax_stackup[idx].fill_betweenx(_y,0,_width,color=_color,step="post")
        ax_stackup[idx].invert_yaxis()
        ax_stackup[idx].invert_xaxis()
        ax_stackup[idx].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: f"{int(x*100)}"))
        ax_stackup[idx].set_ylim(max(_y),0)
        ax_stackup[idx].set_xlim(max(_width),0)

    hm_arr = {}
    profile_hm = {}
    vlims = {}
    for idx,(hm,lms) in enumerate(zip(hmss,limss)):
        if fillmissing:
            X = hm[:]
            missing = ~np.isfinite(X)
            mu = np.nanmean(X, axis=0, keepdims=True) # axis 0 or 1 - rows or columns ?!
            hm_arr[idx] = np.where(missing, mu, X)
        else:
            hm_arr[idx] = hm[:]
        profile_hm[idx] = np.nanmean(hm_arr[idx],axis=0)
        vlims[idx] = lms

    for idx,cmap in enumerate(cmps):
        ax_profile[idx+len(extra_plots)].plot(profile_hm[idx],"k-")
        if norms is not None:
            ax_profile[idx+len(extra_plots)].set_yscale("linear" if norms[idx] is None else "log")
        stack_hm = ax_stackup[idx+len(extra_plots)].imshow(
                          hm_arr[idx],
                          norm=None if norms is None else norms[idx],
                          aspect="auto",
                          vmin=limss[idx][0],
                          vmax=limss[idx][1],
                          cmap=cmap,
                          interpolation=interpolation
        )
        # beautify ...
        center_bin = hm_arr[idx].shape[1]/2 - .5
        ax_profile[idx+len(extra_plots)].set_ylim(vlims[idx])
        ax_profile[idx+len(extra_plots)].set_title(titles[idx])
        # ax_profile[idx+len(extra_plots)].axvline(center_bin, color="grey")
        # ax_stackup[idx+len(extra_plots)].axvline(center_bin, color="grey")
        if binsizes is not None:
            flank_in_kb = int((center_bin+.5)*binsizes[idx]/1000)
            flank_ticks = [0-.5,center_bin,hm_arr[idx].shape[1]-.5]
            flank_ticklabels = [-flank_in_kb,0,flank_in_kb]
            ax_profile[idx+len(extra_plots)].set_xticks(flank_ticks)
            ax_profile[idx+len(extra_plots)].set_xticklabels(flank_ticklabels)
            ax_stackup[idx+len(extra_plots)].set_xticks(flank_ticks)
            ax_stackup[idx+len(extra_plots)].set_xticklabels(flank_ticklabels)
        plt.colorbar(stack_hm,cax=ax_cbar[idx+len(extra_plots)],orientation="horizontal")
        
    for idx in range(1, num_stackups):
        ax_stackup[idx].set_yticks([])

    plt.savefig(fname)


def generate_random_bed(num_intervals, db, footprint=1_000, trunc=10_000):
    """
    generates a BED file with random
    intervals around the genome, where
    number of intervals is ~ chrom size.
    
    num_intervals - int
        number of intervals to generate
    db - str
        name of the genome assembly (bioframe pull)
    footprint - int
        interval footprint, should be < trunc
    trunc - int
        prevent intervals going over chrom length 
    """
    sortedchromsizes = bioframe.fetch_chromsizes(db, as_bed=True).reset_index(drop=True)
    # calculate proportions of intervals per chromosome (~ chrom length)
    weights = sortedchromsizes["end"].values/sortedchromsizes["end"].values.sum()
    # number of intervals per chrom, based on proportions (weights)
    number_per_chrom = (weights*num_intervals).astype(int)
    
    df = []
    for i,chrom,start,end in sortedchromsizes.itertuples():
        _number = number_per_chrom[i]
        _chrlen = end
        _intervals_starts = np.random.randint(0, _chrlen - trunc, _number)
        df.append(
            pd.DataFrame({
                "chrom": chrom,
                "start": _intervals_starts,
            })
        )
    df = pd.concat(df)
    df["end"] = df["start"] + footprint
    df = df.sort_values(["chrom","start"])
    df = df[~(df["chrom"] == "chrM")]
    return df.reset_index(drop=True)