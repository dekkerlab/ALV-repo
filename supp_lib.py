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
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker
from more_itertools import chunked

def rstacks(fname,regs,flank=250_000,binsize=50000,fill_missing=.0):
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
        missing=fill_missing,
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
                  extra_order = None,
                  hmss_order = None,
                  fillmissing=False,
                  interpolation="nearest"
                 ):
    """
    plot a buch of stackups ...
    """
    if extra_plots is None:
        extra_plots = []
    # how many stackups are there
    num_stackups = len(hmss)
    num_rows = num_stackups + len(extra_plots) 
    # let's figure out - how tall is this stackup
    # get heights of stackups from each groups
    stackup_items = len(hmss[0])
    stackup_height = stackup_items*12/10_000
    figure_height = stackup_height + 2.5
    fig = plt.figure(
        figsize=(3.5*num_rows, figure_height),
        facecolor="white",
        constrained_layout=True
    )
    gs = fig.add_gridspec(
        1+2, # average plot, stackup and colorbar ...
        num_rows,
        width_ratios=[1]*num_rows,
        height_ratios = \
            [0.95*2.5/figure_height] + \
            [(figure_height-2.5)/figure_height] + \
            [0.05*2.5/figure_height]
    )

    ax_profile = {}
    ax_stackup = {}
    ax_xtra = OrderedDict()
    ax_cbar = {}
    # let's define order
    if extra_order is None:
        extra_order = list( range(len(extra_plots)) )
    if hmss_order is None:
        hmss_order = list( range(len(extra_plots), num_rows) )
    # replace following with the pre-defined column indexes ...
    for idx in hmss_order:
        ax_profile[idx] = fig.add_subplot(gs[0,idx])
        ax_stackup[idx] = fig.add_subplot(gs[1,idx])
        ax_cbar[idx] = fig.add_subplot(gs[2,idx])
    for idx in extra_order:
        ax_xtra[idx] = fig.add_subplot(gs[1,idx])

    hm_arr = {}
    profile_hm = {}
    for idx, hm in zip(hmss_order, hmss):
        if fillmissing:
            X = hm[:]
            missing = ~np.isfinite(X)
            mu = np.nanmean(X, axis=0, keepdims=True) # axis 0 or 1 - rows or columns ?!
            hm_arr[idx] = np.where(missing, mu, X)
        else:
            hm_arr[idx] = hm[:]
        profile_hm[idx] = np.nanmean(hm_arr[idx],axis=0)
    # turning some of the input parameters into "oredered" or labeled dicts ...
    if norms is None:
        norms = { _i:None for _i in hmss_order}
    else:
        norms = { _i:norms[i] for i,_i in enumerate(hmss_order)}
    vlims = { _i:limss[i] for i,_i in enumerate(hmss_order)}
    titles = { _i:titles[i] for i,_i in enumerate(hmss_order)}
    if binsizes is None:
        binsizes = { _i:1 for _i in hmss_order}
    else:
        binsizes = { _i:binsizes[i] for i,_i in enumerate(hmss_order)}

    for idx, cmap in zip(hmss_order, cmps):
        ax_profile[idx].plot(profile_hm[idx],"k-")
        ax_profile[idx].set_yscale("linear" if norms[idx] is None else "log")
        stack_hm = ax_stackup[idx].imshow(
                          hm_arr[idx],
                          norm=norms[idx],
                          aspect="auto",
                          vmin=vlims[idx][0],
                          vmax=vlims[idx][1],
                          cmap=cmap,
                          interpolation=interpolation
        )
        # beautify ...
        center_bin = hm_arr[idx].shape[1]/2 - .5
        ax_profile[idx].set_ylim(vlims[idx])
        ax_profile[idx].set_title(titles[idx])
        if binsizes is not None:
            flank_in_kb = int((center_bin+.5)*binsizes[idx]/1000)
            flank_ticks = [0-.5,center_bin,hm_arr[idx].shape[1]-.5]
            flank_ticklabels = [-flank_in_kb,0,flank_in_kb]
            ax_profile[idx].set_xticks(flank_ticks)
            ax_profile[idx].set_xticklabels(flank_ticklabels)
            ax_stackup[idx].set_xticks(flank_ticks)
            ax_stackup[idx].set_xticklabels(flank_ticklabels)
        plt.colorbar(stack_hm,cax=ax_cbar[idx],orientation="horizontal")
        # remove ticks
        ax_stackup[idx].set_yticks([])

    return ax_xtra

def plot_stackups_sets(
                  extra_plots,
                  hmss, # will become a dictionary now (or list) ...
                  titles,
                  limss,
                  cmps,
                  norms=None,
                  binsizes=None,
                  extra_order = None,
                  hmss_order = None,
                  fillmissing=False,
                  interpolation="nearest"
                 ):
    """
    plot a buch of stackups ...
    """
    # rewrite everyhting assuming hmss is a dict of stackup groups !
    # groups are plotted on top of each other ...
    
    if extra_plots is None:
        extra_plots = []
    # regardless - claculate number of axes for stackups ...
    num_stackup_groups = len(hmss)
    # pick in every stackup group and see how many are there
    num_stackups = max(len(hmss[k]) for k in hmss)
    num_rows = num_stackups + len(extra_plots)
    # let's figure out - how tall is this stackup
    # get heights of stackups from each groups
    stackup_group_heights = [ len(hmss[k][0]) for k in hmss ]
    stackup_height = sum(stackup_group_heights)*12/10_000
    service_height = 3.5
    width_per_col = 3.7
    figure_height = stackup_height + service_height
    fig = plt.figure(
        figsize=(width_per_col*num_rows, figure_height + 2*.8 + 2*.8),
        facecolor="white",
        tight_layout=False,
        constrained_layout=False,
    )
    gs = fig.add_gridspec(
        num_stackup_groups+2,
        num_rows,
        width_ratios=[1]*num_rows,
        height_ratios = \
            [0.95*service_height/figure_height] + \
            [(_h/sum(stackup_group_heights))*(figure_height-service_height)/figure_height for _h in stackup_group_heights] + \
            [0.05*service_height/figure_height],
        left=0,
        right=1,
        top=1 - .8/(figure_height + 2*.8 + 2*0.8),
        bottom=.8/(figure_height + 2*.8 + 2*0.8),
        hspace=((len(stackup_group_heights)+2)/figure_height)*.8,
        wspace=0.07,
    )

    ax_profile = {}
    ax_stackup = {}
    ax_xtra = OrderedDict()
    ax_cbar = {}
    # let's define order
    if extra_order is None:
        extra_order = list( range(len(extra_plots)) )
    if hmss_order is None:
        hmss_order = list( range(len(extra_plots), num_rows) )
    # replace following with the pre-defined column indexes ...
    for idx in hmss_order:
        ax_profile[idx] = fig.add_subplot(gs[0,idx])
        ax_stackup[idx] = [fig.add_subplot(gs[_i+1,idx]) for _i in range(num_stackup_groups)] # stackup groups ...
        ax_cbar[idx] = fig.add_subplot(gs[-1,idx])
    for idx in extra_order:
        ax_xtra[idx] = [fig.add_subplot(gs[_i+1,idx]) for _i in range(num_stackup_groups)] # stackup groups ...

    # turning some of the input parameters into "oredered" or labeled dicts ...
    if norms is None:
        norms = { _i:None for _i in hmss_order}
    else:
        norms = { _i:norms[i] for i,_i in enumerate(hmss_order)}
    vlims = { _i:limss[i] for i,_i in enumerate(hmss_order)}
    titles = { _i:titles[i] for i,_i in enumerate(hmss_order)}
    if binsizes is None:
        binsizes = { _i:1 for _i in hmss_order}
    else:
        binsizes = { _i:binsizes[i] for i,_i in enumerate(hmss_order)}

    hm_arr = {}
    profile_hm = {}
    # for each group of stackups (vertically set)
    for group_id, k in enumerate(hmss):
        hm_arr[group_id] = {}
        profile_hm[group_id] = {}
        # for every stackup in each group (horizontal set)
        for idx, hm in zip(hmss_order, hmss[k]):
            if fillmissing:
                X = hm[:]
                missing = ~np.isfinite(X)
                mu = np.nanmean(X, axis=0, keepdims=True) # axis 0 or 1 - rows or columns ?!
                hm_arr[group_id][idx] = np.where(missing, mu, X)
            else:
                hm_arr[group_id][idx] = hm[:]
            if norms[idx] is None:
                profile_hm[group_id][idx] = np.nanmean(hm_arr[group_id][idx],axis=0)
            else:
                profile_hm[group_id][idx] = np.exp(np.nanmean(np.log(hm_arr[group_id][idx]),axis=0))

    for idx, cmap in zip(hmss_order, cmps):
        # plot profiles from every group on a single common axis for profiles...
        for _i in range(num_stackup_groups):
            ax_profile[idx].plot(profile_hm[_i][idx],linewidth=4)
        ax_profile[idx].set_yscale("linear" if norms[idx] is None else "log")
        # stackups for every group ...
        for _i in range(num_stackup_groups):
            stack_hm = ax_stackup[idx][_i].imshow(
                              hm_arr[_i][idx],
                              norm=norms[idx],
                              aspect="auto",
                              vmin=vlims[idx][0] if norms[idx] is None else None,
                              vmax=vlims[idx][1] if norms[idx] is None else None,
                              cmap=cmap,
                              interpolation=interpolation,
            )
        # beautify ...
        group_id_beautify = 0
        first_bin = 0-.5
        center_bin = hm_arr[group_id_beautify][idx].shape[1]/2 - .5
        last_bin = hm_arr[group_id_beautify][idx].shape[1]-.5
        ax_profile[idx].set_xlim([first_bin, last_bin])
        ax_profile[idx].set_ylim(vlims[idx])
        ax_profile[idx].set_title(titles[idx])
        ax_profile[idx].tick_params(axis="y", length=0, direction="in", pad=-5)
        ax_profile[idx].tick_params(axis="x", length=6)
        ax_profile[idx].set_yticks(vlims[idx])
        ax_profile[idx].set_yticklabels(vlims[idx],fontsize=50)
        for _tidx, tick in enumerate(ax_profile[idx].yaxis.get_majorticklabels()):
            tick.set_horizontalalignment("left")
            if _tidx == 0:
                tick.set_verticalalignment("bottom")
            elif _tidx == 1:
                tick.set_verticalalignment("top")
        # human readable kb stuff:
        flank_in_kb = int((center_bin+.5)*binsizes[idx]/1000)
        flank_ticks = [first_bin, center_bin, last_bin]
        flank_ticklabels = [-flank_in_kb, "", flank_in_kb]
        ax_profile[idx].set_xticks(flank_ticks)
        ax_profile[idx].set_xticklabels(flank_ticklabels,fontsize=50)
        for _tidx, tick in enumerate(ax_profile[idx].xaxis.get_majorticklabels()):
            if _tidx == 0:
                tick.set_horizontalalignment("left")
            elif _tidx == 2:
                tick.set_horizontalalignment("right")
            else:
                tick.set_horizontalalignment("center")
        for _i in range(num_stackup_groups-1):
            ax_stackup[idx][_i].set_xticks([])
            ax_stackup[idx][_i].set_xticklabels([])
            ax_stackup[idx][_i].set_yticks([])
            ax_stackup[idx][_i].set_yticklabels([])
        # bottom one - show ticks for now ...
        _i = num_stackup_groups-1
        ax_stackup[idx][_i].set_xticks(flank_ticks)
        ax_stackup[idx][_i].set_xticklabels(flank_ticklabels,fontsize=50)
        ax_stackup[idx][_i].tick_params(axis="x", length=6)        
        ax_stackup[idx][_i].set_yticks([])
        ax_stackup[idx][_i].set_yticklabels([])
        for _tidx, tick in enumerate(ax_stackup[idx][_i].xaxis.get_majorticklabels()):
            if _tidx == 0:
                tick.set_horizontalalignment("left")
            elif _tidx == 2:
                tick.set_horizontalalignment("right")
            else:
                tick.set_horizontalalignment("center")
        plt.colorbar(stack_hm,
                    cax=ax_cbar[idx],
                    orientation="horizontal",
                    ticks=vlims[idx])
        ax_cbar[idx].tick_params(axis="x", length=6)
        ax_cbar[idx].set_xticklabels(vlims[idx],fontsize=50)
        for _tidx, tick in enumerate(ax_cbar[idx].xaxis.get_majorticklabels()):
            if _tidx == 0:
                tick.set_horizontalalignment("left")
            elif _tidx == 1:
                tick.set_horizontalalignment("right")

    return ax_xtra

def get_4cstackup(clr, features_df, binsize=None, flank=5000, fill_missing=np.nan, ignore_diags=2):
    """
    extract 4C profiles (3-bin averaged values along the rows of Hi-C matrix) for
    "small" genomic intervals in features_df ("small" ~<  binsize).
    
    This gets tricky very quickly ! (mostly because of edge-cases
    - features near chromosome starts and ends ...)
    
    To be safe, robust and reliable we do it in several steps:
     - align features with the bins of the cooler. 1 feature -> 1 bin
     - expand features according to flank, minding the chromosomal limits !
     - for every expanded feature extract a 4C profile and assign to mat separately
       for left and right "halfs" to avoid edge-case type of issues ...
       
    the center column of the returned mat corresponds to the main diagonal of Hi-C matrix
    """
    if binsize is None:
        binsize = clr.binsize
    # flank - should be multiples of binsize ...
    assert flank % binsize == 0
    flank_bins = flank // binsize
    
    # align features with the bins
    anchor_center_aligned = .5*(features_df["start"]+features_df["end"])
    anchor_center_aligned = binsize * (anchor_center_aligned/binsize).astype(int)
    #
    anchor_df = copy(features_df) # we'll be modifying it ...
    anchor_df["start"] = anchor_center_aligned
    anchor_df["end"]   = anchor_center_aligned + binsize # edge effects isn't fully solved
    
    # prepare an output matrix:
    mat = np.full((2*flank_bins+1, len(anchor_df)), fill_missing, "float")
    
    # dealing with features close to chromosome ends and starts using bioframe:
    expand_anchor_df = bioframe.expand(
        copy(anchor_df), # importan ! otherwise anchor_df is modified inplace !
        flank,
        limits=clr.chromsizes.to_dict(),
    )
    
    for n in expand_anchor_df.itertuples(index=True):
        # coordinates of an expanded anchor:
        i,c,s,e = n
        # coordinates of the true-central bin:
        _c,_s,_e = anchor_df.loc[i] # need this to deal with edge-effects
        # extract bin ids of those:
        icenter = clr.offset((_c,_s,_e))
        ifrom, ito = clr.extent((c,s,e))
        # extract a row of a Hi-C matrix (3-rows average, actually) for the expanded range:
        mat_row = np.nanmean(clr.matrix()[icenter-1:icenter+2, ifrom:ito], axis=0)
        # now carefully assign that to the "empty" mat - dealing with edge cases
        left_half = (ifrom - icenter + flank_bins, 0 + flank_bins) # relative to mat
        right_half = (0 + flank_bins, ito - icenter + flank_bins) # relative to mat
        # that's it:
        mat[slice(*left_half),i] = mat_row[:(icenter - ifrom)]
        mat[slice(*right_half),i] = mat_row[(icenter - ito):]
    # ignore diags ...
    for i in range(ignore_diags):
        mat[flank_bins + i] = np.nan
        mat[flank_bins - i] = np.nan
    #return output
    return mat.T


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