import matplotlib
from copy import copy
# now let's show a stackup of flipped and ordered EV1-s ...
cmap1 = copy(matplotlib.cm.get_cmap("RdBu_r"))
cmap1.set_bad(color='lightgrey')
cmap2 = copy(matplotlib.cm.get_cmap("RdYlBu_r"))
cmap2.set_bad(color='lightgrey')
cmap3 = copy(matplotlib.cm.get_cmap("Reds"))
cmap3.set_bad(color='lightgrey')

cmapYlGnBu_r = copy(matplotlib.cm.get_cmap("YlGnBu_r"))
cmapYlGnBu_r.set_bad(color='lightgrey')

# it's nice sometimes to exclude sex chromosomes and mito for analysis
autosomal_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
 'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
 'chr18','chr19','chr20','chr21','chr22']

# describe samples ....
samples = {}
# each sample would have several attributres related to stackups ... - fname, binsize, flank, datarange, colormap
samples["ctcf"] = {"fname": "NT-CTCF_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,2.2),
                   "cmap":cmap3,
                  }
samples["noctcf"] = {"fname":"IAA-CTCF_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,2.2),
                   "cmap":cmap3,
                  }
samples["rad21_CTCF"] = {"fname":"NT-RAD21_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "cmap":cmap3,
                  }
samples["rad21_noCTCF"] = {"fname":"IAA-RAD21_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "cmap":cmap3,
                  }
# best G4 track that we found ...
samples["G4a_r2"] = {"fname":"GSM2876095_B_REP2.SLX-12320.K562_P9_Async_a_701_517.rmdup.clean.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,30),
                   "cmap":cmap3,
                  }
# add ddx for fun ...
samples["ddx_CTCF"] = {"fname":"NT-DDX55_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "cmap":cmap3,
                  }
samples["ddx_noCTCF"] = {"fname":"IAA-DDX55_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "cmap":cmap3,
                  }
# add RNA-PolII
samples["polII"] = {"fname":"PolIIChipSeq-Hap1.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.7),
                   "cmap":cmap3,
                  }

ins_binsize = 2_000
ins_diamond = 20_000
ins_binsize_human = f"{int(ins_binsize/1000)}kb"
ins_diamond_human = f"{int(ins_diamond/1000)}kb"
# deeper versions - combined for nice insulation stackups (todo - make sure CTCF and noCTCF - are at comparable depth)...
samples["ins_CTCF"] = {"fname":f"CkoCT442_NT_pool.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "cmap":cmap2,
                  }
samples["ins_noCTCF"] = {"fname":f"CkoCT442_IAA_pool.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "cmap":cmap2,
                  }

# add compartment status signal here in the samples ...
ev_binsize = 25_000
ev_binsize_human = f"{int(ev_binsize/1000)}kb"
samples["ev1_CTCF"] = {"fname": f"allcools.EV1.{ev_binsize_human}.bw",
                   "binsize":ev_binsize,
                   "flank":4*ev_binsize,
                   "crange":(-1.,1.),
                   "cmap":cmap1,
                  }
samples["ev1_noCTCF"] = {"fname": f"noCTCF_IAA_combined.EV1.{ev_binsize_human}.bw",
                   "binsize":ev_binsize,
                   "flank":4*ev_binsize,
                   "crange":(-1.,1.),
                   "cmap":cmap1,
                  }
