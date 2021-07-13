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
                   "crange":(0,3.2),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["noctcf"] = {"fname":"IAA-CTCF_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,3.2),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["rad21_CTCF"] = {"fname":"NT-RAD21_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.51),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["rad21_noCTCF"] = {"fname":"IAA-RAD21_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.51),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["rad21_wt1"] = {"fname":"GSM4625024_2594.hg19.cpm.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.5),
                   "norm":None,
                   "cmap":cmap3,
                  }
# RAD21 from Rao et al 2017
samples["rad21_rad21"] = {"fname":"GSM2809609_Rao-2017-CHIP001-RAD21-untreated.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,2),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["rad21_norad21"] = {"fname":"GSM2809610_Rao-2017-CHIP002-RAD21-treated.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,2),
                   "norm":None,
                   "cmap":cmap3,
                  }

#
#
# best G4 track that we found ...
samples["G4a_r2"] = {"fname":"GSM2876095_B_REP2.SLX-12320.K562_P9_Async_a_701_517.rmdup.clean.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(2.5,28.),
                   "norm":None,
                   "cmap":cmap3,
                  }
# DRIP seq based R-loop detection in K562 ...
samples["Rloop_K562"] = {"fname":"DRIP_K562_Rloops.bw",
                   "binsize":200,
                   "flank":10_000,
                   "crange":(0,10.),
                   "norm":None,
                   "cmap":cmap3,
                  }
# add ddx for fun ...
samples["ddx_CTCF"] = {"fname":"NT-DDX55_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.52),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["ddx_noCTCF"] = {"fname":"IAA-DDX55_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.52),
                   "norm":None,
                   "cmap":cmap3,
                  }
# add TAF5L for fun ...
samples["taf5l_CTCF"] = {"fname":"NT-TAF5L_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.21),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["taf5l_noCTCF"] = {"fname":"IAA-TAF5L_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.21),
                   "norm":None,
                   "cmap":cmap3,
                  }
#############################
# ADD some replicate2-s:
samples["ddx_CTCF_r2"] = {"fname":"NT-DDX55_R1.mLb.clN.Replicate2.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.52),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["ddx_noCTCF_r2"] = {"fname":"IAA-DDX55_R1.mLb.clN.Replicate2.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.52),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["taf5l_CTCF_r2"] = {"fname":"NT-TAF5L_R1.mLb.clN.Replicate2.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.21),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["taf5l_noCTCF_r2"] = {"fname":"IAA-TAF5L_R1.mLb.clN.Replicate2.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.21),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["rad21_CTCF_r2"] = {"fname":"NT-RAD21_R2.mLb.clN.Replicate2.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["rad21_noCTCF_r2"] = {"fname":"IAA-RAD21_R1.mLb.clN.Replicate2.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["input_CTCF_r2"] = {"fname":"NT-input_R2.mLb.clN.Replicate2.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["input_noCTCF_r2"] = {"fname":"IAA-input_R1.mLb.clN.Replicate2.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["input_CTCF_r3"] = {"fname":"NT-input_R1.mLb.clN.Replicate3.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["input_noCTCF_r3"] = {"fname":"IAA-input_R2.mLb.clN.Replicate3.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "norm":None,
                   "cmap":cmap3,
                  }

#############################
# add input-control for fun ...
# input controls for RAD21 R1, DDX55 R1, TAF5L R1
samples["input_CTCF"] = {"fname":"NT-input_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["input_noCTCF"] = {"fname":"IAA-input_R1.mLb.clN.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.4),
                   "norm":None,
                   "cmap":cmap3,
                  }
# add RNA-PolII
samples["polII"] = {"fname":"PolIIChipSeq-Hap1.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.9),
                   "norm":None,
                   "cmap":cmap3,
                  }
# add H3K4 public
# HAP1_H3K4me3_WT_K4_1.bw  HAP1_H3K4me3_WT_K4_2.bw  HAP1_Input_WT_K4-control.bw
samples["h3k4_r1"] = {"fname":"HAP1_H3K4me3_WT_K4_1.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,33.),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["h3k4_r2"] = {"fname":"HAP1_H3K4me3_WT_K4_2.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,33.),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["h3k4_input"] = {"fname":"HAP1_Input_WT_K4-control.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,10.),
                   "norm":None,
                   "cmap":cmap3,
                  }

# some default h3k4me3 Chip-Seq bigWigs from ENCODE:
# https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&target.label=H3K4me3&biosample_ontology.term_name=HCT116
# https://www.encodeproject.org/files/ENCFF176NSX/@@download/ENCFF176NSX.bigWig
samples["h3k4_hct"] = {"fname":"HCT116_H3K4_2014.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,220.),
                   "norm":None,
                   "cmap":cmap3,
                  }
# ...
samples["ctcf_hct_nt1"] = {"fname":"GSM2809613_Rao-2017-CHIP005-CTCF-untreated.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,10.),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["ctcf_hct_nt2"] = {"fname":"GSM2809615_Rao-2017-CHIP007-CTCF-untreated.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,10.),
                   "norm":None,
                   "cmap":cmap3,
                  }



# add H3K27ac public
# GSM2828241_HAP1_ChIP-seq_H3K27ac_WT_1.bigWig GSM2828242_HAP1_ChIP-seq_H3K27ac_WT_2.bigWig GSM2828243_HAP1_ChIP-seq_IgG_1.bigWig
samples["h3k27ac_r1"] = {"fname":"GSM2828241_HAP1_ChIP-seq_H3K27ac_WT_1.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.1),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["h3k27ac_input"] = {"fname":"GSM2828243_HAP1_ChIP-seq_IgG_1.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,.1),
                   "norm":None,
                   "cmap":cmap3,
                  }


# add public H3K27ac for HCT - for comparisons with Rao 2017 ...
samples["h3k27ac_HCTfold"] = {"fname":"HCT-H3K27ac_foldchange.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,5),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["h3k27ac_HCTpval"] = {"fname":"HCT-H3k27ac-signalPval.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,1),
                   "norm":None,
                   "cmap":cmap3,
                  }



# add DNAse-seq public
samples["dnase-hap1"] = {"fname":"HAP1-DNAse.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,1),
                   "norm":None,
                   "cmap":cmap3,
                  }


# add DNAse-seq public
samples["dnase-HCT-r1"] = {"fname":"DNase-seq-HCT-r1.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,1.3),
                   "norm":None,
                   "cmap":cmap3,
                  }
samples["dnase-HCT-r2"] = {"fname":"DNase-seq-HCT-r2.bigWig",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,1.3),
                   "norm":None,
                   "cmap":cmap3,
                  }




# K562 data to compare with the G4 as well ...
# K562_ChIP_seqMerge_H3K4me3_v2
# K562_ChIP_exoMerge_H3K4me3.bw
samples["h3k4_k562"] = {"fname":"K562_ChIP_seqMerge_H3K4me3_v2.bw",
                   "binsize":200,
                   "flank":5_000,
                   "crange":(0,190.),
                   "norm":None,
                   "cmap":cmap3,
                  }

ins_binsize = 5_000
ins_diamond = 100_000
ins_binsize_human = f"{int(ins_binsize/1000)}kb"
ins_diamond_human = f"{int(ins_diamond/1000)}kb"
# deeper versions - combined for nice insulation stackups (todo - make sure CTCF and noCTCF - are at comparable depth)...
samples["ins_CTCF"] = {"fname":f"CkoCT442_NT_pool.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_noCTCF"] = {"fname":f"CkoCT442_IAA_pool.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
# insulation in PolII degron system ...
samples["ins_polII"] = {"fname":f"PolII-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_nopolII"] = {"fname":f"PolII-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
# insulation in RAD21 degron system ...
samples["ins_rad21"] = {"fname":f"RAD21-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_norad21"] = {"fname":f"RAD21-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }

# insulation for mutants ...
samples["ins_mutCtr_CTCF"] = {"fname":f"mutControl-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_mutCtr_noCTCF"] = {"fname":f"mutControl-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_mutDDX_CTCF"] = {"fname":f"mutDDX55-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_mutDDX_noCTCF"] = {"fname":f"mutDDX55-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_mutTAF_CTCF"] = {"fname":f"mutTAF5L-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_mutTAF_noCTCF"] = {"fname":f"mutTAF5L-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }

samples["ins_siCtr_CTCF"] = {"fname":f"siControl-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_siCtr_noCTCF"] = {"fname":f"siControl-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_siDDX_CTCF"] = {"fname":f"siDDX55-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_siDDX_noCTCF"] = {"fname":f"siDDX55-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_siTAF_CTCF"] = {"fname":f"siTAF5L-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_siTAF_noCTCF"] = {"fname":f"siTAF5L-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }



# insulation for PlaB splicing inhibitor ...
samples["ins_CtrPlaB_CTCF"] = {"fname":f"CtrlPlaB-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_CtrPlaB_noCTCF"] = {"fname":f"CtrlPlaB-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_PlaB_CTCF"] = {"fname":f"PlaB-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_PlaB_noCTCF"] = {"fname":f"PlaB-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }


# Controls ... for CTCF-degron - downsampled together to be comparable ...
samples["ins_noTIR1_500M"] = {"fname":f"Ctrl500M-noTIR1.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_WT_500M"] = {"fname":f"Ctrl500M-wtHAP1.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_CTCF_500M"] = {"fname":f"Ctrl500M-CT442-NT.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }
samples["ins_noCTCF_500M"] = {"fname":f"Ctrl500M-CT442-IAA.hg19.{ins_binsize_human}.{ins_diamond_human}.bw",
                   "binsize":ins_binsize,
                   "flank":2*ins_diamond,
                   "crange":(-.5,.5),
                   "norm":None,
                   "cmap":cmap2,
                  }




# add compartment status signal here in the samples ...
ev_binsize = 25_000
ev_binsize_human = f"{int(ev_binsize/1000)}kb"
samples["ev1_CTCF"] = {"fname": f"CkoCT442_NT_pool.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1.,1.),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_noCTCF"] = {"fname": f"CkoCT442_IAA_pool.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1.,1.),
                   "norm":None,
                   "cmap":cmap1,
                  }
# ploII degron and 
samples["ev1_polII"] = {"fname": f"PolII-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1.,1.),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_nopolII"] = {"fname": f"PolII-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1.,1.),
                   "norm":None,
                   "cmap":cmap1,
                  }
# RAD21 degron and 
samples["ev1_rad21"] = {"fname": f"RAD21-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1.,1.),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_norad21"] = {"fname": f"RAD21-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1.,1.),
                   "norm":None,
                   "cmap":cmap1,
                  }

# EV1 for mutants ...
samples["ev1_mutCtr_CTCF"] = {"fname":f"mutControl-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1.,1.),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_mutCtr_noCTCF"] = {"fname":f"mutControl-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_mutDDX_CTCF"] = {"fname":f"mutDDX55-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_mutDDX_noCTCF"] = {"fname":f"mutDDX55-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_mutTAF_CTCF"] = {"fname":f"mutTAF5L-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_mutTAF_noCTCF"] = {"fname":f"mutTAF5L-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }

samples["ev1_siCtr_CTCF"] = {"fname":f"siControl-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_siCtr_noCTCF"] = {"fname":f"siControl-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_siDDX_CTCF"] = {"fname":f"siDDX55-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_siDDX_noCTCF"] = {"fname":f"siDDX55-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_siTAF_CTCF"] = {"fname":f"siTAF5L-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_siTAF_noCTCF"] = {"fname":f"siTAF5L-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }

# EV1 for PlaB splicing inhibitor ...
samples["ev1_CtrPlaB_CTCF"] = {"fname":f"CtrlPlaB-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1.,1.),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_CtrPlaB_noCTCF"] = {"fname":f"CtrlPlaB-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_PlaB_CTCF"] = {"fname":f"PlaB-NT.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1.,1.),
                   "norm":None,
                   "cmap":cmap1,
                  }
samples["ev1_PlaB_noCTCF"] = {"fname":f"PlaB-IAA.hg19.EV.{ev_binsize_human}.cis.bw",
                   "binsize":ev_binsize,
                   "flank":6*ev_binsize,
                   "crange":(-1,1),
                   "norm":None,
                   "cmap":cmap1,
                  }


# rna-seq ...
samples["mrna_ctcf"] = {"fname": "S442_NT.bw",
                   "binsize":200,
                   "flank":5000,
                   "crange":(2.,22.),
                   "norm":matplotlib.colors.LogNorm(vmin=1.,vmax=20.),
                   "cmap":cmap3,
                  }
samples["mrna_noctcf"] = {"fname": "S442_IAA.bw",
                   "binsize":200,
                   "flank":5000,
                   "crange":(2.,22.),
                   "norm":matplotlib.colors.LogNorm(vmin=1.,vmax=20.),
                   "cmap":cmap3,
                  }