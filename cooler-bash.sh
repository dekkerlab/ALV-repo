
# request a long interactive job on the cluster
bsub -q interactive -W 04:00 -n 4 -R "rusage[mem=16000]" -R "span[hosts=1]" -Is bash

# CTCF degron no auxin - which I erroneously called wt_like_coolers_to_combine/all_combined - ~1B interactions
uris_NT=( \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC44-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC44-NT-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC44-IAA48H-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC44-IAA48H-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC442-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC442-NT-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-G1-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-NT-si-ctrl-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-NT-si-ctrl-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200805_G1sorted_442_replicate/results/coolers_library/CkoCT442-G1-NT-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
)

# noCTCF_IAA_coolers_to_combine
uris_IAA=( \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-IAA-si-ctrl-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-IAA-si-ctrl-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC442-IAA48H-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC442-IAA48H-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200805_G1sorted_442_replicate/results/coolers_library_group/G1-IAA.hg19.mapq_30.1000.mcool::/resolutions/1000" \
)

# HAP1-mother-cell
uris_WT=( \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-NT__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-IAA4H__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-IAA24H__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-IAA48H__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-NT-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-IAA4H-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-IAA48H-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-IAA4HW4-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-IAA4HW24-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200417_remap_polIIdegron/results/coolers_library/Hap1-WT-IAA4HW48-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
)

# this should work now ! assuming correct locations ...
cooler merge xxx.1kb.cool ${uris_NT[@]} # merge an array of cooler URIs
cooler zoomify -p 4 --balance -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000 -o xxx.mcool xxx.1kb.cool

# some relevant cooler commands as examples
# coarsing - i.e. turning high-res (small bin) into lower res (larger binsize) and copying it into existing MCOOL...
cooler coarsen -k 2 -p 3 -o all_combined.mcool::/resolutions/2000 all_combined.1kb.cool
cooler balance -p 4  all_combined.mcool::/resolutions/2000 # typical balance
cooler cp xxx.1kb.cool yyy.mcool::resolutions/1000 # copy standalone 1kb cooler into an mcool @resolutions/1000
