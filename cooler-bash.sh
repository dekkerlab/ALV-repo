
# request a long interactive job on the cluster
bsub -q interactive -W 04:00 -n 4 -R "rusage[mem=16000]" -R "span[hosts=1]" -Is bash

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


uris_IAA=( \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-IAA-si-ctrl-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-IAA-si-ctrl-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC442-IAA48H-R1-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoC442-IAA48H-R2-T1__hg19.hg19.mapq_30.1000.mcool::/resolutions/1000" \
"/nl/umw_job_dekker/users/av90w/distiller-run/20200805_G1sorted_442_replicate/results/coolers_library_group/G1-IAA.hg19.mapq_30.1000.mcool::/resolutions/1000" \
)

# this should work now ! assuming correct locations ...
cooler merge xxx.1kb.cool ${uris_NT[@]} # merge an array of cooler URIs
cooler zoomify -p 4 --balance -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000 -o xxx.mcool xxx.1kb.cool
