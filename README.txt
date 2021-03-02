##What do we miss for the paper?

R-loop boundaries NO CTCF NO G4 do not show cohesin accumulation. Are they compartment boundaries? or are they boundaries that do not show RAD21 accumulation?

## BED File locations

# The different categories used in the paper
R-loops NO CTCF G4
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/diff-cat/intersect-Rloops-NO-all-NT-CTCF-G4-Mao
R-loops NO CTCF NO G4 (####### THIS IS THE CATEGORY WE WOULD KNOW IF THEY ARE COMPARTMENTS BOUNDARIES)
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/diff-cat/intersect-Rloops-NO-all-NT-CTCF-NO-G4-Mao
CTCF NO R-loops G4
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/diff-cat/intersect-all-NT-CTCF-NO-Rloops-G4-Mao
CTCF NO R-loops NO G4
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/diff-cat/intersect-all-NT-CTCF-NO-Rloops-NO-G4-Mao
CTCF R-loops G4
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/diff-cat/intersect-all-NT-CTCF-Rloops-G4-Mao
CTCF R-loops NO G4
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/diff-cat/intersect-all-NT-CTCF-Rloops-NO-G4-Mao

# ALL R-loops
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/GSM1720619_K562_DRIP_peaks-sort.bed

# ALL G4
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/G4-Mao-sort

# ChIP seq peaks
CTCF peaks
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/NT-CTCF-narrowPeaks-sort-merge
DDX55 peaks
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/NT-DDX55-narrowPeaks-sort-merge
RAD21 peaks
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/NT-RAD21-narrowPeaks-sort-merge

# splicing events
# In CTCF degron after CTCF depletion
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/S44_NT_S442_IAA.output_events_all_select_hg19.txt
# In DDX55 clone1
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/AAVS1_sg24_NT_DDX55_sg2B_NT.output_events_all_select_hg19.txt

## coolers
#CTCF degron asynchronous
NT
/nl/umw_job_dekker/users/av90w/distiller-run/20200325_mapping_for_dot_calling_cispercent/results/coolers_library/CkoC442-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool
IAA
/nl/umw_job_dekker/users/av90w/distiller-run/20200325_mapping_for_dot_calling_cispercent/results/coolers_library/CkoC442-IAA48H-R1-T1__hg19.hg19.mapq_30.1000.mcool

#DDX55 clones
DDX55sg2B NT
/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-DDX55sg2-B-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool
DDX55sg2B IAA
/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-DDX55sg2-B-IAA-R1-T1__hg19.hg19.mapq_30.1000.mcool
DDX55sg27 NT
/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-DDX55sg2-7-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool
DDX55sg27 IAA
/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-DDX55sg2-7-IAA-R1-T1__hg19.hg19.mapq_30.1000.mcool

#AAVS1 clone
NT
/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-AAVS1sg2-4-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool
IAA
/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library/CkoCT442-AAVS1sg2-4-IAA-R1-T1__hg19.hg19.mapq_30.1000.mcool

# libraries to merge together to call compartments at a higher resolution (10 libraries together)
/nl/umw_job_dekker/users/av90w/distiller-run/20200220_mapping_for_dot_calling/results/coolers_library
CkoC44-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool
CkoC44-NT-R2-T1__hg19.hg19.mapq_30.1000.mcool
CkoC44-IAA48H-R1-T1__hg19.hg19.mapq_30.1000.mcool
CkoC44-IAA48H-R2-T1__hg19.hg19.mapq_30.1000.mcool
CkoC442-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool
CkoC442-NT-R2-T1__hg19.hg19.mapq_30.1000.mcool
CkoCT442-G1-NT-R1-T1__hg19.hg19.mapq_30.1000.mcool
CkoCT442-NT-si-ctrl-R1-T1__hg19.hg19.mapq_30.1000.mcool
CkoCT442-NT-si-ctrl-R2-T1__hg19.hg19.mapq_30.1000.mcool
/nl/umw_job_dekker/users/av90w/distiller-run/20200805_G1sorted_442_replicate/results/coolers_library
CkoCT442-G1-NT-R2-T1__hg19.hg19.mapq_30.1000.mcool

## compartment files

# CTCF degron non treated
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/CTCFdegron-442-nontreated-rep1.250kb.eigs.cis.vecs.txt
#AAVS1 clone
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/AAVS1sg24-nontreated.250kb.eigs.cis.vecs.txt
# DDX55sg2B clone
/nl/umw_job_dekker/users/av90w/projects/for-paper/for-higlass/DDX55sg2B_NT-nontreated.250kb.eigs.cis.vecs.txt

## insulation files

#for asynch CTCF degron, siRNA and clones (on shadow)
/data/alv/CTCF_degron/analysis/siRNA/insulation

## ChIP seq

#RAD21 bigwig
rep1
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/test-ChIP2/new-mapping/results/bwa/mergedLibrary/bigwig/NT-RAD21_R1.mLb.clN.bigWig
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/test-ChIP2/new-mapping/results/bwa/mergedLibrary/bigwig/IAA-RAD21_R1.mLb.clN.bigWig
rep2
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP4-bis/results/bwa/mergedLibrary/bigwig/NT-RAD21_R2.mLb.clN.bigWig
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP3/results/bwa/mergedLibrary/bigwig/IAA-RAD21_R1.mLb.clN.bigWig

#CTCF bigWig
rep1
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP3/results/bwa/mergedLibrary/bigwig/NT-CTCF_R1.mLb.clN.bigWig
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP3/results/bwa/mergedLibrary/bigwig/IAA-CTCF_R1.mLb.clN.bigWig
rep2
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP4-bis/results/bwa/mergedLibrary/bigwig/NT-CTCF_R2.mLb.clN.bigWig
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP4-bis/results/bwa/mergedLibrary/bigwig/IAA-CTCF_R2.mLb.clN.bigWig

#DDX55 bigWig
rep1
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/test-ChIP2/new-mapping/results/bwa/mergedLibrary/bigwig/NT-DDX55_R1.mLb.clN.bigWig
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/test-ChIP2/new-mapping/results/bwa/mergedLibrary/bigwig/IAA-DDX55_R1.mLb.clN.bigWig
rep2
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP3/results/bwa/mergedLibrary/bigwig/NT-DDX55_R1.mLb.clN.bigWig
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP3/results/bwa/mergedLibrary/bigwig/IAA-DDX55_R1.mLb.clN.bigWig

#input control bigwig
#input controls for IAA RAD21 R2, CTCF R1, DDX55 R2
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP3/results/bwa/mergedLibrary/bigwig/NT-input_R1.mLb.clN.bigWig 
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP3/results/bwa/mergedLibrary/bigwig/IAA-input_R1.mLb.clN.bigWig 

#input controls for NT RAD21 R2, CTCF R2
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP4-bis/results/bwa/mergedLibrary/bigwig/NT-input_R2.mLb.clN.bigWig
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/ChIP4-bis/results/bwa/mergedLibrary/bigwig/IAA-input_R2.mLb.clN.bigWig

#input controls for RAD21 R1, DDX55 R1
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/test-ChIP2/new-mapping/results/bwa/mergedLibrary/bigwig/NT-input_R1.mLb.clN.bigWig
/nl/umw_job_dekker/users/av90w/ChIPseq/fastq/test-ChIP2/new-mapping/results/bwa/mergedLibrary/bigwig/IAA-input_R1.mLb.clN.bigWig

## RNA seq
/nl/umw_job_dekker/users/av90w/RNAseq/data/siRNA/report2546/bigwig_star/S442_NT.bw

