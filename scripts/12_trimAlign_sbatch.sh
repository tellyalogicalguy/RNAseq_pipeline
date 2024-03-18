#!/bin/bash
#SBATCH --account=rpp-jfcote11
#SBATCH --mail-user=poorani.subramani@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4000
#SBATCH --time=0:30:0
#SBATCH --tmp=100G
#SBATCH --output="./slurm_outputs/%x_slurm-%j.out"

set -e
set -u
set -o pipefail

SLRM_LOCAL_TMP_DIR=$SLURM_TMPDIR/$RANDOM/${SAMPLE_ID}/
mkdir -p ${SLRM_LOCAL_TMP_DIR}

#---------------------------------------------------------#
RAW_FASTQ_FILES=($(awk -v S_ID="${SAMPLE_ID}" '$1==S_ID {print $3}' ${DATA_INFO_FILE}))
declare -a RAW_FASTQ
RAW_FASTQ[0]=$(find ${RAW_FASTQ_DIR} -name "${RAW_FASTQ_FILES[0]}")
RAW_FASTQ[1]=$(find ${RAW_FASTQ_DIR} -name "${RAW_FASTQ_FILES[1]}")

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 1. Starting FASTQ filtering with fastp"
module load fastp/0.23.1

FILT_FASTQ_R1=${SLRM_LOCAL_TMP_DIR}/${SAMPLE_ID}_R1.filteredFastq.gz
FILT_FASTQ_R2=${SLRM_LOCAL_TMP_DIR}/${SAMPLE_ID}_R2.filteredFastq.gz

fastp --in1 ${RAW_FASTQ[0]} --in2 ${RAW_FASTQ[1]} \
--out1 $FILT_FASTQ_R1 --out2 $FILT_FASTQ_R2 \
--detect_adapter_for_pe --thread 8 \
--html ${FILT_FASTQ_DIR}/${SAMPLE_ID}.html --json /dev/null

cp $FILT_FASTQ_R1 ${FILT_FASTQ_DIR}/
cp $FILT_FASTQ_R2 ${FILT_FASTQ_DIR}/

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 1. Done FASTQ filtering."
#---------------------------------------------------------#
echo $(date -Iseconds)" SPG_RNAseq_pipeline: 2. Starting aligning with STAR"

module load star/2.7.9a

STAR --runThreadN "$SLURM_CPUS_PER_TASK" \
	--genomeDir "$GENOME_DIR" \
	--readFilesIn $FILT_FASTQ_R1 $FILT_FASTQ_R2 \
	--readFilesCommand gunzip -c \
	--outSAMtype BAM Unsorted \
	--outFileNamePrefix "${ALIGNED_BAM_DIR}/${SAMPLE_ID}_"

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 2. Done STAR alignment."
#---------------------------------------------------------#
# if spike-in genome present, align to spike-in genome
if [[ $SPIKEIN_GENOME != "none" ]]
then

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 2. Starting aligning to $SPIKEIN_GENOME genome with STAR"

# "--seedPerWindowNmax 10" option is necessary for spike-in because this determines the amount of time STAR looks for better matches. The default is 50 and becuase the % of mapped reads is low for splike-in, a lot of reads spend a lot of time looking for better alignment.
STAR --runThreadN "$SLURM_CPUS_PER_TASK" \
        --genomeDir "$SPIKEIN_GENOME_DIR" \
        --readFilesIn $FILT_FASTQ_R1 $FILT_FASTQ_R2 \
        --readFilesCommand gunzip -c \
        --outSAMtype BAM SortedByCoordinate \
	--seedPerWindowNmax 10 \
        --outFileNamePrefix "${ALIGNED_BAM_DIR}/${SAMPLE_ID}_${SPIKEIN_GENOME}_"

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 2. Done STAR alignment to $SPIKEIN_GENOME genome."

fi
#---------------------------------------------------------#
#making symlinks
for SYML_DIR in FILT_FASTQ_SYML_DIR ALIGNED_BAM_SYML_DIR
do
	SCRATCH_DIR_VARNAME="$(sed "s:_SYML::" <(echo ${SYML_DIR}))"
	eval SCRATCH_DIR='$'"${SCRATCH_DIR_VARNAME}"
	echo $SCRATCH_DIR_VARNAME" = "$SCRATCH_DIR
		for file in $(ls ${SCRATCH_DIR})
		do
		if [ ! -L "${!SYML_DIR}${file}" ]; then
			ln -s ${SCRATCH_DIR}${file} ${!SYML_DIR}${file}
		fi
	done

find ${!SYML_DIR} -xtype l -delete

done
#---------------------------------------------------------#
# update progress file and start the next step in the pipeline if appropriate
echo $(date -Iseconds)" SPG_RNAseq_pipeline: Updating BAM progress file"


exec 200>${SCRIPT_CUR_DIR}/slurm_outputs/.lock.file || exit 1
# set flock on the BAM_PROG_FILE ?
flock 200 || exit 1

gawk -i inplace -v S_ID="${SAMPLE_ID}" 'BEGIN{OFS="\t"} $1==S_ID {$2="Done"}'1 ${BAM_PROG_FILE}

cat ${BAM_PROG_FILE}

while IFS=$'\t' read -r s_id status; do
        if [[ $status != "Done" ]]
        then
                echo "Waiting for "$s_id
                exit 0
        fi
done < ${BAM_PROG_FILE}

flock -u 200

#---------------------------------------------------------#
# if spike-in, do deseq2 normalization
if [[ $SPIKEIN_GENOME != "none" ]]
then
	module load gcc/9.3.0
	module load r/4.2.1

	echo $(date -Iseconds)" SPG_RNAseq_pipeline: 2. Getting normFactors from spike-in DESeq2.."
	Rscript ${SCRIPT_CUR_DIR}/getNormFactorsFromSpikeIn.R
	echo $(date -Iseconds)" SPG_RNAseq_pipeline: 2. Done getting normFactors."
fi
#---------------------------------------------------------#
echo $(date -Iseconds)" SPG_RNAseq_pipeline: Starting second pass with STAR"
#for each read pair, send the variables to sbatch

SAMPLE_IDS=($(awk 'NR>1 {print $1}' ${DATA_INFO_FILE} | uniq))

for SAMPLE_ID in "${SAMPLE_IDS[@]}"
do
        echo "Submitting: "$SAMPLE_ID
        export SAMPLE_ID
        SLRM_NAME="34_star2ndPassBW_${SAMPLE_ID}"
	echo "Starting paired-end alignment procedure..."
	sbatch --job-name "$SLRM_NAME" 34_star2ndPassBW_sbatch.sh
        sleep 0.2
done

