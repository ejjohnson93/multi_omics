#!/bin/bash 

HISAT2_INDEX=/path/goes/here/hg19_hisat2_index/hg19
SALMON_INDEX=/path/goes/here/hg19/hg19_salmon_index/hg19

# Set hard file path for files to align if required
DATA_IN=/path/goes/here/ageing_datasets/neurons_martens_PRJEB44542/
cd ${DATA_IN}

# Create directory for output data
HISAT2_OUTPUT=${DATA_IN}/hisat2
SALMON_OUTPUT=${DATA_IN}/salmon_quant

mkdir -p ${HISAT2_OUTPUT} ${SALMON_OUTPUT}

# Insert a .txt or .tsv file here with the IDs of the .fastq files to align
# Note: if script isn't working make sure there aren't dos linebreaks in the input file

for i in $(cut -f 1 just_RNA_seq.tsv); do

	echo "${i}"
 
  # Alignment, HISAT2 and salmon
  
  # If the current ID has a file that ends with _2 then the reads are paired...
  if [ -e ${i}_2.fastq.gz ]; then
      echo "${i} is paired-end, beginning alignments with HISAT2 and salmon..."
      echo
      
      # hisat2 alignment
      
      echo "Aligning reads for ${i} with HISAT2..."
      hisat2 -x ${HISAT2_INDEX} \
          -p $SLURM_CPUS_PER_TASK \
          -1 "${DATA_IN}/${i}_1.fastq.gz" \
          -2 "${DATA_IN}/${i}_2.fastq.gz" \
          | samtools view -b | samtools sort > ${HISAT2_OUTPUT}/${i}_hisat2.bam
  
      
      # salmon alignment 
  
      echo "Salmon quasi-alignment for ${i}..."
          salmon quant -i ${SALMON_INDEX} \
              --validateMappings \
              -l A \
              -1 "${DATA_IN}/${i}_1.fastq.gz" \
              -2 "${DATA_IN}/${i}_2.fastq.gz" \
              -o "${SALMON_OUTPUT}/${i}_quant" \
              -p 4
  
      
  # If the current ID doesn't have a file ending with _2 but does have a file ending with _1 then the reads are single read...
  
  elif [ -e ${i}_1.fastq.gz ]; then
      echo "${i} is single-read, beginning alignment with HISAT2 and salmon..."
      echo
      
      # hisat2 alignment
  
      echo "Aligning reads for ${i} with HISAT2..."
      hisat2 -x ${HISAT2_INDEX} \
          -p $SLURM_CPUS_PER_TASK \
          -U "${DATA_IN}/${i}_1.fastq.gz" \
          | samtools view -b | samtools sort > ${HISAT2_OUTPUT}/${i}_hisat2.bam
  
      
      # salmon alignment 
  
      echo "Salmon quasi-alignment for ${i}..."
          salmon quant -i ${SALMON_INDEX} \
              --validateMappings \
              -l A \
              -r "${DATA_IN}/${i}_1.fastq.gz" \
              -o "${SALMON_OUTPUT}/${i}_quant" \
              -p 4
  
  else 
  
  # No files ending in _1 or _2 suggest the files are incorrectly named or the target directory is wrong, script determines no reads available..
  
      echo "No reads available..."
  
  fi 

 
done

