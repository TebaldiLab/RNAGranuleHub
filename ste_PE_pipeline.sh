# rename the files if not _R1_001.fastq and _R2_001.fastq using the ste_PE_rename.sh script

#!/usr/bin/env bash
set -euo pipefail

# Top-level of the project
ROOT_DIR="/data/stefan/stress_granules/paired/test2"
GENOME_INDEX="/data/stefan/genome/hs_STAR_index_v48"
TRIMMOMATIC_JAR="/data/stefan/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="/data/stefan/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"

# Function to rename FASTQ files to standard pattern
rename_fastq_files() {
    local sample_dir="$1"
    local sample_name="$2"
    
    # Try to find R1/R2 files with common patterns
    local r1_patterns=("_R1_" "_1_" "_f1" "_forward" "_read1" "_1")
    local r2_patterns=("_R2_" "_2_" "_r2" "_reverse" "_read2" "_2")
    
    for r1_file in "$sample_dir"/*; do
        [[ -f "$r1_file" ]] || continue
        
        # Check for R1 patterns
        for pattern in "${r1_patterns[@]}"; do
            if [[ "$r1_file" == *"$pattern"* ]]; then
                local r2_file="${r1_file/$pattern/${r2_patterns[0]}}" # Replace with first R2 pattern
                r2_file="${r2_file%.*}_001.fastq" # Ensure consistent suffix
                
                # If R2 exists with any pattern, rename both
                for r2_candidate in "$sample_dir"/*; do
                    for r2pat in "${r2_patterns[@]}"; do
                        if [[ "$r2_candidate" == *"$r2pat"* ]] && 
                           [[ "${r2_candidate/$r2pat/}" == "${r1_file/$pattern/}" ]]; then
                            # Rename files to standard pattern
                            local new_r1="${sample_dir}/${sample_name}_R1_001.fastq"
                            local new_r2="${sample_dir}/${sample_name}_R2_001.fastq"
                            
                            echo "    Renaming files to standard pattern:"
                            echo "      $r1_file -> $new_r1"
                            echo "      $r2_candidate -> $new_r2"
                            
                            mv "$r1_file" "$new_r1"
                            mv "$r2_candidate" "$new_r2"
                            return 0
                        fi
                    done
                done
            fi
        done
    done
    
    return 1
}

# Loop through all sample folders
for dataset in "$ROOT_DIR"/*_dataset; do
  for experiment in "$dataset"/*; do
    [[ -d "$experiment" ]] || continue
    echo "Processing experiment: $experiment"
    
    # Create BAM folder 
    BAM_DIR="$experiment/bam"
    rm -rf "$BAM_DIR"
    mkdir -p "$BAM_DIR"

    for sample_dir in "$experiment"/*; do
      [[ -d "$sample_dir" ]] || continue
      [[ "$(basename "$sample_dir")" == "bam" ]] && continue  # Skip bam directory
      sample_name=$(basename "$sample_dir")
      echo "  Sample: $sample_name"

      # Prepare metrics file in BAM directory
      METRICS="$BAM_DIR/${sample_name}_metrics.tsv"
      echo -e "step\treads" > "$METRICS"
      # Raw read counts (using R1)
      R1_RAW_FILES=("$sample_dir"/*_R1_001.fastq)
      if [[ ! -f "${R1_RAW_FILES[0]}" ]]; then
        echo "    Counting raw reads"
      fi
     
      # First try standard filenames
      R1=("$sample_dir"/*_R1_001.fastq)
      R2=("$sample_dir"/*_R2_001.fastq)
      
      # If standard names not found, try renaming
      if [[ ! -f "${R1[0]}" || ! -f "${R2[0]}" ]]; then
        echo "    Standard R1/R2 files not found, attempting to rename..."
        if rename_fastq_files "$sample_dir" "$sample_name"; then
          R1=("$sample_dir"/*_R1_001.fastq)
          R2=("$sample_dir"/*_R2_001.fastq)
        else
          echo "    >> ERROR: Could not find or rename R1/R2 FASTQ files"
          continue
        fi
      fi

      # Count raw reads
      raw_lines=$(wc -l < "${R1[0]}")
      raw_reads=$(( raw_lines / 4 )) # divide by 4 because 1 read takes up 4 lines
      echo -e "raw_reads\t$raw_reads" >> "$METRICS"

      # For gzipped files, decompress first
      if [[ "${R1[0]}" == *.gz ]]; then
        echo "    Decompressing gzipped FASTQ files..."
        gunzip "${R1[0]}"
        gunzip "${R2[0]}"
        R1=("$sample_dir"/*_R1_001.fastq)
        R2=("$sample_dir"/*_R2_001.fastq)
      fi

      # Deduplicate with FastUniq
      R1_RAW="${sample_dir}/${sample_name}_dup_R1.fq"
      R2_RAW="${sample_dir}/${sample_name}_dup_R2.fq"
      cat "${R1[@]}" > "$R1_RAW"
      cat "${R2[@]}" > "$R2_RAW"

      LIST="$sample_dir/${sample_name}_fq_list.txt"
      echo -e "$R1_RAW\n$R2_RAW" > "$LIST"

      /storage/shared/fastuniq \
        -i "$LIST" \
        -t q \
        -o "$sample_dir/${sample_name}_dedup_R1.fq" \
        -p "$sample_dir/${sample_name}_dedup_R2.fq" \
        -c 0

      rm "$R1_RAW" "$R2_RAW" "$LIST"
      # Count deduplicated reads
      dedup_lines=$(wc -l < "$sample_dir/${sample_name}_dedup_R1.fq")
      dedup_reads=$(( dedup_lines / 4 ))
      echo -e "dedup_reads\t$dedup_reads" >> "$METRICS"

      # Compress for Trimmomatic
      gzip "$sample_dir/${sample_name}_dedup_R1.fq"
      gzip "$sample_dir/${sample_name}_dedup_R2.fq"

      # Adapter + quality trimming
      java -jar "$TRIMMOMATIC_JAR" PE -phred33 \
        "$sample_dir/${sample_name}_dedup_R1.fq.gz" "$sample_dir/${sample_name}_dedup_R2.fq.gz" \
        "$sample_dir/${sample_name}_trim_forward_paired.fq.gz" "$sample_dir/${sample_name}_trim_forward_unpaired.fq.gz" \
        "$sample_dir/${sample_name}_trim_reverse_paired.fq.gz" "$sample_dir/${sample_name}_trim_reverse_unpaired.fq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 \
        -threads 4
         
      # Merge unpaired
      cat "$sample_dir/${sample_name}_trim_forward_unpaired.fq.gz" \
          "$sample_dir/${sample_name}_trim_reverse_unpaired.fq.gz" \
        > "$sample_dir/${sample_name}_trim_unpaired.fq.gz"

      # Count trimmed reads
      trimmed_paired_lines=$(zcat "$sample_dir/${sample_name}_trim_forward_paired.fq.gz" | wc -l)
      trimmed_paired=$(( trimmed_paired_lines / 4 ))
      echo -e "trimmed_paired\t$trimmed_paired" >> "$METRICS"
      

      # STAR alignment — paired
      STAR --runThreadN 8 \
        --genomeDir "$GENOME_INDEX" \
        --readFilesIn "$sample_dir/${sample_name}_trim_forward_paired.fq.gz" "$sample_dir/${sample_name}_trim_reverse_paired.fq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$BAM_DIR/${sample_name}_paired_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFilterMultimapNmax 5000

      # STAR alignment — unpaired
      STAR --runThreadN 8 \
        --genomeDir "$GENOME_INDEX" \
        --readFilesIn "$sample_dir/${sample_name}_trim_unpaired.fq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$BAM_DIR/${sample_name}_unpaired_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFilterMultimapNmax 5000

      # Rename and index BAMs
      cd "$BAM_DIR"
      mv "${sample_name}_paired_Aligned.sortedByCoord.out.bam"   "${sample_name}_paired.bam"
      mv "${sample_name}_unpaired_Aligned.sortedByCoord.out.bam" "${sample_name}_unpaired.bam"

      samtools index "${sample_name}_paired.bam"
      samtools index "${sample_name}_unpaired.bam"
  
      # Count aligned paired reads
      paired_count=$(samtools view -c -F 0xD04 -f 0x40 "$BAM_DIR/${sample_name}_paired.bam")
      echo -e "aligned_paired\t$paired_count" >> "$METRICS"
      
      # Count aligned unpaired reads
      unpaired_count=$(samtools view -c -F 0xD04 "$BAM_DIR/${sample_name}_unpaired.bam")
      echo -e "aligned_unpaired\t$unpaired_count" >> "$METRICS"     

      # Cleanup intermediates
      rm -f "$sample_dir/${sample_name}_trim_forward_unpaired.fq.gz" \
             "$sample_dir/${sample_name}_trim_reverse_unpaired.fq.gz" 
      
      
      echo "    Metrics written to $METRICS"
      cd - > /dev/null
    done
  done
done

echo "All samples processed. Output BAMs, count files and Metrics are stored in each experiment's bam/ directory."