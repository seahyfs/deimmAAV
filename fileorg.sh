#!/bin/bash

master_folder="/mnt/volume2/data/Cas_only/trial"  # Replace with the actual path to the master folder
consolidated_folder="${master_folder}/consolidated"  # Replace with the name of the consolidated folder

# Create the consolidated folder if it doesn't exist
mkdir -p "$consolidated_folder"

# Iterate over all directories within the master folder
for directory in "${master_folder}"/*/; do
  if [[ -d "$directory" ]]; then  # Check if it is a directory

    # Extract the directory name
    dir_name=$(basename "$directory")

    # Change to the directory
    cd "$directory" || { echo "Failed to change directory: $directory"; exit 1; }

     # Run the zcat command and redirect output to the corresponding output files
    zcat *_1.fq.gz > "${dir_name}_R1.fastq"
    zcat *_2.fq.gz > "${dir_name}_R2.fastq"

	#zip everything
	gzip *_R1.fastq
	gzip *_R2.fastq

	#move to consolidated folder
	mv "${dir_name}_R1.fastq.gz" "$consolidated_folder/"
	mv "${dir_name}_R2.fastq.gz" "$consolidated_folder/"

    # Display a message indicating the completion of the task
    echo "Processed directory ${directory}"
  fi
done

# Display a message indicating the completion of the task
echo "Consolidation completed."
