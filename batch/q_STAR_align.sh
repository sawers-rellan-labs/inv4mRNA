#!/usr/bin/csh

# Set the path to the sample list file
set sampleList=$1

  # Loop over the sample names in the sample list file
  foreach sample (`cat ${sampleList}`)
  
  # Construct the bsub command to run the alignment job for the current sample
    # set par="-q sara -n 12 -R 'span[hosts=1]' -R 'rusage[mem=32GB]' -W 4:00 -o %J.stdout -e %J.stderr" 
    set  par="-q sara -n 8 -W 4:00 -o %J.stdout -e %J.stderr" 
    echo $sample 
    echo $par
    # Submit the bsub command to the job scheduler
    bsub -R 'span[hosts=1]' -R 'rusage[mem=24GB]' $par ./align_sample.sh ${sample}
  end
