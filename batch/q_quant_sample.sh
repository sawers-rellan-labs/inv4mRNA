#!/usr/bin/csh

# Set the path to the sample list file
set sampleList=$1
set index=Zm-B73-REFERENCE-NAM-5.0.cdna.all.idx
set transcripts=/rsstu/users/r/rrellan/sara/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa

  # Loop over the sample names in the sample list file
  foreach sample (`cat ${sampleList}`)
  
  # Construct the bsub command to run the alignment job for the current sample
    # set par="-q sara -n 12 -R 'span[hosts=1]' -R 'rusage[mem=32GB]' -W 4:00 -o %J.stdout -e %J.stderr" 
    set  par="-n 1 -W 4:00 -o %J.stdout -e %J.stderr" 
    echo $sample 
    echo $par
    # Submit the bsub command to the job scheduler
    bsub -R 'rusage[mem=2GB]' ${par} ./quant_sample.sh ${sample}
  end
