#/usr/bin/csh

set STAR=/usr/local/usrapps/maize/STAR-2.7.10b/bin/Linux_x86_64_static/STAR
set ref_idx_dir=
set ref_fatsa_dir=
set GTF_file=

$STAR 

$STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir $idx_dir \
--genomeFastaFiles $ref_fatsa_dir \
--sjdbGTFfile $GTF_file \
--sjdbOverhang 99
