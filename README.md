
### 1. Fully mask genome
#### a. customize a repeat library if do not have one
$ BuildDatabase -name dbname -engine wublast  refSeq.fasta

$ RepeatModeler -pa 32 -engine wublast -database dbname

#### b. Identify segmental duplications in reference genome
$ WGAC_SD.pl -i refSeq.fasta -l repeat.lib.fasta -t 48

#### c. thoroughly mask reference genome and get a unique genome

$ thorough_mask.pl -i ref.fasta -r TEseq.fasta -t 48 -s SD_interval.bed

### 2. Mapping reads to unique genome

mrfast

### 3. calculate and normalize read depth for each sample


### 4. call CNVs and calulate gene copy number

### 5. compare allele frequency
