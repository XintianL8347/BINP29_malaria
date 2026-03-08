# Malaria Case Study

This project is created by Xintian Liu during the course BINP29 at Lund University. Following this project, we will examnine the phylogenetic relationships among 8 parasite species causing malaria in different hosts (6 Plasmodium species, Haemoproteus tartakovskyi, and Toxoplasma gondii). We wanted to investigate if Plasmodium falciparum is related to the other mammalian parasites, or if it originates from a bird malaria parasite which has changed its host. 

## Installation

All the tools used in this project were either already present in the Server or installed using conda. Example commands look like this:

```bash
conda create -n malaria
conda activate malaria
conda install -c bioconda proteinortho
conda install -c bioconda clustalo raxml
```

## Usage
1. Processing of Ht data

1.1 Clean genome sequence

All the analysis were carried out on the Server. Firstly, a filtration using a dowloaded python script was performed on the raw Ht genome file, with a maximum GC content of 30% and minimum scaffold length of 3000: 

```bash
python Scripts/removeScaffold.py Raw_Data/Haemoproteus_tartakovskyi.raw.genome 30 Haemoproteus_tartakovskyi.filter30.genome 3000
```
30% GC content was chosen as it falls in between the GC range of parasite (19-42%) while is distinct from bird GC (42%). A few other thresholds were tested (37%, 40%, 41%, etc.), and cutting at 41% and 30% didn't make dramatic differences in the number of scaffolds kept. This suggests that 30% GC cutoff returned a cleaner parasite set while keeping a sufficient number of representative sequences. However, there might still be bird genome remained in the file. 

1.2 Gene prediction

Next, we conducted a gene prediction using: 

```bash
nohup gmes_petap.pl --ES --min_contig 3000 --sequence ../Filtered/Haemoproteus_tartakovskyi.filter30.genome &
```
This output several files including genemark.gtf, which we renamed to Ht.gtf. Then we reformatted it and ran a downloaded Gff parser: 

```bash
cat Ht.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht2.gtf
gffParse.pl -i ../Filtered/Haemoproteus_tartakovskyi.filter30.genome -g Ht2.gtf -f gene -c -p
```
We then performed a blastp search on the parsed protein file, to identify the remaining bird scaffolds after the filtering step: 

```bash
blastp -query Prediction/gffParse.faa -db SwissProt -evalue 1e-10 -out Blast/Ht.blastp -num_descriptions 10 -num_alignments 5 -num_threads 10 
```
After getting the BLAST result, we parsed the result using downloaded script and dat files. Then we removed the bird sequences based on their scaffold IDs using bash. First grabbed the headers and reads into a new file: 

```bash
python datParser.py Ht.blastx Ht.fna taxonomy.dat uniprot_sprot.dat > scaffolds.txt
grep -f scaffolds.txt -A1 --no-group-separator Filtered/Haemoproteus_tartakovskyi.filter30.genome > scaffoldreads.txt
```
And then removed both headers and reads and wrote the cleaned results into a new file: 

```bash
grep -vf scaffoldreads.txt Filtered/Haemoproteus_tartakovskyi.filter30.genome > Ht_no_birds.fna
```
Following all that, a second prediction was run on the no-bird Ht fasta file: 

```bash
nohup gmes_petap.pl --ES --min_contig 3000 --cores 10 --sequence ../Ht_no_birds.fna &
```

2. Phylogenetic trees

After processing Ht data, we ran gffParse.pl on the gff/gtf files from all species. Some gff/gtf files needed extra formatting to be compatable with the script, such as replacing tabs with spaces. After all the formatting, example commands look like: 

```bash
cat Plasmodium_berghei.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Pb.gtf  # depending on the original fomat of gff/gtf files, this step is optional 
gffParse.pl -i ../Raw_Data/Plasmodium_berghei.genome -g Pb.gtf -f gene -c -p -b Pb    # if want to regenerate outfiles using the same basename, a flag -F was needed
```
We added corresponding group names in the headers in the generated faa files. When running ProteinOrtho for the first time, it complained about unrecognizable symbols, therefore a cleanning step was run on all species' faa files, such as: 

```bash
nohup proteinortho6.pl ../prePhylo/{Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa &  # didn't work
sed -i -E '/^>/! s/[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]//g; /^$/d' ../prePhylo/Tg.faa
```
Then reran the ProteinOrtho: 
```bash
nohup proteinortho6.pl -force -cpus=100 ./{Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa &
```
Next, the output file "myproject.proteinortho.tsv" was examined, and only the orthologs with one gene per species were exracted to a new file: 
```bash
head -n 1 myproject.proteinortho.tsv > filtered.tsv
grep "^8  8" myproject.proteinortho.tsv >> filtered88.tsv   # yielded 9 ortholog groups
# the following command can be run instead if one wants to also allow orthologs with one extra gene in one of the species to be used in the following tree building. (more data, more complex processings)
grep "^8  [89]" myproject.proteinortho.tsv >> filtered89.tsv 
```
Aftering filtering, we ran another proteinortho script "proteinortho_grab_protein.pl" to get the protein sequences according to the ortholog IDs: 

```bash
proteinortho_grab_protein.pl -tofiles filtered88.tsv ../{Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa
```
A script were written to loop over all the ortholog groups and build one alignment for each file. After alignment files where generated, we replaced the header again with only group names. For example: 

```bash
bash ~/Malaria/Scripts/clustalo.sh  # run clustalo on all the ortho groups
sed -E 's/[0-9]+_g//' filtered88.tsv.OrthoGroup0.fasta_aligned.faa > OrthoGroup0_aligned_cleaned.faa    # replace headers (there are of course smarter ways of doing it or possibilities of doing it in between other steps, but here I just did it in the hard way as I was only working with 9 ortho files)
```
Similarly, we wrote another script to loop over all the aligned fasta and build one tree per alignmnet: 

```bash
bash ~/Malaria/Scripts/raxml.sh
```
Finally, we concatenated all the single-ortholog trees into one file called "intree", and merged the trees using consense from the phylip package: 

```bash
cat RAxML_result.*.tre > intree
consense
```
The output tree was then visualized at: at: http://itol.embl.de/

A BUSCO analysis was also run for each sample, but no further steps (i.e., writing a script to compare busco ids) were carried out due to time limit. 
```bash
busco -i ../prePhylo/Pb.faa -o Pb -m prot -l apicomplexa 
```
