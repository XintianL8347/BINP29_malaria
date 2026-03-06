for file in *OrthoGroup*.fasta; do
  clustalo -i $file -o ../Align/$(basename $file)_aligned.faa -v; 
done
