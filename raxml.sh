for file in OrthoGroup*_aligned_cleaned.faa; do
  out=$(sed "s/_aligned_cleaned.faa//" <<< $file); raxmlHPC -s $file -n $out.tre -o Tg -m PROTGAMMABLOSUM62 -p 12345; 
done
