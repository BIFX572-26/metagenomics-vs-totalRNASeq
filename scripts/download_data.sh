# download data associated with this analysis
# run in the root directory of the repository (after linking data/ to where you want the data downloaded)

TotalRNA_IDs=("SRR24451015" "SRR24451017" "SRR24451016")
MetaGen_IDs=("SRR24451016" "SRR24451018" "SRR24451020")

cd data
mkdir -p totalRNA
cd totalRNA

for item in "${TotalRNA_IDs[@]}"; do
  prefetch $item
  fasterq-dump $item --split-files -e 8 -p
  gzip $item
done

cd ..
mkdir -p MetaGen
cd MetaGen

for item in "${MetaGen_IDs[@]}"; do
  prefetch $item
  fasterq-dump $item --split-files -e 8 -p
  gzip $item
done

cd ../..