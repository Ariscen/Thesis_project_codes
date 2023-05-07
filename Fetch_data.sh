
# GenotypeSamples files
cd ~/projects/sceQTLsim/data/onek1k/GenotypeSamples

for ((i=1; i<=39; i ++))
do
  echo https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5899nnn/GSM5899`expr $i + 872`/suppl/GSM5899`expr $i + 872`_OneK1K_scRNA_Sample${i}_GenotypeSamples.txt.gz >> GenotypeSamples_urls.txt
done
for ((i=40; i<=64; i ++))
do
  echo https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5899nnn/GSM5899`expr $i + 872`/suppl/GSM5899`expr $i + 872`_OneK1K_scRNA_Sample`expr $i + 1`_GenotypeSamples.txt.gz >> GenotypeSamples_urls.txt
done
for ((i=65; i<=75; i ++))
do
  echo https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5899nnn/GSM5899`expr $i + 872`/suppl/GSM5899`expr $i + 872`_OneK1K_scRNA_Sample`expr $i + 2`_GenotypeSamples.txt.gz >> GenotypeSamples_urls.txt
done

wget -i GenotypeSamples_urls.txt #75 files in total

# Individual Barcodes files
cd ~/projects/sceQTLsim/data/onek1k/Individual_Barcodes

for ((i=1; i<=39; i ++))
do
  echo https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5899nnn/GSM5899`expr $i + 872`/suppl/GSM5899`expr $i + 872`_OneK1K_scRNA_Sample${i}_Individual_Barcodes.csv.gz >> Individual_Barcodes_urls.txt
done
for ((i=40; i<=64; i ++))
do
  echo https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5899nnn/GSM5899`expr $i + 872`/suppl/GSM5899`expr $i + 872`_OneK1K_scRNA_Sample`expr $i + 1`_Individual_Barcodes.csv.gz >> Individual_Barcodes_urls.txt
done
for ((i=65; i<=75; i ++))
do
  echo https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5899nnn/GSM5899`expr $i + 872`/suppl/GSM5899`expr $i + 872`_OneK1K_scRNA_Sample`expr $i + 2`_Individual_Barcodes.csv.gz >> Individual_Barcodes_urls.txt
done

wget -i Individual_Barcodes_urls.txt #75 files in total

# RawCounts files
cd ~/projects/sceQTLsim/data/onek1k/RawCounts

for ((i=1; i<=39; i ++))
do
  echo https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5899nnn/GSM5899`expr $i + 872`/suppl/GSM5899`expr $i + 872`_OneK1K_Pool${i}_RawCounts.csv.gz >> RawCounts_urls.txt
done
for ((i=40; i<=64; i ++))
do
  echo https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5899nnn/GSM5899`expr $i + 872`/suppl/GSM5899`expr $i + 872`_OneK1K_Pool`expr $i + 1`_RawCounts.csv.gz >> RawCounts_urls.txt
done
for ((i=65; i<=75; i ++))
do
  echo https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5899nnn/GSM5899`expr $i + 872`/suppl/GSM5899`expr $i + 872`_OneK1K_Pool`expr $i + 2`_RawCounts.csv.gz >> RawCounts_urls.txt
done

wget -i RawCounts_urls.txt #75 files in total

