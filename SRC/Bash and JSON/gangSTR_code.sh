module av
module load GangSTR/2.4.2

GangSTR --bam ~/path_to_bamfile/x.bam --ref ~/path_to_reference/GRCh38/GRCh38Decoy_no_alt.fa --regions ~/path_to_genomic_regions/GangSTR/GangSTR_regions/HTT_GRCh38.bed --out ~/path_to_output_folder/file_output

Compress vcf > vcf.gz
bgzip ~/re_gecip/neurology/Valentina/GangSTR/GangSTR_output/x.vcf

Decompress vcf.gz > vcf
bgzip -d

Extract GangSTR info from .vcf files
for r in *.vcf; do grep cag $r; done

Save as txt/csv
for r in *.vcf; do grep cag $r &>> GangSTR_HTT_output.txt; done

Save as txt/csv output + name file
for r in *.vcf; do grep 'cag\|LP' $r &>> GangSTR_HTT_output_alldata.txt; done

Extract LP only
grep -o "LP.......-DNA_..." GangSTR_HTT_output_alldata.csv | sort --unique

List all VCF files
find *.vcf>list_GangSTR_vcf.txt

Cut LP and data
grep '^#CHROM' GangSTR_HTT_output_alldata.vcf | cut -f10
grep -v '^#' GangSTR_HTT_output_alldata.vcf | cut -f10 | cut -d':' -f1,2,3,4

Change permission of file
chmod +rwx retriving_fields_from_vcf.sh (r = read, w = write, x = execute)

Execute sh file bash
./file.sh
