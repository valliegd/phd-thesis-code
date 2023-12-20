module load bio/illumina/ExpansionHunter/3.2.2

# Please ensure there are no headers in the FILE

cat all_families_FP_38.txt | while read i; do 
	ExpansionHunter --reads $i  --reference /public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa --variant-catalog /re_gecip/neurology/Valentina/Running_EHv3/variant_catalog_GRCh38_ATXN3.json --output-prefix /re_gecip/shared_allGeCIPs/AD_VGD/ATXN3_EHv3/hg38/$(basename ${i} .bam)"_ATXN3" ;
done
