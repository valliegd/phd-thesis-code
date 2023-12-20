module load bio/illumina/ExpansionHunter/3.2.2

# Please ensure there are no headers in the FILE

cat all_families_FP_37.txt | while read i; do 
	ExpansionHunter --reads $i  --reference /public_data_resources/reference/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --variant-catalog /re_gecip/neurology/Valentina/Running_EHv3/variant_catalog_GRCh37_ATXN3.json --output-prefix /re_gecip/shared_allGeCIPs/AD_VGD/ATXN3_EHv3/hg37/$(basename ${i} .bam)"_ATXN3" ;
done
