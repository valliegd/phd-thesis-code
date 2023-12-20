module load vis/matplotlib/3.0.2-foss-2018b-Python-3.6.6
module load lang/Python/3.7.2-GCCcore-8.2.0
module load bio/Pysam/0.15.1-foss-2018b-Python-3.6.6

python3 /path_to_alignmentviewer/GraphAlignmentViewer/GraphAlignmentViewer.py --variant_catalog /path_to_json/variant_catalog_GRCh37_HTT_allconfig.json --read_align /path_to_bamfile/x.bam --output_prefix /path_to_output_folder/file_output --reference_fasta /path_to_reference/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --locus_id "HTT" --file_format v3
