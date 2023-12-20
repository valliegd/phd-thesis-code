INPUT_LIST=list_GangSTR_vcf.txt
OUTPUT_FILE=GangSTR_output_final.csv

cat $INPUT_LIST | while read line; do

    IFS=';/,' read -ra NAMES <<< "$line"
    INPUT_FILE=${NAMES[0]}
    echo ${INPUT_FILE}
    PLATEKEY=`grep '^#CHROM' ${INPUT_FILE} | cut -f10`
    echo ${PLATEKEY}
    RESULTS=`grep -v '^#' ${INPUT_FILE} | cut -f10 | cut -d ':' -f1,2,3,4,5`
    echo ${RESULTS}

    echo "${PLATEKEY}, ${RESULTS}" >> ${OUTPUT_FILE}

done
