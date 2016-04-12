DATA_FOLDER=$1
CASE=$2
OUTPUT_FILE=$3

DATA_DIR=$(dirname `pwd`)"/results/${DATA_FOLDER}"
RESULT_FILES=${DATA_DIR}"/"${CASE}"_run_*.out"

# get the position of the file-path name where the run number starts
PREFIX_SIZE=${#DATA_DIR}+${#CASE}
PREFIX_SIZE=$(($PREFIX_SIZE + 7))
echo ${DATA_DIR}

#PREFIX_SIZE2=$(expr length "--executor-cores x")

# summarize results for case with varying number of cores (fixed blocking factor)
echo ""
echo "run, cores, blocking" > $CASE
grep 'Command:' $RESULT_FILES | awk -F' ' -v var1=$PREFIX_SIZE '{print substr($1, var1, 2)", "$4", "$12}' | sort >> $CASE

# consolidate results in a CSV file
join <(sort $CASE) <(grep 'rank: 0, SVD time:' $RESULT_FILES | awk -F': ' -v var=$PREFIX_SIZE '{print substr($1, var, 2) ", " $3}' | sort) -t$', ' | sort -t $',' -k 1,1 -V > $OUTPUT_FILE
