MEM=40g

#create directory for logs (if it doesn't already exist)
LOG_DIR="$HOME/projects/multiseq_analysis/results/roc/logs/"
if [ ! -d "${LOG_DIR}" ]; then
    mkdir ${LOG_DIR}
fi


for i in `echo "7"`; do
	PROCESS_NAME=`echo "dnase_$i"`
	echo "/data/tools/R-3.1.1/bin/Rscript all_tf_comparison_bp_all_methods_fix.R $i" | \
	qsub -l h_vmem=${MEM} -v PATH -cwd -N ${PROCESS_NAME} \
        -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"
done

