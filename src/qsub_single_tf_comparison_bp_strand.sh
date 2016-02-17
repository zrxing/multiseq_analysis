MEM=20g

#create directory for logs (if it doesn't already exist)
LOG_DIR="$HOME/projects/multiseq_analysis/results/roc/logs/"
if [ ! -d "${LOG_DIR}" ]; then
    mkdir ${LOG_DIR}
fi

#run multiseq on each locus in list_loci using qsub
PROCESS_NAME=`echo "dnase_roc"`
echo "/data/tools/R-3.1.1/bin/Rscript single_tf_comparison_bp_strand.R" | \
qsub -l h_vmem=${MEM} -v PATH -cwd -N ${PROCESS_NAME} \
        -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"

