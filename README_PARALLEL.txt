1-29-2015 Sasidhar Madugula Stanford University

The drug_analysis_edited_par.R script can be run in parallel on any 
Stanford Corn account with minimal adjustment. The only parameters 
which need to be changed are those in run_drug_analysis_par.sh. The 
user's pwd needs to be changed, and the directory they would like to 
put the output files in needs to be changed as well. The run_drug_analysis_par.sh
script needs to be in the same directory as the drug_analysis_edited_par.R script.
After the jobs have finished running (confirmed by typing "qstat"), the 
script process_qsub_output.sh should be run. Remember to match the
the output directory specified in the run_drug_analysis_par.sh 
script with that specified in the process_qsub_output.sh script.
