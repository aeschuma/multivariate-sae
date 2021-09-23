#!/bin/sh
#
#

# request bash shell as shell for job
#$ -S /bin/bash
#
# set working directory
cd ~/Dropbox/dissertation_2/survey-csmf/out/stage-2-simulation

# set local variables

# execute command(s)
Rscript-R-3.6.1 --no-save --no-restore --verbose ~/Desktop/survey-csmf/stage-2-simulation/stage-2-sim-parallel.R $a $b $c $d $e $f $g $rr $s > stage-2-sim-parallel_$t.Rout 2>&1

# cp stage-2-sim-parallel_$a_$b_$c_$d_$e_$f_$g_$rr_$s.Rout ~/Dropbox/dissertation_2/survey-csmf/out/stage-2-simulation
# rm -f stage-2-sim-parallel_$a_$b_$c_$d_$e_$f_$g_$rr_$s.Rout

rm -f /home/users/aeschuma/s2sim_$t*
rm -f /home/students/aeschuma/qsub-stage-2-sims.sh*

# notify when script has finished
echo "Outfile in /home/users/aeschuma/Dropbox/dissertation_2/survey-csmf/stage-2-simulation/out" | mail -v -s "qsub-stage-2-sims-TEST.sh $a $b $c $d $e $f $g $rr $s finished" aeschuma@uw.edu
