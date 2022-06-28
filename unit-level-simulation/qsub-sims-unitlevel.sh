#!/bin/sh
#
#

# request bash shell as shell for job
#$ -S /bin/bash
#
# set working directory
cd ~/Dropbox/dissertation_2/survey-csmf/out/unit-level-simulation

# set local variables

# execute command(s)
Rscript --no-save --no-restore --verbose ~/Desktop/survey-csmf/unit-level-simulation/sim-unitlevel-parallel.R $a $rr $s > sim-unitlevel-parallel_$t.Rout 2>&1

# cp stage-2-sim-parallel_$a_$b_$c_$d_$e_$f_$g_$rr_$s.Rout ~/Dropbox/dissertation_2/survey-csmf/out/stage-2-simulation
# rm -f stage-2-sim-parallel_$a_$b_$c_$d_$e_$f_$g_$rr_$s.Rout

rm -f /home/users/aeschuma/ul_sim_$t*
rm -f /home/students/aeschuma/qsub-sims-unitlevel.sh*

# notify when script has finished
# echo "Outfile in /home/users/aeschuma/Dropbox/dissertation_2/survey-csmf/stage-2-simulation/out" | mail -v -s "qsub-stage-2-sims-TEST.sh $a $b $c $d $e $f $g $rr $s finished" aeschuma@uw.edu
