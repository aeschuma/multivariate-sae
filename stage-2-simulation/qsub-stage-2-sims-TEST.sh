#!/bin/sh
#
#

# request bash shell as shell for job
#$ -S /bin/bash
#
# set working directory
cd /scratch

# set local variables

# execute command(s)
Rscript --no-save --no-restore --verbose /home/students/aeschuma/Desktop/survey-csmf/stage-2-simulation/stage-2-sim-parallel-TEST.R $a $b $c $d $e > stage-2-sim-parallel-TEST_$a_$b_$c_$d_$e.Rout 2>&1

cp stage-2-sim-parallel-TEST_$a_$b_$c_$d_$e.Rout /home/students/aeschuma/Dropbox/dissertation_2/survey-csmf/out/stage-2-simulation
rm -f stage-2-sim-parallel-TEST_$a_$b_$c_$d_$e.Rout

rm -f /home/students/aeschuma/qsub-stage-2-sims-TEST.sh*
rm -f /home/students/aeschuma/s2sim*

# notify when script has finished
# echo "Outfile in /home/students/aeschuma/Dropbox/dissertation_2/survey-csmf/stage-2-simulation/out" | mail -v -s "qsub-stage-2-sims-TEST.sh $a $b $c $d $e $f $g $h $ii $j $k $l $m $n $o $p $q $r $s $t $u $v $w $x finished" aeschuma@uw.edu
