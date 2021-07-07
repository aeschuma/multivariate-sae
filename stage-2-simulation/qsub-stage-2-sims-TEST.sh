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
Rscript --no-save --no-restore --verbose /home/students/aeschuma/Desktop/survey-csmf/stage-2-simulation/stage-2-sim-parallel-TEST.R $a $b $c $d $e $f $g $h $ii $j $k $l $m $n $o $p $q $r $s > stage-2-sim-parallel-TEST_$a_$b_$c_$d_$e_$f_$g_$h_$ii_$j_$k_$l_$m_$n_$o_$p_$q_$r_$s.Rout 2>&1

cp stage-2-sim-parallel-TEST_$a_$b_$c_$d_$e_$f_$g_$h_$ii_$j_$k_$l_$m_$n_$o_$p_$q_$r_$s.Rout /home/students/aeschuma/Dropbox/dissertation_2/survey-csmf/stage-2-simulation/out
rm -f stage-2-sim-parallel-TEST_$a_$b_$c_$d_$e_$f_$g_$h_$ii_$j_$k_$l_$m_$n_$o_$p_$q_$r_$s.Rout
rm -f /home/students/aeschuma/qsub-stage-2-sims-TEST.sh*
rm -f /home/students/aeschuma/s2sim_$a_$b_$c_$d_$e_$f_$g_$h_$ii_$j_$k_$l_$m_$n_$o_$p_$q_$r_$s*

# notify when script has finished
echo "Outfile in /home/students/aeschuma/Dropbox/dissertation_2/survey-csmf/stage-2-simulation/out" | mail -v -s "qsub-stage-2-sims-TEST.sh $a $b $c $d $e $f $g $h $ii $j $k $l $m $n $o $p $q $r $s $t $u $v $w $x finished" aeschuma@uw.edu
