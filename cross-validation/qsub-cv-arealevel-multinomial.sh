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
Rscript --no-save --no-restore --verbose /home/students/aeschuma/Desktop/survey-csmf/cross-validation/ccv-fit-arealevel-multinomial.R

# notify when script has finished
# echo "Outfile in /home/students/aeschuma/Dropbox/dissertation_2/survey-csmf/stage-2-simulation/out" | mail -v -s "qsub-stage-2-sims-TEST.sh $a $b $c $d $e $f $g $h $ii $j $k $l $mm $n $o $p $q $r $rr $s $t $u $v $w $x finished" aeschuma@uw.edu
