#!/bin/bash -x
#SBATCH --partition=intr
#SBATCH --account=st48872
#SBATCH --time=1:00:00
####SBATCH --qos=small
#SBATCH --job-name=multiNode
#SBATCH -o SLURM.log
#SBATCH -e SLURM.err
#SBATCH --nodes=4
####SBATCH --ntasks-per-node=16

#--------------------------------------------------------------------------
# User variables
# @@@ WARNING @@@ Modify ~/.bashrc to allow environment variables in non interactive session
NPROCSHARED=8
targetExtension="*.xyz"
historyDir="history.optg"
command="M3C-adf.optg ../../PBE.optg-ADF.inp $NPROCSHARED"

# NPROCSHARED=4
# targetExtension="*.xyz"
# historyDir="history.optg"
# command="M3C-nwchem.optg ../b3lyp.optg-NWCHEM.inp $NPROCSHARED"
#--------------------------------------------------------------------------

nNodes=$SLURM_NNODES
workDir=$SLURM_SUBMIT_DIR
nodeName=( `scontrol show hostnames $SLURM_JOB_NODELIST` )

cd $workDir

nElem=`ls $targetExtension | wc -w`
nElemPerNode=`echo "$nElem $nNodes" | awk '{v=$1/$2; if(v==int(v)) print int(v); else print int(v+1)}'`

i=0; nodeId=-1
for f in `ls $targetExtension`
do
	if (( $i%$nElemPerNode==0 ))
	then
		nodeId=$(($nodeId+1))
		mkdir ${nodeName[$nodeId]}
	fi
	
	echo "$f --> ${nodeName[$nodeId]}"
	cp $f ${nodeName[$nodeId]}/
	i=$(($i+1))
done

for (( i=0; i<${#nodeName[@]}; i++ ))
do
	ssh -o StrictHostKeyChecking=no ${nodeName[$i]} "cd $workDir/${nodeName[$i]}; $command" &
done

wait

if [ -d "$historyDir" ]
then
        echo "@@@ WARNING @@@ There is already a history directory (history.freqs). Backup copy will be generated"
        mv $historyDir $historyDir.backup-`date +%Y%m%d`
fi
mkdir $historyDir

for (( i=0; i<${#nodeName[@]}; i++ ))
do
	mv ${nodeName[$i]}/$historyDir/* $historyDir/
	mv ${nodeName[$i]}/* .
	rm -rf ${nodeName[$i]}
done

