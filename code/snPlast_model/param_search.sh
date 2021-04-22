#!/bin/sh
for DA in $(seq 0.00075 0.0002 0.00275)
do 
	for ratio in $(seq 0.0 0.02 0.54)
	do
		ACh=$(bc <<< "scale=9;$DA*$ratio") 
		echo "running sim DA: $DA ACh:$ACh"
		python run_model.py -cr 30 -it 100 -tr 200 -sw 81 -eD $DA -eA $ACh -out "./model_results"
	done
done

