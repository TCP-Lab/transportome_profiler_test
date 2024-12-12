#!/bin/bash

# Run me by 'source' from the 'transportome_profiler' root directory
# and kill me by ctrl+C when you feel it's enough...
# Then you can use 'echo ${n_sets[@]}' to check the final array.

#set -e           # "exit-on-error" shell option
set -u           # "no-unset" shell option
set -o pipefail  # exit on within-pipe error

# Array to store the number of genesets made by Ariadne at each run of make_genesets.py
n_sets=(0 0)

while [[ ${n_sets[-2]} == 0 || ${n_sets[-1]} == ${n_sets[-2]} ]]; do

	# Do the cleaning
	if [[ -f ./data/genesets_repr.txt ]]; then
		printf "\nRemoving old genesets...\n"
		rm data/genesets*
	fi

	# Run Ariadne
	python ./src/modules/make_genesets.py ./data/MTPDB.sqlite \
		./data/in/config/gene_lists/basic.json \
		./data/genesets.json ./data/genesets_repr.txt \
		--min_pop_score 0.2 \
		--min_recurse_set_size 0 \
		--prune_direction "bottomup" \
		--prune_similarity 0.95
		#--no_prune

	n_report=$(wc -l ./data/genesets_repr.txt)
	echo $n_report
	
	n=$(echo $n_report | cut -f1 -d' ')
	n_sets+=($n)
done
