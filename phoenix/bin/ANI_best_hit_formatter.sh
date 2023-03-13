#!/bin/bash -l

#
# Description: Script to format ANI output to include more information on the top line
#
# Usage: ./ANI_best_hit_formatter.sh -e explicit_path_to_isolate_folder
#
# Output location: /sample_name/fastANI/
#
# Modules required: None
#
# V1.0 (06/03/2022)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#


#  Function to print out help blurb
show_help () {
	echo "./ANI_best_hit_formatter.sh -a ani_file -n sample_name"
	echo "required: -a = ani file"
}

# Parse command line options
options_found=0
while getopts ":h?a:n:t:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		a)
			echo "Option -a triggered, argument = ${OPTARG}"
			ani_file=${OPTARG};;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
		t)
			echo "Option -t triggered"
			terra=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

# set the correct path for bc/wget - needed for terra
if [[ $terra = "terra" ]]; then
	bc_path=/opt/conda/envs/phoenix/bin/bc
else
	bc_path=bc
fi

if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

if [[ ! -f ${ani_file} ]]; then
	echo "ani file does not exist, exiting"
	exit
fi

sorted_ani=${ani_file//.txt/.sorted.txt}

sort -k3 -n -r -o "${sorted_ani}" "${ani_file}"

best=$(head -n 1 "${sorted_ani}")
#Creates an array from the best hit
IFS='	' read -r -a def_array <<< "${best}"
best_file=${def_array[1]}
best_file=$(echo "${best_file}" | rev | cut -d'/' -f1 | rev)
best_percent=${def_array[2]}
fragment_matches=${def_array[3]}
total_fragments=${def_array[4]}
best_percent=$(echo "scale=2; $best_percent / 1" | $bc_path -l)
best_coverage=$(echo "scale=2; 100 * $fragment_matches / $total_fragments" | $bc_path -l)

# Pulling taxonomy from filename which was looked up. Can possibly be out of date. REFSEQ file will ALWAYS be current though
best_genus=$(echo "${best_file}" | cut -d'_' -f1)
best_species=$(echo "${best_file}" | cut -d'_' -f2)
best_organism_guess="${best_genus} ${best_species}"


#Creates a line at the top of the file to show the best match in an easily readable format that matches the style on the MMB_Seq log
echo -e "% ID	% Coverage	Organism	Source File" > "${sample_name}.fastANI.txt"
echo -e "${best_percent}	${best_coverage}	${best_organism_guess}	${best_file}" >> "${sample_name}.fastANI.txt"

### Add headers to file for Splunk integration
# sed 1i 'Isolate_Assembly_File	RefSEQ_Assembly_File	ANI_value	Mtaching_fragments	Total_fragments' "${sorted_ani}"


end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "ENDed ANI at ${end_time}"

#Script exited gracefully (unless something else inside failed)
exit 0
