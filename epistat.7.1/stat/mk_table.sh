#!/bin/bash
#This script join one column files into a tab delimited table
#Arguments <dir_to_search> <.files'_extension>

function join_by { local IFS="$1"; shift; echo "$*"; }

SFOLDER=$1
FEXT=$2

id_list=()
file_list=()
while IFS=  read -r -d $'\0'; do
	file_list+=("$REPLY")
	fname=`basename "$REPLY"`
	fname="${fname%%.*}"
	fname="${fname%%_fake*}"
    id_list+=("${fname}")
done < <(find ${SFOLDER} -name "*${FEXT}" -print0)
id_list[0]="#SAMPLE_ID:${id_list[0]}"
join_by $'\t' "${id_list[@]}"
printf '%s\n' "${file_list[@]}" | xargs paste