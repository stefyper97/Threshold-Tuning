#!/bin/bash

#read the path

path_to_folders=.
if [[ -n $1 ]]; then
	path_to_folders=$3
fi
stave_id=$1
Vbb=$2

#echo ${path_to_folders}
#ls ${path_to_folders}
current_path=$(pwd)

outfile_name="stave_${stave_id}_Vbb_${Vbb}_filelist.txt"

for d in ${path_to_folders}/*
	do 
		echo "$stave_id $Vbb $current_path${d##.}" >> "${outfile_name}"
	done