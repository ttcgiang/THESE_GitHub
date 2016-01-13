#!/bin/bash
for file in inputTemp/*
do
	echo $file
	newfile="outputTemp/"${file%%.*}"_out."${file#*.}
	echo newfile
	sh epis.sh $file $newfile
	sys
done
