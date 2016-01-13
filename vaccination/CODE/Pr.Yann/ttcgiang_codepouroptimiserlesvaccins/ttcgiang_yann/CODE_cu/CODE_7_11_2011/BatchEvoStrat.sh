#!/bin/bash
for file in Input/*
do
	echo $file
	newfile="Output/"${file%%.*}"_out."${file#*.}
	echo $newfile
#	sh epis.sh $file $newfile
	./EvoStrat -EvoStratInput $file -EvoStratOutput $newfile

done
