#!/bin/sh

for n in 50 200 400
do

for pathway in map00564 map02010
do

qsub -cwd -e iotrash/ -o iotrash/ ./$1  $n $pathway

done
done
