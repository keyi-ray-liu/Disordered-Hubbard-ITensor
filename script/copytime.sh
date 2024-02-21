#! /bin/zsh

for file in $(ls -d sd*)
do
	cp timescale $file
done
