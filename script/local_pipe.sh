#! /bin/zsh

read -A dirs
for dir in $dirs
do
	cd $dir
	rm -r work
	julia -t 4 top.jl 15
	julia top.jl 7 occ
	julia top.jl 7 EE
	cd ..
done
./move.sh
python3 qedynachg.py 20 $dirs 
