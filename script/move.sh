#! /bin/zsh

for dir in $(ls -d sd*)
do
	cp ${dir}/work/occ $dir
	cp ${dir}/work/occ ${dir}/expN
	cp ${dir}/work/EE $dir
	cp ${dir}/work/bonds $dir
done
