file=./hap_config.txt 

cat $file | while read line
do 
	 a=($line) 
	 r1=${a[0]} 
	 r2=${a[1]} 
	 fastqName=`basename $r2` 
	 prefix=${fastqName%%.*} 
	 cur_path=$(pwd)
	 path="$cur_path/res/$prefix"
	 if [[ ! -d $path ]];then
		mkdir -p $path 
	 fi 
	docker run --rm -v /NAS/sl:/NAS/sl pkrusche/hap.py /opt/hap.py/bin/hap.py $cur_path/$r1 $cur_path/$r2 -r /NAS/sl/hg38/Homo_sapiens_assembly38.fasta -T $cur_path/HG002_GRCh38_GIAB.bed -o $path/$prefix --engine vcfeval 
done

