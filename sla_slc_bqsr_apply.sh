#bam=""
curDir=$(pwd)
cd $(dirname $0)

#source ../../auto_test_task/fun.sh

exe="/data1/new_cmd/"
path=/data2/output
tmpPath=/data2/fastq_tmp
#baseTmp=$tmpPath
baseTmp=""
base=$path
if [[ ! -d md5 ]];then
	mkdir -p md5
fi
md5db=md5/md5.db
file=../fastq_test.txt
#file=all.txt
#file=fastq.txt
#file=fastq2.txt
host=`hostname`
#checkFile=${file}.ok
inputFileFlag=0

if [[ $1 == "async" ]] && [[ ! -z $2 ]];then
inputFileFlag=1
file=$2
fi
if [[ ! -d $path ]];then
	mkdir -p $path
fi
if [[ ! -d log/std ]];then
	mkdir -p log/std
fi
runNum=1
bamSplit=/data1/bam.split
fastaPath=/data/ref/Bomo_genome_assembly/

fasta=/data1/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta
ksBin=$fastaPath/jiacan_practical_ks.binary
cmdLog=cmd.tmp.txt
trace=trace.log
echo "------------------" >> $cmdLog
echo "------------------" >> $trace
cnt=3
i=0

cat $file | while read line
do
	a=($line)
	r1=${a[0]}
	r2=${a[1]}
	bed=${a[2]}
	#tag=`echo $r1 | awk -F '/' '{print $2}'`
	fastqName=`basename $r1`
	output_awsPath=`basename $(dirname $r1)`
	awstmpPath=`dirname $r1`
	#awsPath=$(echo $awstmpPath | cut -f 4 -d /)"/RUN${runNum}/"$(echo $awstmpPath | cut -f 5-100 -d /)
	awsPath=$(echo $awstmpPath | cut -f 4 -d /)"/"$(echo $awstmpPath | cut -f 5-100 -d /)


	#prefix=${fastqName%%.*}
	prefix=${fastqName%%_R*}
	#prefix=$(echo $awsPath | awk -F "/" '{print $NF}')
	#bed=${fastqPath}/$(basename $fastqPath).bed
	#bed="/file1_data2/lhwes/sl_wes.bed"
	if [[ ! -z $bed ]];then
		bed="-B $bed"
	fi
	path=${base}/${output_awsPath}/$prefix
	#path=${base}/${awsPath}/$prefix
	if [[ ! -d $path ]];then
		mkdir -p $path
	fi

	if [[ ! -z $baseTmp ]];then

	inputPath=$baseTmp/${prefix}
	if [[ ! -d $inputPath ]];then
		mkdir -p $inputPath
	fi
	cmd="cp $r1 $r2 $inputPath"
	eval $cmd
		if [[ $? -ne 0 ]];then
			echo "ERR, cmd:$cmd" >> $trace
			continue;
		fi
	r1Name=`basename $r1`
	r2Name=`basename $r2`
	r1=$inputPath/$r1Name
	r2=$inputPath/$r2Name
	fi

	
		slaArgus=" --filter-enable --min-len 45 --cut-head --cut-tail --r1-adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --r2-adapter AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --max-n-base 10"
		#slbcArgus=""
		#slbcArgus="--with-umi"
		slbcArgus=""
		stdLog=log/std/${prefix}_std.log
		metrics=${path}/${prefix}.metrics
		tableName=${path}/${prefix}.table
cmd="${exe}/sla -R \"@RG\tID:TEST\tSM:${prefix}\tPL:ILLUMINA\" -B $bamSplit --filter-enable --min-len 45 --cut-head --cut-tail --r1-adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --r2-adapter AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --max-n-base 10 -i $metrics --bqsr-known-site /data1/Homo_sapiens_assembly38/know_sites2 /data1/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta --bqsr-tab $tableName $r1 $r2 >> $stdLog 2>&1"
		echo $cmd >> $cmdLog
		echo $cmd >> $stdLog
		bts=`date +%s`
		eval $cmd
		retCode=$?
		if [[ $retCode -ne 0 ]];then
			if [[ $retCode -eq 20 ]];then
				slaLoopFlag=1
				for slaLoop in {1..3}
				do
					eval $cmd
					if [[ $? -eq 0 ]];then
						slaLoopFlag=0
						break;
					else
						echo "ERR, cmd:$cmd" >> $trace
					fi
				done
				
				if [[ $slaLoopFlag -eq 1 ]];then
					continue
				fi
			else
				continue;
			fi
			
			
		fi
		cmd="rm -rf $r1 $r2"
		eval $cmd
		ets=`date +%s`
		slats=`expr $ets - $bts`

		bamName=${path}/${prefix}.bam
		mret=`md5sum $metrics`
		metricsMd5=`echo $mret | awk '{print $1}'`
		echo $mret >> $md5db
		echo "$metricsMd5  ${prefix}.metrics"> ${metrics}".md5"
		mret=`md5sum $tableName`
		tableMd5=`echo $mret | awk '{print $1}'`
		echo $mret >> $md5db
		echo "$tableMd5  ${prefix}.table" > ${tableName}".md5"
		#tableArr[$i]=$ret_md5



		gvcf=${path}/${prefix}.gvcf.gz
		cmd="${exe}/slc -e GVCF --keep-split --no-mark-supplementary$slbcArgus -P $bamSplit -b $bamName -R $fasta -o $gvcf -z $tableName $bed >> $stdLog 2>&1"
		echo $cmd >> $stdLog
		echo $cmd >> $cmdLog
		bts=`date +%s`
		eval $cmd
		retCode=$?
		if [[ $retCode -ne 0 ]];then
			echo "ERR, cmd:$cmd" >> $trace
			continue;
		fi
		ets=`date +%s`
		slcgvcfts=`expr $ets - $bts`

		mret=`md5sum $gvcf`
		gvcfMd5=`echo $mret | awk '{print $1}'`
		echo $mret >> $md5db
		echo "$gvcfMd5  ${prefix}.gvcf.gz" >> ${gvcf}".md5"

		vcf=${path}/${prefix}.vcf.gz
		cmd="${exe}/slgtvcf -i $gvcf -o $vcf >> $stdLog 2>&1" 
		echo $cmd >> $cmdLog
		echo $cmd >> $stdLog
		bts=`date +%s`
		eval $cmd
		if [[ $retCode -ne 0 ]];then
			echo "ERR, cmd:$cmd" >> $trace
			continue;
		fi
		ets=`date +%s`
		slcts=`expr $ets - $bts`

		mret=`md5sum $vcf`
		vcfMd5=`echo $mret | awk '{print $1}'`
		echo $mret >> $md5db
		echo "$vcfMd5  ${prefix}.vcf.gz" > ${vcf}".md5"

		mret=`md5sum $bamName`	
		bamMd5=`echo $mret | awk '{print $1}'`
		echo $mret >> $md5db
		echo "$bamMd5  ${prefix}.bam" > ${bamName}".md5"

		bamBqsrName=$path/${prefix}_bqsr.bam
		cmd="${exe}/slapplybqsr --compression-level 1 --tab $tableName --input $bamName --output $bamBqsrName >> $stdLog 2>&1"
		echo $cmd >> $stdLog
		echo $cmd >> $cmdLog
		bts=`date +%s`
		eval $cmd
		if [[ $? -ne 0 ]];then
			echo "ERR, cmd:$cmd" >> $trace
			continue;
		fi
		ets=`date +%s`
		bqsrbamts=`expr $ets - $bts`
		mret=`md5sum $bamBqsrName`
		bqsrbamMd5=`echo $mret | awk '{print $1}'`
		echo $mret >> $md5db
		echo "$bqsrbamMd5  ${prefix}_bqsr.bam" > ${bamBqsrName}".md5"
		total=`expr $slats + $slcts + $slcgvcfts + $bqsrbamts`
		echo "germline,$prefix,$slats,$slcts,$slcgvcfts,$bqsrbamts,$total,$bamMd5,$bqsrbamMd5,$vcfMd5,$gvcfMd5,$tableMd5,$metricsMd5" >> ret.csv
	#cmd="nohup bash upload.sh $path $awsPath $prefix &"
	cmd="bash upload.sh $path $awsPath $prefix"
	eval $cmd
	echo $cmd >> $cmdLog
	if [[ ! -z $baseTmp ]];then
		cmd="rm -rf $inputPath"
		eval $cmd
		echo $cmd >> $cmdLog
	fi

	#fi

	#i=`expr $i + 1`
	#left=`expr $i % $cnt`





#echo $prefix >> $checkFile
done
if [[ $inputFileFlag -eq 1 ]];then
	#date
	rm $(cat $file)
fi

#sleep 600
cd $curDir

