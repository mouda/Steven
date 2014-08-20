Fidelity="0.6 0.65 0.7 0.75 0.8"
NumNode=50
NumHead=13
Correlation="0.477"

for item in $Fidelity ; do
  for (( i = 1; i <= 50; i++ )); do
    echo "=========== Fidelity : $item, Index : $i ========="
    echo "
    ### Job name
    #PBS -N genCluster_N50_${i}_F${item}
    ### out files
    #PBS -e ./log/genCluster_N50_F${item}_${i}.err
    #PBS -o ./log/genCluster_N50_F${item}_${i}.log
    ### put the job to which queue (qwork)
    #PBS -q qwork" > ./run_temp.sh
    echo ' 
    echo Working directory is $PBS_O_WORKDIR
    cd $PBS_O_WORKDIR
    echo running on host `hostname`
    echo Start time is `date`
    time1=`date +%s`
    echo Directory is `pwd`' >> ./run_temp.sh
    echo " 
    ../../../cluster/genCluster -n ${NumNode} -H ${NumHead} -q 32 -p 25 -b 180.0 -t 0.5 -I 10.0 -m  \
    ../../mapFile/mapFile_uni_50_r150/mapFile_uniR150_N50_$i.txt -f ${item} -c 0.477  -T  1.0 \
    -F MinPower -N 3  -C data/CS_N50_H13_${i}_F${item}.out
    " >> ./run_temp.sh
    echo '
    echo End time is `date`
    time2=`date +%s`
    echo Computing time is `echo $time2-$time1 | bc` sec
    ' >> ./run_temp.sh
    qsub run_temp.sh
  done
done


