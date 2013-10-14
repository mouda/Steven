
    ### Job name
    #PBS -N ULSA4b2_DC_runHN_N150_23
    ### out files
    #PBS -e ./log/ULSA4b2_DC_runHN_N150_23_50.err
    #PBS -o ./log/ULSA4b2_DC_runHN_N150_23_50.log
    ### put the job to which queue (qwork)
    #PBS -q qwork
 
    echo Working directory is $PBS_O_WORKDIR
    cd $PBS_O_WORKDIR
    echo running on host `hostname`
    echo Start time is `date`
    time1=`date +%s`
    echo Directory is `pwd`
 
    ../../../ULSA4b2_DC/ULSA4b2_DC 150 23 8 0 180 ../../mapFile/mapFile_uni_150_r500/mapFile_uniR500_N150_50.txt 0.477 2 1 1 0.95 0 20000 
    

    echo End time is `date`
    time2=`date +%s`
    echo Computing time is `echo $time2-$time1 | bc` sec
    
