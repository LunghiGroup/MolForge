#! /bin/bash

export WORK_DIR=/ichec/home/users/lunghia/VOdmit2/1x/IR  #### Directory where all the following files are and where new folders will be created
export data_file=VOdmit2_1x_opt     #### It needs to end as .xyz
export step=0.010000                #### Displacement in Angstrom
export nats_unit=432                #### Number of atoms in the unit-cell. Atoms in data_file can be more if simulationg a super-cell
export INP_FILE=VOdmit2_1x_IR       #### It needs to end as .inp
export numperjob=18		    #### number of atom distorions per job


mkdir -p ${WORK_DIR}

################### Generate job files in WORK_DIR/

export totjobs=$(( ${nats_unit}/${numperjob} ))

echo Total number of jobs: ${totjobs}

export i=1

for (( l=1 ; l<=${totjobs} ; l++ ))
do

cat > ${WORK_DIR}/job_${l} << endmsg
#!/bin/sh

#SBATCH --time=72:00:00   # 1 day and 3 hours
#SBATCH --nodes=1         # 16 cores
#SBATCH -A tcphy143b
#SBATCH -p ProdQ       # partition name
#SBATCH -J VO_${l}

module load cp2k/gfortran/7.1
source /ichec/packages/cp2k/gfortran/7.1/setup

ulimit -s unlimited
endmsg

for (( v=1 ; v<=${numperjob} ; v++ ))
do

 cat >> ${WORK_DIR}/job_${l} << endmsg2 
 export INP_DIR=${WORK_DIR}/${i}
 export OUT_DIR=${WORK_DIR}/${i}

 cd \${INP_DIR}

 for ll in {1..6}
 do

  export INP=${INP_FILE}_${i}_\${ll}.inp
  export OUT=${INP_FILE}_${i}_\${ll}.out

  mpirun -n 40 cp2k.popt -i \${INP_DIR}/\${INP} -o \${OUT_DIR}/\${OUT}
  
 done
endmsg2

 export i=$(( ${i} +1 ))
 done

done

################### Generate inpute_files and distortions in WORK_DIR/XXX


export nats=$(  head -n 1 ${data_file}.xyz )

echo 'Reading xyz coordinate file '$1
echo 'Looking for '$nats' atmos'

export offset=2

for (( i=1 ; i<=${nats_unit} ; i++ ))
do

export s=1

mkdir -p  ${WORK_DIR}/${i}

export after=$(sed -n $((${offset}+${i}))p ${data_file}.xyz | awk '{printf "%s %1.16f %1.16f %1.16f \n", $1,$2+'''${step}''',$3,$4}' )
export before=$(sed -n $((${offset}+${i}))p ${data_file}.xyz  )
sed "s/${before}/${after}/g" ${data_file}.xyz > ${WORK_DIR}/${i}/${data_file}_${i}_${s}.xyz

sed "s/XXX/${i}/g" ${INP_FILE}.inp > ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed "s/XXX/${i}/g" ${INP_FILE}.inp > ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g" ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g" ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
export s=$(($s+1))


export after=$(sed -n $((${offset}+${i}))p ${data_file}.xyz | awk '{printf "%s %1.16f %1.16f %1.16f \n", $1,$2-'''${step}''',$3,$4}' )
export before=$(sed -n $((${offset}+${i}))p ${data_file}.xyz  )
sed "s/${before}/${after}/g" ${data_file}.xyz >  ${WORK_DIR}/${i}/${data_file}_${i}_${s}.xyz


sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
export s=$(($s+1))


export after=$(sed -n $((${offset}+${i}))p ${data_file}.xyz | awk '{printf "%s %1.16f %1.16f %1.16f \n", $1,$2,$3+'''${step}''',$4}' )
export before=$(sed -n $((${offset}+${i}))p ${data_file}.xyz  )
sed "s/${before}/${after}/g" ${data_file}.xyz >  ${WORK_DIR}/${i}/${data_file}_${i}_${s}.xyz


sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
export s=$(($s+1))


export after=$(sed -n $((${offset}+${i}))p ${data_file}.xyz | awk '{printf "%s %1.16f %1.16f %1.16f \n", $1,$2,$3-'''${step}''',$4}' )
export before=$(sed -n $((${offset}+${i}))p ${data_file}.xyz  )
sed "s/${before}/${after}/g" ${data_file}.xyz >  ${WORK_DIR}/${i}/${data_file}_${i}_${s}.xyz


sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
export s=$(($s+1))

export after=$(sed -n $((${offset}+${i}))p ${data_file}.xyz | awk '{printf "%s %1.16f %1.16f %1.16f \n", $1,$2,$3,$4+'''${step}'''}' )
export before=$(sed -n $((${offset}+${i}))p ${data_file}.xyz  )
sed "s/${before}/${after}/g" ${data_file}.xyz >  ${WORK_DIR}/${i}/${data_file}_${i}_${s}.xyz


sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
export s=$(($s+1))

export after=$(sed -n $((${offset}+${i}))p ${data_file}.xyz | awk '{printf "%s %1.16f %1.16f %1.16f \n", $1,$2,$3,$4-'''${step}'''}' )
export before=$(sed -n $((${offset}+${i}))p ${data_file}.xyz  )
sed "s/${before}/${after}/g" ${data_file}.xyz >  ${WORK_DIR}/${i}/${data_file}_${i}_${s}.xyz


sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed "s/XXX/${i}/g" ${INP_FILE}.inp >  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
sed -i "s/YYY/${s}/g"  ${WORK_DIR}/${i}/${INP_FILE}_${i}_${s}.inp
export s=$(($s+1))

###sed "s/XXX/${i}/g" job >  ${WORK_DIR}/${i}/job_${i}

echo 'atom '$i' done'

done


