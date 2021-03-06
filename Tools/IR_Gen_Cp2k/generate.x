#! /bin/bash

export WORK_DIR=~/Nanoribbon/IR_3x  #### Directory where all the following files are and where new folders will be created
export data_file=GNR_Nit_3x_opt     #### It needs to end as .xyz
export step=0.010000                #### Displacement in Angstrom
export nats_unit=120                #### Number of atoms in the unit-cell. Atoms in data_file can be more if simulationg a super-cell
export INP_FILE=GNR_Nit_3x_IR       #### It needs to end as .inp

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

sed "s/XXX/${i}/g" job >  ${WORK_DIR}/${i}/job_${i}

echo 'atom '$i' done'

done



