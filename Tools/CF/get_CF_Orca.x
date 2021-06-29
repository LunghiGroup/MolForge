#! /bin/bash
 
 export CF_DIR=XXXSETDIRXXX
 export orca_output=$1
 export Nj=$2
 export lmax=$3

 export alpha=$4
 export beta=$5
 export gamma=$6


 export soc_size=$( grep 'Dim(SO)' ${orca_output} | awk '{print $3}' | tail -n 1 ) 

 export num_blocks=$(( ${soc_size}/6   ))

 if [ $(( ${soc_size} % 6   )) -ne 0 ]
 then  
  export num_blocks=$(( ${num_blocks}+1 ))
 fi

 export block_size=$(( ((${soc_size}+1) * ${num_blocks}) + 1 ))

 echo Reading ORCA output file ${orca_output}
 echo Projecting the lowest ${Nj} states of ${soc_size} CI solutions on a CF operator of order ${lmax}

 rm -f Os.dat
 
 grep -A $(( ${block_size} +2 )) 'SZ' ${orca_output} | tail -n ${block_size} > Sz.dat
 grep -A $(( ${block_size} +2 )) 'SX' ${orca_output} | tail -n ${block_size} > Sx.dat
 grep -A $(( ${block_size} +2 )) 'LZ' ${orca_output} | tail -n ${block_size} > Lz.dat
 grep -A $(( ${block_size} +2 )) 'LX' ${orca_output} | tail -n ${block_size} > Lx.dat
 grep -A $(( ${block_size} +2 )) 'SOC MATRIX (A.' ${orca_output} | tail -n ${block_size} > SOCR.dat
 grep -A $(( ${block_size}*2 +3 )) 'SOC MATRIX (A.' ${orca_output} | tail -n ${block_size} > SOCI.dat
 sed -i "s/-/ -/g" SOCR.dat
 sed -i "s/-/ -/g" SOCI.dat
 sed -i "s/-/ -/g" Sz.dat
 sed -i "s/-/ -/g" Sx.dat
 sed -i "s/-/ -/g" Lz.dat
 sed -i "s/-/ -/g" Lx.dat

 if [ ${alpha:-1000} -eq "1000" ] 
 then

  ${CF_DIR}/CF.x -JMult ${Nj} -lmax ${lmax} -CISIZE ${soc_size} 

 else

  ${CF_DIR}/CF.x -JMult ${Nj} -lmax ${lmax} -CISIZE ${soc_size} -rot ${alpha} ${beta} ${gamma}

 fi


