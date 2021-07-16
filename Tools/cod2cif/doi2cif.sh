#! /bin/bash

export DOI=$1

echo select file from data where doi like \"${DOI}\" > cod_query.tmp

mysql -u cod_reader -h sql.crystallography.net < cod_query.tmp cod > cif_ids.dat
rm cod_query.tmp

for (( i=2 ; i<=$( wc -l cif_ids.dat | awk '{print $1}' ) ; i++ ))
do
 export cif_id=$( sed -n ${i}p cif_ids.dat)
 curl -s http://www.crystallography.net/cod/${cif_id}.cif > ${cif_id}.cif
done
