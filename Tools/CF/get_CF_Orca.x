#! /bin/bash

 export orca_output=$1

 rm Os.dat

 grep -A 2670 'SZ' ${orca_output} | tail -n 2668 > Sz.dat
 grep -A 2670 'SX' ${orca_output} | tail -n 2668 > Sx.dat
 grep -A 2670 'LZ' ${orca_output} | tail -n 2668 > Lz.dat
 grep -A 2670 'LX' ${orca_output} | tail -n 2668 > Lx.dat
 grep -A 2670 'SOC MATRIX (A.' ${orca_output} | tail -n 2668 > SOCR.dat
 grep -A 5339 'SOC MATRIX (A.' ${orca_output} | tail -n 2668 > SOCI.dat
 sed -i "s/-/ -/g" SOCR.dat
 sed -i "s/-/ -/g" SOCI.dat
 sed -i "s/-/ -/g" Sz.dat
 sed -i "s/-/ -/g" Sx.dat
 sed -i "s/-/ -/g" Lz.dat
 sed -i "s/-/ -/g" Lx.dat

 ./CF.x

