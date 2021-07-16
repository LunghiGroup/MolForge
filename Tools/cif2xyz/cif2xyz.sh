#! /bin/bash

export MOLFORGE=~/ERC/MolForge

export cifname=$1

export name=$( echo $cifname | sed 's/.cif//g' )

############## Get Lattice Parameters from Cif file

export vec_a=$( grep 'cell_length_a' ${cifname} | awk '{printf "%1.4f \n",$2}' )
export vec_b=$( grep 'cell_length_b' ${cifname} | awk '{printf "%1.4f \n",$2}' )
export vec_c=$( grep 'cell_length_c' ${cifname} | awk '{printf "%1.4f \n",$2}' )

export vec_al=$( grep 'cell_angle_alpha' ${cifname} | awk '{printf "%1.4f \n",$2}' )
export vec_be=$( grep 'cell_angle_beta' ${cifname} | awk '{printf "%1.4f \n",$2}' )
export vec_ga=$( grep 'cell_angle_gamma' ${cifname} | awk '{printf "%1.4f \n",$2}' )

echo ${vec_a} ${vec_b} ${vec_c} ${vec_al} ${vec_be} ${vec_ga} > ${name}.abc

${MOLFORGE}/bin/abc2cell.x ${name}.abc ${name}.cell

####################################

rm -f ${name}.xyz
atomsk $cifname -remove-doubles 0.6 ${name}.xyz

${MOLFORGE}/bin/find_mols.x -xyz ${name}.xyz -cell ${name}.cell -reorder_mols -remap_mols > ${name}_remap.xyz





