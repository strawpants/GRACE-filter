#!/bin/bash
#Roelof Rietbroek, 24 April 2020

#input approximate normal matrix
nggin=../Normal/ngg_bd_002-120.bin

#note that the program SH_createDDK is now part of RLFtlbx  
# factor are resp. 5, 20 and 50 times weaker than DDK8
#for fac in 1e9 2.5e8 1e8
#do
    #outf=Wbd_2-120.${fac}p_4
    #SH_createDDK -p${fac},4 -f$outf  $nggin
#done


#also create filters which are even smoother than ddk1 (by a factor 5 and 10 respectively)
for fac in 5e14 1e15 5e15
do
    outf=Wbd_2-120.${fac}p_4
    SH_createDDK -p${fac},4 -f$outf  $nggin
done
