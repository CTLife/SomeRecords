filelist=`ls  allwigs/*.wig`
i=1
for file in $filelist
do 
    echo $i
    let  i=i+1
    echo $file
done
