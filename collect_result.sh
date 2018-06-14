echo '' > est.csv
echo '' > conv.csv
echo '' > success.csv
dirName="session"
i="1"
while [ $i -lt 16 ]
do
    name=$dirName$i
    cat $name/est.csv >> est.csv
    cat $name/conv.csv >> conv.csv
    cat $name/success.csv >> success.csv
    i=$[$i+1]
done

R CMD BATCH clean_data.r 
