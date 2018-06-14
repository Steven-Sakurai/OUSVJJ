
dirName="session"
i="1"
while [ $i -lt 16 ]
do 
    name=$dirName$i
    echo $name
    cd $name 
	nohup R CMD BATCH test.R &
	cd ..
    i=$[$i+1]
done
