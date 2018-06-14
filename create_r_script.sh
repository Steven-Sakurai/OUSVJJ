i="1"
dirName="session"
while [ $i -lt 16 ]
do
	sentence="set.seed($RANDOM)"
	name=$dirName$i
	(echo $sentence && cat test.R) > $name/test.R
	i=$[$i+1]
done