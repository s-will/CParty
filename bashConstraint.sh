# input=$1
# for ((j = 1 ; j < 2 ; j=$((j+2)))); do

file="/home/mgray7/output2/constraintT/structures/rnafold.txt";

exec 5<$file
seq=""
i=0
while read line1 <&5; do

if ((i % 2 == 1))  
then

/usr/bin/time -o out.txt -f "%e\t%M" ./build/CParty -d2 -r $line1 $seq

cat "/home/mgray7/CParty/out.txt" >> "/home/mgray7/output2/constraintT/time/rnafold.txt"

else

seq=$line1
# echo $seq

fi 
i=$((i+1));
done

echo "first is $i";
exec 5</dev/null

# done



# for ((j = 1 ; j < 2 ; j=$((j+2)))); do

# file="/home/mgray7/output2/constraintT/structures/$j.txt";
# echo $j

# exec 5<$file
# seq=""
# i=0
# while read line1 <&5; do

# if ((i % 2 == 1))  
# then

# /usr/bin/time -o out.txt -f "%e\t%M" ./build/CParty -d2 -r $line1 $seq

# cat "/home/mgray7/CParty/out.txt" >> "/home/mgray7/output2/constraintT/time/$j.txt"

# else

# seq=$line1

# fi 
# i=$((i+1));
# done

# echo "first is $i";
# exec 5</dev/null

# done

# input=$1
# file="/home/mgray7/output2/fasta/$input.txt";
# file2="/home/mgray7/output2/structures/$input.txt";
# exec 5<$file
# exec 6<$file2
# while read line1 <&5 && read line2 <&6; do

# if ((i % 2 == 1))  

# then
# echo $line1 > out1.in
# echo $line2 >> out1.in


# ./build/CParty -d2 -n 20 $line1 > "/home/mgray7/output2/constraintTest/$i.txt"


# fi
# i=$((i+1));
# done
# echo "first is $i";
# #   echo "${line}";
# exec 5</dev/null
# exec 6</dev/null
