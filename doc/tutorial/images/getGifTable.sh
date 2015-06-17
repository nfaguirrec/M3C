#!/bin/bash

nColumns=$1
highlight=$2  # 1,4,5

echo "% ---------------------------------------------------------------------------------------------"
echo "% This table has been generated with the following command:"
echo "% $0 $*"
echo "% ---------------------------------------------------------------------------------------------"
echo ""

declare -A highlighted

for item in `echo $highlight | awk 'BEGIN{RS=","}{print $1}'`
do
	highlighted[$item]=1
done

files=`ls *.gif`
nFiles=`echo $files | wc -w`

echo "\begin{figure}[ht]"
echo "\centering"
echo "\begin{tabular}{|"
for i in `seq 1 $nColumns`
do
	echo ">{\centering\arraybackslash}p{1.6cm}|"
done
echo "}"
echo "\hline"

i=1
for f in $files
do
	if (( $i%$nColumns != 0 && i != $nFiles ))
	then
		if [ "${highlighted[$i]}" = 1 ]
		then
# 			echo "\includegraphics[bb = 14 14 115 115, scale=0.3]{images/table1/$f} \textcolor{red}{ \tiny{${f%.*}} } &"
			echo "\includegraphics[bb = 14 14 115 115, scale=0.3]{images/table1/$f} \textcolor{red}{\tiny{$i \hspace{5pt}${f%.*}}} &"
		else
# 			echo "\includegraphics[bb = 14 14 115 115, scale=0.3]{images/table1/$f} \tiny{${f%.*}} &"
			echo "\includegraphics[bb = 14 14 115 115, scale=0.3]{images/table1/$f} \tiny{$i \hspace{5pt} ${f%.*}} &"
		fi
	else
		if [ "${highlighted[$i]}" = 1 ]
		then
# 			echo "\includegraphics[bb = 14 14 115 115, scale=0.3]{images/table1/$f} \textcolor{red}{\tiny{${f%.*}} } "
			echo "\includegraphics[bb = 14 14 115 115, scale=0.3]{images/table1/$f} \textcolor{red}{\tiny{$i \hspace{5pt}${f%.*}}} "
		else
# 			echo "\includegraphics[bb = 14 14 115 115, scale=0.3]{images/table1/$f} \tiny{${f%.*}} "
			echo "\includegraphics[bb = 14 14 115 115, scale=0.3]{images/table1/$f} \tiny{$i \hspace{5pt} ${f%.*}} "
		fi
		
		if (( i != $nFiles ))
		then
			echo "\\\\\hline"
		else
			echo "\\\\\cline{1-$(($i%$nColumns))}"
		fi
	fi
	
	i=$(( $i+1 ))
done

# echo "\\\\\hline"
echo "\end{tabular}"
echo "\caption{\footnotesize{"
echo "}}"
echo "% \label{}"
echo "\end{figure}"
echo ""
echo "% ---------------------------------------------------------------------------------------------"