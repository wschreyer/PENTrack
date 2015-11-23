#!/bin/bash
#This script will perform some basic preliminary checks of the input files of PENTrack before launching
#The script should be executed from PENTrack's main directory (ie one level higher than the /in and /out 
#directory's default location
#Sanmeet Chahal
#July 10, 2015

#### Possible user inputs ####
minJobNum=$1
maxJobNum=$2
jobName=$3
inDirName=$4
outDirName=$5
pbsfileName=$6
#wall time in hours or minutes (specify hours by h and minutes by m so 1h20m is interpreted as 1:20:00
#maximum is 72 hours, minimum is 15 minutes (usually the minimum time required to run PENTrack)
#wallTimeReq=$4 #saving this for future work

if [[ $inDirName == "" ]]; then 
	inDirName="in"
fi

if [[ $outDirName == "" ]]; then 
	outDirName="out"
fi

##Always regenerate the same template 
if [[ $pbsfileName == "" ]]; then 
	pbsfileName="batch.pbs"
fi

if [ -e $pbsfileName ]; then 
	rm $pbsfileName
fi

printf "#!/bin/bash\n\
#PBS -S /bin/bash\n\
#PBS -l mem=1gb  \n\
#PBS -l procs=1\n\
#PBS -l walltime=1:30:00 \n\
#PBS -t 1-10 \n\
#PBS -N SampleRun\n\
#PBS -e /home/username/error.txt  \n\
#PBS -o /home/username/output.txt \n\
#### For getting emails: PBS -M yourName@email.com #for receiving emails when certain tasks are complete (specified by PBS -m option): \
#PBS -m be #emails when the job is aborted(a),begins execution(b), is terminate(e)\n\n\
# Script testing out running multiple instances of a program\n\
cd \$PBS_O_WORKDIR\n\
echo \"Current working directory is \`pwd\`\"\n\
echo \$PBS_JOBNAME\n\
echo \"Starting run at: \`date\`\"\n\
./PENTrack \$PBS_ARRAYID \"./$inDirName\" \"./$outDirName\"\n\
#also put log files .o / .e in a subfolder in same output folder\n\
echo \"Program diffuse finished with exit code \$? at: \`date\`\"\n\n" > $pbsfileName

printf "******************************************************************************************************************************************************\n" 
printf "                                                                   Pre Run Errror Check Log                                                           \n" 
printf "******************************************************************************************************************************************************\n"

#################################################################
####Check all STL file referenced in geometry.in exist ##########
#################################################################

# print geometry.in -> sed: remove everything before the [GEOMETRY] string -> sed: remove everything after [SOURCE]
# -> grep: remove all lines that begin with a # (comment) -> grep: remove all lines with spaces -> don't include any 
#line that has 'ignored' in it 
if [ -e $inDirName/geometry.in ]; then 
	allSTL=`cat $inDirName/geometry.in | sed -e '1,/\[GEOMETRY\]/d' | sed -e '/\[SOURCE\]/,$d' | grep -v '^#' | grep -v '^ *$' | grep -v 'ignored' | awk '{print $2}'`

	#create a file that will hold the line number of all the files referenced
	cat -n $inDirName/geometry.in | sed -e '1,/\[GEOMETRY\]/d' | sed -e '/\[SOURCE\]/,$d' | awk '{print $1 " " $3}' > temp.txt

	#create array to loop over files
	allSTLArray=($allSTL) 

	#so that the message is only printed when necessary
	for filePath in "${allSTLArray[@]}"; do
		if [ ! -f $filePath ]; then #-f is for file exists
			printf "Paths specified in geometry.in for which no file exists:\n\n"
			printf "%-150s | %-10s \n" "File Path" "Line Number"
			break
		fi
	done

	#loop over the file paths in the geometry.in and check that they exist
	for filePath in "${allSTLArray[@]}"; do
		#check if file is missing
		if [ ! -f $filePath ]; then 
			lineNumNotExit=`cat temp.txt | grep $filePath | awk '{print $1}'`
			printf "%-150s | %+10s \n" $filePath $lineNumNotExit	 
		fi  
	done
fi

##############################################################################
####Check STL files in a *STL* folder are referenced in geometry.in ##########
##############################################################################

if [ -d $inDirName/*STL*/ ]; then 
	printf "\nFiles in the STL folder that are not referenced in geometry.in:\n"
	stlFolderPath=`ls -d $inDirName/*STL*/`
	allStlFromFolder=`ls $inDirName/*STL*/ | awk '{ print "'"$stlFolderPath"'" $0 }'`
	allStlFromFolderArray=( $allStlFromFolder )

	for stlFile in "${allStlFromFolderArray[@]}"; do
		c=`echo $allSTL | grep $stlFile`
		
		if [ -z "$c" ]; then
			printf "$stlFile\n"
		fi
	done
fi

##########################################################################
#######Set the starting and ending job numbers in the batch.pbs file #####
##########################################################################

if [[ $minJobNum != "" ]] && [[ $maxJobNum != "" ]];  then
	cat $pbsfileName | sed 's/#PBS -t.*/#PBS -t '$minJobNum'-'$maxJobNum' /g'  > toto.txt
	mv toto.txt $pbsfileName
	
	printf "\nFinished setting initial and final job numbers in $pbsfileName.\n"
fi

#############################################################
##### Set output file paths in batch.pbs file  ##############
#############################################################

outputFilesPath=`echo $PWD/$outDirName`
#print $pbsfileName -> find the line with #PBS -e (labelling it with a line number) -> delete everything after the line number
lineNumErrorPath=`cat $pbsfileName | grep -n "#PBS -e" | sed 's/[:].*//'`
lineNumOutputPath=`cat $pbsfileName | grep -n "#PBS -o" | sed 's/[:].*//'`

#replace the error path (need to use a | as the separator because  $outputFilesPath contains a / char
#need to use double quotes because we are replacing with something that is a variable
cat $pbsfileName | sed "s|#PBS -e.*|#PBS -e $outputFilesPath/error.txt |g" > toto.txt
mv toto.txt $pbsfileName

#replace the output path
cat $pbsfileName | sed "s|#PBS -o.*|#PBS -o $outputFilesPath/output.txt|g" > toto.txt
mv toto.txt $pbsfileName	

printf "Finished setting path of the error and output files in $pbsfileName.\n"

##########################################################################
#####Check and set the memory set for PENTrack in the batch.pbs file #####
##########################################################################

#always set the memory requested to 5000mb
#For the future: 
#convert the du -ch in/ | grep "tot" to a memory required
#multiply that number by 10 safety factor
cat $pbsfileName | sed 's/#PBS -l pmem=.*/#PBS -l pmem=1gb/g' > toto.txt
mv toto.txt $pbsfileName

printf "Finished setting the memory requested in $pbsfileName.\n"

####################################################
#######Set the job name if there is user input #####
####################################################

if [[ $jobName != "" ]];  then
	cat $pbsfileName | sed 's/#PBS -N.*/#PBS -N '$jobName' /g' > toto.txt
	mv toto.txt $pbsfileName
	
	printf "Finished setting the job name in $pbsfileName to $jobName.\n"
fi

########################################################################################################################
####### Check if number of jobs running + number of jobs specified in batch.pbs > 2880 (the max allowed by Jasper) #####
########################################################################################################################

jobNumberStart=`grep '^#PBS -t' $pbsfileName | awk '{print $3}' | sed "s/-.*$//"` 
jobNumberEnd=`grep '^#PBS -t' $pbsfileName | awk '{print $3}' | sed "s/^[^-]*-//"`
jobsRequested=$(($jobNumberEnd-$jobNumberStart)) 
jobsRunning=$(qstat -tu $USER | grep $USER | wc -l) 
totJobNum=$(($jobsRequested+jobsRunning))

if [[ $totJobNum -ge 2880 ]]; then 
	printf "\nWARNING!!! You will have $totJobNum jobs in queue if you submit this job. This is more than the 2880 jobs allowed on the Jasper cluster.\n\n"
else
	printf "\nSubmitting this job would exceed your job quota of 2880 jobs on the Jasper cluster.\n\n"
fi  

##################################################################
###Check if job numbers in batch.pbs and /out folder conflict ####
##################################################################
files=$(ls $outDirName/output.txt* 2> /dev/null | wc -l)

if [ "$files" != 0 ]; then
	outputNums=`ls $outDirName/output.txt* | sed 's/output.txt-//g' | sed "s|$outDirName/||g"`
fi

files=$(ls $outDirName/error.txt* 2> /dev/null | wc -l)

if [ "$files" != 0 ]; then
	errorNums=`ls $outDirName/error.txt* | sed 's/error.txt-//g' | sed "s|$outDirName/||g"`
fi

files=$(ls $outDirName/*.out | grep -v "BFCut" 2> /dev/null | wc -l)

if [  "$files" != 0 ]; then 
	outLogNums=`ls $outDirName/*.out | sed "s|$outDirName/||g" | sed 's/0//g' | sed 's/.out//g' | sed 's/proton//g' | sed 's/neutron//g' | sed 's/electron//g' | sed 's/out\///g' \
		| sed 's/end//g' | sed 's/snapshot//g' | sed 's/hit//g' | sed 's/track//g' | sed 's/spin//g'`
fi 

#check that jobStartNum is bigger than any *Nums and jobEndNum is too
allNums=`printf '%s\n' $outputNums $errorNums $outLogNums | sort -nu` #sort: -n for numerical, -u for unique
lastLineEnter=0

#so that it is only printed when necessary
for n in $allNums; do
	if [ $jobNumberStart -le $n ] && [ $jobNumberEnd -ge $n ]; then
		lastLineEnter=1
		printf "Cases of possible output overwrite! Job numbers with output in out/ and within the range specified by $pbsfileName:\n"
		break
	fi
done

#print out the job numbers that would overwrite the existing data
for n in $allNums; do 		
	if [ $jobNumberStart -le $n ] && [ $jobNumberEnd -ge $n ]; then
		printf "$n " 
	fi
done

if [ $lastLineEnter -eq 1 ]; then
	printf "\n"
fi

rm temp.txt



