#!/bin/bash
#This script will perform some basic error checks after a PENTrack simulation has run.
#This script should be executed from the main PENTrack directory (ie one directory above the /out and /in folders)
#Sanmeet Chahal
#Nov 2, 2015

####Input parameters#####
expectedNumJobs=$1
outDirName=$2

if [[ $outDirName == "" ]]; then
	$outDirName=out
fi

if [[ $expectedNumJobs != "" && -d $outDirName ]]; then 
	printf "**************************************************************************************************************************************************\n" 
	printf "                                                               Post Run Errror Check Log                                                          \n" 
	printf "**************************************************************************************************************************************************\n"

	cd $outDirName

	###Move any jobs that were in the error folder back to the out folder and reprocess them##
	if [ -d ErrorJobs ]; then
		if [  "$(ls -A ErrorJobs)" ]; then 
			mv ./ErrorJobs/* ./
		fi
	
		rmdir ErrorJobs
	fi

	##### Count the number of the different log files produced ##########
	logs=( neutronend neutronsnapshot neutrontrack neutronhit neutronspin )
	logsForOutput=( end snapshot track hit spin )
	logCount=( 0 0 0 0 0 )

	for i in `seq 0 $((${#logs[@]}-1))`; do #to loop from 0 to size of log array -1
		temp=( *${logs[$i]}.out ) 
			
		if [ -e $temp ]; then #check if any of the logs of a particular type exist
			logCount[$i]=`ls *${logs[$i]}.out | wc -l`
		fi
	
		printf "Number of ${logsForOutput[$i]} logs produced: ${logCount[$i]}.\n" 
	done

	#####Check for non-empty error files#######
	a=`du -h error.txt* | awk '{ if ($1 == "0") print }' | wc -l`
	printf "\nNumber of empty error files produced: $a.\n" 

	#####Check for the number of non-empty error files#######
	b=`du -h error.txt* | awk '{ if ($1 != "0") print }' | wc -l` 
	printf "Number of non-empty error files produced: $b.\n\n" 

	errJobs=""
	notFinishJobs=""
	errExitStatusJobs=""

	####List the job number of the error files that are not empty####
	if [ $b != "0" ]; then 
		printf "The job number(s) of the non-empty error files are:\n"  
		errJobs=`du -h error.txt* | awk '{ if ($1 != "0") print }' | awk '{ print $2 }' | sed 's/error.txt-//g'`
		echo $errJob
	fi
	
	printf "\n"
	
	#######Count the number of output files that had an exit code, ie that reached completion######
	exitCodesCount=`cat output.txt* | grep -c "exit code"` 
	printf "Number of simulations that finished (had an exit code): $exitCodesCount.\n" 
	numZeroCodes=`cat output.txt* | grep -c "exit code 0"`
	printf "Number of simulations that finished normally (had an exit code of 0): $numZeroCodes.\n\n"

	#####Check for and list the output files that did not have an exit code (i.e. did not finish########
	if [ $exitCodesCount -lt $expectedNumJobs ]; then 
		printf "The jobs that did not finish (i.e. did not have an exit code):\n" 
	
		for outFile in $(ls output.txt-*); do
			exitCodeCount=`cat $outFile | grep -c "exit code"`
			jobNum=`echo $outFile | sed 's/output.txt-//g'`
		
			if [ $exitCodeCount -lt 1 ]; then 
				notFinishJobs=`echo "$notFinishJobs $jobNum"` #get a list of jobs with non-zero exit code
			fi
		done

		printf "$notFinishJobs\n" | sed -e 's/^[[:space:]]*//' #remove leading spaces
	fi

	#####Check for and list the output files that have non-zero exit codes (exit code is assigned by qsub job system)########
	if [ $numZeroCodes -lt $exitCodesCount ]; then 
		printf "\nThe jobs with non-zero qsub exit codes are:\n" 
	
		for outFile in $(ls output.txt*); do
			exitCode=`cat $outFile | grep "exit code" | awk '{print $7}'`
			jobNum=`echo $outFile | sed 's/output.txt-//g'`
			didJobFinsh=`echo $notFinishJobs | grep $jobNum`
		
			if [[ $exitCode != "0" && $didJobFinsh == "" ]]; then 
				errExitStatusJobs="$errExitStatusJobs $jobNum" #get a list of jobs with non-zero exit code
				printf "Job $jobNum has exit code $exitCode.\n" 
			fi
		done
	fi

	###Move all the jobs with non-zero exit code and non-empty error file to ErrorJobs folder####
	mkdir ErrorJobs
	errTot=`printf '%s\n' $errJobs $errExitStatusJobs $notFinishJobs | sort -nu` #get a unique list of error jobs

	if [[ $errTot != "" ]]; then
		printf "\nMoving all jobs with potential errors to $PWD/ErrorJobs directory:\n"

		for f in $errTot; do
			mv -t ./ErrorJobs output.txt-$f error.txt-$f 
		
			#list all output files -> remove all non numbers -> only print if the string has all zeros at the beginning and then the job number 
			#-> take the first one in case there are multiple output logs 
			properFileNamePrefix=`ls *.out | sed 's/[^0-9]//g' | grep '^0*'$f'$' | head -1`
		
			#check if each of the snapshot logs exist and move them to the ErrorJobs folder if they do
			if [[ -e ${properFileNamePrefix}neutronsnapshot.out ]]; then 
				mv ${properFileNamePrefix}neutronsnapshot.out ./ErrorJobs
			fi
	
			if [ -e ${properFileNamePrefix}neutronend.out ]; then 
				mv ${properFileNamePrefix}neutronend.out ./ErrorJobs
			fi
	
			if [ -e ${properFileNamePrefix}neutrontrack.out ]; then 
				mv ${properFileNamePrefix}neutrontrack.out ./ErrorJobs
			fi
	
			if [ -e ${properFileNamePrefix}neutronhit.out ]; then 
				mv ${properFileNamePrefix}neutronhit.out ./ErrorJobs
			fi
	
			if [ -e ${properFileNamePrefix}neutronspin.out ]; then 
				mv ${properFileNamePrefix}neutronspin.out ./ErrorJobs
			fi
	
			printf "Moved job $f output to $PWD/ErrorJobs folder.\n"
		done
	fi

	if [ -d out/ErrorJobs ] && [ "$(ls -A ErrorJobs)" ]; then
		rmdir ErrorJobs #delete the error jobs directory if empty
	fi
	
	##Find the average simulation time of all the neutrons as well as the min and max times####
	simTimes=`cat output.txt* | grep Simulation: | awk '{ print $4 }' | sed 's/[^0-9,.]*//g'` #store in an array
	simTimesArr=( $simTimes )
	minTime=`echo "${simTimes[@]}" | sort -n | head -1` 
	maxTime=`echo "${simTimes[@]}" | sort -n | tail -1`
	minTimeJob=`grep -H "Simulation: $minTime*s" output.txt-* | sed 's/:Init.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxTimeJob=`grep -H "Simulation: $maxTime*s" output.txt-* | sed 's/:Init.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totSimTime=0

	#echo simTimes: $simTimes
	#echo minTim: $minTime, max time: $maxTime, min time job: $minTimeJob, max time job: $maxTimeJob 

	##Find the total number of neutrons for each stopID####
	neutAbsSurf=`cat output.txt* | grep "neutron(s) were absorbed on a surface" | awk '{ print $2 }'` #stopID 2
	neutAbsSurfMat=( $neutAbsSurf ) 
	minNeutAbsSurf=`echo "${neutAbsSurf[@]}"  | sort -n | head -1`
	maxNeutAbsSurf=`echo "${neutAbsSurf[@]}" | sort -n | tail -1`
	minNeutAbsSurfJob=`grep -H "\s$minNeutAbsSurf    neutron(s) were absorbed on a surface" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutAbsSurfJob=`grep -H "\s$maxNeutAbsSurf    neutron(s) were absorbed on a surface" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutAbsSurf=0

	#echo neutAbsSurf: $neutAbsSurf
	#echo min: $minNeutAbsSurf, max: $maxNeutAbsSurf, min job: $minNeutAbsSurfJob, max job: $maxNeutAbsSurfJob

	neutAbsMat=`cat output.txt* | grep "neutron(s) were absorbed in a material" | awk '{ print $2 }'` #stopID 1
	neutAbsMatMat=( $neutAbsMat ) #name stands for: NEUTron ABSorbed in MATerial MATrix
	minNeutAbsMat=`echo "${neutAbsMat[@]}" | sort -n | head -1`
	maxNeutAbsMat=`echo "${neutAbsMat[@]}" | sort -n | tail -1`
	minNeutAbsMatJob=`grep -H "\s$minNeutAbsMat    neutron(s) were absorbed in a material" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutAbsMatJob=`grep -H "\s$maxNeutAbsMat    neutron(s) were absorbed in a material" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutAbsMat=0

	neutNotCat=`cat output.txt* | grep "neutron(s) were not categorized" | awk '{ print $2 }'` #stopID 0
	neutNotCatMat=( $neutNotCat ) 
	minNeutNotCat=`echo "${neutNotCat[@]}" | sort -n | head -1`
	maxNeutNotCat=`echo "${neutNotCat[@]}" | sort -n | tail -1`
	minNeutNotCatJob=`grep -H "\s$minNeutNotCat    neutron(s) were not categorized" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutNotCatJob=`grep -H "\s$maxNeutNotCat    neutron(s) were not categorized" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutNotCat=0

	neutNotFinish=`cat output.txt* | grep "neutron(s) did not finish" | awk '{ print $2 }'` #stopID -1
	neutNotFinishMat=( $neutNotFinish )
	minNeutNotFinish=`echo "${neutNotFinish[@]}" | sort -n | head -1`
	maxNeutNotFinish=`echo "${neutNotFinish[@]}" | sort -n | tail -1`
	minNeutNotFinishJob=`grep -H "\s$minNeutNotFinish    neutron(s) did not finish" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutNotFinishJob=`grep -H "\s$maxNeutNotFinish    neutron(s) did not finish" output.txt-* | sed 's/:.*//g'  | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutNotFinish=0

	neutHitBound=`cat output.txt* | grep "neutron(s) hit outer boundaries" | awk '{ print $2 }'` #stopID -2
	neutHitBoundMat=( $neutHitBound ) 
	minNeutHitBound=`echo "${neutHitBound[@]}" | sort -n | head -1`
	maxNeutHitBound=`echo "${neutHitBound[@]}" | sort -n | tail -1`
	minNeutHitBoundJob=`grep -H "\s$minNeutHitBound    neutron(s) hit outer boundaries" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutHitBoundJob=`grep -H "\s$maxNeutHitBound    neutron(s) hit outer boundaries" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutHitBound=0

	neutProdIntErr=`cat output.txt* | grep "neutron(s) produced integration error" | awk '{ print $2 }'` #stopID -3
	neutProdIntErrMat=( $neutProdIntErr )
	minNeutProdIntErr=`echo "${neutProdIntErr[@]}" | sort -n | head -1`
	maxNeutProdIntErr=`echo "${neutProdIntErr[@]}" | sort -n | tail -1`
	minNeutProdIntErrJob=`grep -H "\s$minNeutProdIntErr    neutron(s) produced integration error" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutProdIntErrJob=`grep -H "\s$maxNeutProdIntErr    neutron(s) produced integration error" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutProdIntErr=0

	neutDecayed=`cat output.txt* | grep "neutron(s) decayed" | awk '{ print $2 }'` #stopID -4
	neutDecayedMat=( $neutDecayed ) 
	minNeutDecayed=`echo "${neutDecayed[@]}" | sort -n | head -1`
	maxNeutDecayed=`echo "${neutDecayed[@]}" | sort -n | tail -1`
	minNeutDecayedJob=`grep -H "\s$minNeutDecayed    neutron(s) decayed" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutDecayedJob=`grep -H "\s$maxNeutDecayed    neutron(s) decayed" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutDecayed=0

	neutNoInitPos=`cat output.txt* | grep "neutron(s) found no initial position" | awk '{ print $2 }'` #stopID -5
	neutNoInitPosMat=( $neutNoInitPos ) 
	minNeutNoInitPos=`echo "${neutNoInitPos[@]}" | sort -n | head -1`
	maxNeutNoInitPos=`echo "${neutNoInitPos[@]}" | sort -n | tail -1`
	minNeutNoInitPosJob=`grep -H "\s$minNeutNoInitPos    neutron(s) found no initial position" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutNoInitPosJob=`grep -H "\s$maxNeutNoInitPos    neutron(s) found no initial position" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutNoInitPos=0
	
	neutEncountCGALErr=`cat output.txt* | grep "neutron(s) encountered CGAL error" | awk '{ print $2 }'` #stopID -6
	neutEncountCGALErrMat=( $neutEncountCGALErr ) 
	minNeutEncountCGALErr=`echo "${neutEncountCGALErr[@]}" | sort -n | head -1`
	maxNeutEncountCGALErr=`echo "${neutEncountCGALErr[@]}" | sort -n | tail -1`
	minNeutEncountCGALErrJob=`grep -H "\s$minNeutEncountCGALErr    neutron(s) encountered CGAL error" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutEncountCGALErrJob=`grep -H "\s$maxNeutEncountCGALErr    neutron(s) encountered CGAL error" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutEncountCGALErr=0

	neutEncountGeomErr=`cat output.txt* | grep "neutron(s) encountered geometry error" | awk '{ print $2 }'` #stopID -7
	neutEncountGeomErrMat=( $neutEncountGeomErr ) 
	minNeutEncountGeomErr=`echo "${neutEncountGeomErr[@]}" | sort -n | head -1`
	maxNeutEncountGeomErr=`echo "${neutEncountGeomErr[@]}" | sort -n | tail -1`
	minNeutEncountGeomErrJob=`grep -H "\s$minNeutEncountGeomErr    neutron(s) encountered geometry error" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	maxNeutEncountGeomErrJob=`grep -H "\s$maxNeutEncountGeomErr    neutron(s) encountered geometry error" output.txt-* | sed 's/:.*//g' | sed 's/output.txt-//g' | sort -n | head -1`
	totNeutEncountGeomErr=0

	###Next steps: add each of them and find the min and max ####
	numTot=`echo ${#simTimesArr[@]}`; #to find the average you should divide by the same number to get the average if every single output that finished had

	#sum up the parameters to be printed
	for i in `seq 0 $((numTot-1))`; do
		totSimTime=`echo "${simTimesArr[$i]}+$totSimTime" | bc`
	
		totNeutAbsSurf=`echo "${neutAbsSurfMat[$i]}+$totNeutAbsSurf" | bc` #stop id = 2
		totNeutAbsMat=`echo "${neutAbsMatMat[$i]}+$totNeutAbsMat" | bc` #stop id = 1
		totNeutNotCat=`echo "${neutNotCatMat[$i]}+$totNeutNotCat" | bc` #stop id = 0
		totNeutNotFinish=`echo "${neutNotFinishMat[$i]}+$totNeutNotFinish" | bc` #stop id = -1
		totNeutHitBound=`echo "${neutHitBoundMat[$i]}+$totNeutHitBound" | bc` #stop id = -2
		totNeutProdIntErr=`echo "${neutProdIntErrMat[$i]}+$totNeutProdIntErr" | bc` #stop id = -3
		totNeutDecayed=`echo "${neutDecayedMat[$i]}+$totNeutDecayed" | bc` #stop id = -4
		totNeutNoInitPos=`echo "${neutNoInitPosMat[$i]}+$totNeutNoInitPos" | bc` #stop id = -5
		totNeutEncountCGALErr=`echo "${neutEncountCGALErrMat[$i]}+$totNeutEncountCGALErr" | bc` #stop id = -6
		totNeutEncountGeomErr=`echo "${neutEncountGeomErrMat[$i]}+$totNeutEncountGeomErr" | bc` #stop id = -7
	done

	avTime=`echo "scale=2; $totSimTime/$numTot" | bc -l`
	#stopID averages
	avAbsSurf=`echo "scale=2; $totNeutAbsSurf/$numTot" | bc -l`
	avAbsMat=`echo "scale=2; $totNeutAbsMat/$numTot" | bc -l`
	avNotCat=`echo "scale=2; $totNeutNotCat/$numTot" | bc -l`
	avNotFinish=`echo "scale=2; $totNeutNotFinish/$numTot" | bc -l`
	avHitBound=`echo "scale=2; $totNeutHitBound/$numTot" | bc -l`
	avProdIntErr=`echo "scale=2; $totNeutProdIntErr/$numTot" | bc -l`
	avDecayed=`echo "scale=2; $totNeutDecayed/$numTot" | bc -l`
	avNoInitiPos=`echo "scale=2; $totNeutNoInitPos/$numTot" | bc -l`
	avEncountCGALErr=`echo "scale=2; $totNeutEncountCGALErr/$numTot" | bc -l`
	avEncountGeomErr=`echo "scale=2; $totNeutEncountGeomErr/$numTot" | bc -l`

	printf "\nStatistics for the Simulation:\n\n"
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Property" "Average" "Minimum" "Maximum" "Job Number for Minimum" "Job number for Maximum" 

	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Simulation time" "$avTime" "$minTime" "$maxTime" "$minTimeJob" "$maxTimeJob" 
	#print stats on stopIds
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons absorbed on a surface"         "$avAbsSurf"          "$minNeutAbsSurf"      "$maxNeutAbsSurf"      "$minNeutAbsSurfJob"        "$maxNeutAbsSurfJob" 
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons absorbed in a material"        "$avAbsMat"           "$minNeutAbsMat"       "$maxNeutAbsMat"       "$minNeutAbsMatJob"         "$maxNeutAbsMatJob" 
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons not categorized"               "$avNotCat"           "$minNeutNotCat"       "$maxNeutNotCat"       "$minNeutNotCatJob"         "$maxNeutNotCatJob" 
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons did not finish"                "$avNotFinish"        "$minNeutNotFinish"    "$maxNeutNotFinish"    "$minNeutNotFinishJob"      "$maxNeutNotFinishJob" 
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons hit outer boundary"            "$avHitBound"         "$minNeutHitBound"     "$maxNeutHitBound"     "$minNeutHitBoundJob"       "$maxNeutHitBoundJob" 
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons produced integration error"    "$avProdIntErr"       "$minNeutProdIntErr"   "$maxNeutProdIntErr"   "$minNeutProdIntErrJob"     "$maxNeutProdIntErrJob"
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons  decayed"                      "$avDecayed"          "$minNeutDecayed"      "$maxNeutDecayed"      "$minNeutDecayedJob"        "$maxNeutDecayedJob" 
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons found no initial position"     "$avNoInitiPos"       "$minNeutNoInitPos"    "$maxNeutNoInitPos"    "$minNeutNoInitPosJob"      "$maxNeutNoInitPosJob" 
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons encountered CGAL error"        "$avEncountCGALErr"   "$minNeutEncountCGALErr"    "$maxNeutEncountCGALErr"    "$minNeutEncountCGALErrJob"      "$maxNeutEncountCGALErrJob" 
	printf "%-45s | %12s | %12s | %12s | %25s | %25s\n" "Neutrons encountered geometry error"    "$avEncountGeomErr"   "$minNeutEncountGeomErr"    "$maxNeutEncountGeomErr"    "$minNeutEncountGeomErrJob"      "$maxNeutEncountGeomErrJob" 
else 
	printf "The $0 script requires that you enter an expected number of jobs as an argument and specify PENTrack output directory or have an \"out\" directory.\n"

fi



