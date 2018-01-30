#!/bin/bash
#
#This script will prepare the STL files in the required format for a PENTrack simulation. 
#It will also remove any spaces from the names of the STL files if they exist. 
# The script takes two optional arguments: 
#
#
# @param argv1  the name of the STL files dir as an argument. If no argument given, it assumes the dir has the name STL. 
#
# The STL files are assumed to have names with the following format:
# 			partName_mater_CuBe_prior_1.stl 
# and if it is a source volume then: 
# 			partName_mater_CuBe_prior_source.stl  
#where CuBe could be an arbitrary material type and 1 represents the priority of the material. 
#The script will then print out an ordered list of stl files in a format that can be copy
#pasted for the [GEOMETRY] section of the geometry.in input file of PENTrack. 
# 
# If no material type is specified the material type is listed as MISSING_MATERIAL and if no priority is specified it is listed 
# with the lowest priority. 
#The script requires you do not use at least one of the following characters in the STL file names: "}", "?", "<", "&", "!", "@"
#
#Sometimes you may have STL files that you want PENTrack to ignore, set the material to "NotUsed" (case insensitive) in the name of these files 
#and they will be ommitted from the output of the script.
#By: Sanmeet Chahal
#August 26, 2015

#the directory in which all the STL files are stored
stlDirName=$1

# If no name is provided it will assume that the STL dir has a name with the word STL in it 
if [ "$stlDirName" == "" ]; then
	stlDirName="STL"
fi

#check to make sure the stl dir exists 
if [ -d $stlDirName ]; then 

	stlDirName=`readlink -f $stlDirName` #make path absolute

	#remove the spaces from the STL file names
	for f in $stlDirName/*.STL;  do 
		a=`echo $f | sed 's/ //g'`
	
		#only necessary to rename the file if it contained a space
		#ie don't rename if the file didn't have a space
		if [ "$f" != "$a" ]; then 
			mv "$f" "$a"
		fi
	done

	#get a sorted list of STL files based on the priority 
	#explanation of command: 
	# list files in STL dir | remove the source STL from list | replace prior_ with a special char | sort numerically the second field in the STL file names 
	# after separating the STL file names based on the special char | replace the special char in the name of the STL files with prior_ to get back the original name of the files
	# In Linux you can have file names with the special chars: "}", "?", "<", "&", "!", "@"
	if [ "`ls $stlDirName | grep -v source | grep "}"`" == "" ]; then 
		sortedSTLFiles=`ls $stlDirName | grep "STL" | grep -v "source" | sed 's/prior_/}/g' | sort -t "}" -n -k 2 | sed 's/}/prior_/g'`
	
	elif [ "`ls $stlDirName | grep -v source | grep "?"`" == "" ]; then 
		sortedSTLFiles=`ls $stlDirName | grep "STL" | grep -v "source" | sed 's/prior_/?/g' | sort -t "?" -n -k 2 | sed 's/?/prior_/g'`
	
	elif [ "`ls $stlDirName | grep -v source | grep "<"`" == "" ]; then
		sortedSTLFiles=`ls $stlDirName | grep "STL" | grep -v "source" | sed 's/prior_/</g' | sort -t "<" -n -k 2 | sed 's/</prior_/g'`
	
	elif [ "`ls $stlDirName | grep -v source | grep "&"`" == "" ]; then
		sortedSTLFiles=`ls $stlDirName| grep "STL" | grep -v "source" | sed 's/prior_/\&/g' | sort -t "&" -n -k 2 | sed 's/\&/prior_/g'`
	
	elif [ "`ls $stlDirName | grep -v source | grep "!"`" == "" ]; then
		sortedSTLFiles=`ls $stlDirName | grep "STL" | grep -v "source" | sed 's/prior_/!/g' | sort -t "!" -n -k 2 | sed 's/!/prior_/g'`
	
	elif [ "`ls $stlDirName | grep -v source | grep "@"`" == "" ]; then
		sortedSTLFiles=`ls $stlDirName | grep "STL" | grep -v "source" | sed 's/prior_/@/g' | sort -t "@" -n -k 2 | sed 's/@/prior_/g'`

	else
		printf "***************\nFILE NAME ERROR\n***************\n\
Your STL file names contains at least one instance of all the following characters: }, ?, <, &, !, @.\n\
This script requires that all file names not contain at least one of the above characters to be executed properly.\n\
Please rename the STL files appropriately and try again. Exiting script.\n" 
		exit 2
	fi

	############################################
	#PRINT THE STL FILES IN THE REQUIRED FORMAT#
	############################################
	
	sortedSTLFilesArray=($sortedSTLFiles)
	i=1 #counter for the files listed in the geometry.in file

	printf "[GEOMETRY]\n#ID	STLfile    material_name    ignore_times\n$i       ignored    default\n" 

	#print out the sorted list of STL files to the console available for copy paste into the geometry.in file
	for f in "${sortedSTLFilesArray[@]}"; do 
		
		matType=`echo $f | sed 's/^.*mater_//' | sed 's/_prior.*//'`
		
		if [ "`echo $f | grep "mater"`" == "" ]; then  #if no material type is specified
			matType="MISSING_MATERIAL"
		fi
		
		if [ "${matType,,}" != "notused" ]; then #check that user wants to use this part in the simulation	
			i=$((i+2))
			printf "$i       $stlDirName/$f    $matType\n" 
		fi
	done

	#print the source volume section
	sourceVolume=`ls $stlDirName | grep "STL" | grep "prior_source"`
        sourceVolumeArray=($sourceVolume)
        printf "\n[SOURCE]\n"

        for g in "${sourceVolumeArray[@]}"; do
                printf "STLvolume       neutron       $stlDirName/$g    0    0\n"
        done
	
else
	printf "******************\nSTL DIR NAME ERROR\n******************\n\
Cannot find a directory of name $stlDirName.\nSpecify the name of the dir \
containing the STL files and try again. \nExiting script.\n" 
	exit 2
fi


