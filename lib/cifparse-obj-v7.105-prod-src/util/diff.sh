#!/bin/sh

# This script does the diff operation on all modules.

# To run the script do the following:
#    diff.sh
#
# Arguments:
#    None.

cd ..

ModFile=./local/modules.txt

echo
echo -------- Diffing the modules --------
echo


while read One Two Three Four; do
    if [ "$One" = "cvs" ];
    then
        Rep=$One
        ModName=$Three
        ModTag=$Four
    else
        if [ "$One" = "svn" ];
        then
            Rep=$One
            ModName=$Three
            ModTag=$Four
        else
            Rep=cvs 
            ModName=$One
            ModTag=$Two
        fi
    fi

    echo -------- Diffing version $ModTag of module $ModName --------

    if [ "$Rep" = "cvs" ];
    then
        DirModName=$ModName
    else
        if [ "$ModTag" = "Latest" ];
        then
            DirModName=$ModName
        else
            DirModName=$ModName-$ModTag
        fi 
    fi

    $Rep diff $DirModName

done < $ModFile

exit 0

