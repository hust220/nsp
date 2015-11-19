#!/bin/sh

# This script cleans test results modules

# To run the script do the following:
#    clean_test.sh
#
# Arguments:
#    None.

cd ..
ModFile=./local/modules.txt

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

    if [ -d $DirModName/test ]
    then
        echo
        echo "------------------------------------------------------------"
        echo "**** Cleaning test results in $DirModName ****"
        echo "------------------------------------------------------------"
        cd $DirModName
        make clean_test
        cd ..
    fi
done < $ModFile

exit 0

