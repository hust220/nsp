#!/bin/sh

# This script builds dependent modules

# To run the script do the following:
#    compile.sh
#
# Arguments:
#    None.

CompType=$1

if [ "$CompType" = "debug" ];
then
    Opt=-g
    DispInfo="Making ($CompType)"
else
    Opt=-O
    DispInfo="Making"
fi

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
            if [ "$ModName" = "etc" ];
            then
                DirModName=$ModName
            else
                DirModName=$ModName-$ModTag
            fi
        fi 
    fi

    echo
    echo "------------------------------------------------------------"
    echo "**** $DispInfo $DirModName ****"
    echo "------------------------------------------------------------"

    cd $DirModName

    case $ModName in 
        etc)
            ./platform.sh
            ;;
        *)
            (make "OPT=$Opt" install) || exit 1
            ;;
    esac
    cd .. 
done < $ModFile

exit 0

