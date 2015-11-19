#!/bin/sh
#
#  File:  initlib.sh
#
#  Purpose: Create the aggregate library if it is not already present.
#
#  Input: The library which is copied to the aggregate library, if aggregate
#    library is to be created.

AGGREGATE_LIB=all.a

if [ -f ${AGGREGATE_LIB} ]
then
#   Just a dummy statement to simply fill in the "then" part
    test ${AGGREGATE_LIB}
else
#   Copy the input file to the aggregate library
    cp $1 ${AGGREGATE_LIB}
fi
