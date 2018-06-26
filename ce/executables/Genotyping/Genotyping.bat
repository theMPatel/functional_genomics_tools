#!/bin/sh
#
# Launcher for executables where the underlying DRM is OpenGridEngine (or associated)
# 
# This script provides a general launcher for binaries. Create a copy
# of this script, set MYNAME, and this script will invoke the
# corresponding binary a proper environment.

CE_STORAGE_CEROOT=/data/ceroot_dev
MYNAME=opensrc_algos/cewrapper/cewrapper.py

export SYSDIR=$CE_STORAGE_CEROOT

#*******************************************************************************
# You should not need to change anything below
#*******************************************************************************
. $SYSDIR/executables/setenv.sh

BINARY=$BINDIR/$MYNAME

echo "Number of slots allocated: $NSLOTS"
if [[ -s $PE_HOSTFILE ]]; then
    echo "The hosts are: "
    /bin/cat $PE_HOSTFILE
fi

# Save the node that this job is run on so if we need to triage
# we know exactly what machine this was run on

hostname >> "output.txt"

# Could check here that the hosts are all on the same machine
# i.e. $PE_HOSTFILE should consist of a single line. CalculationEngine
# executables cannot be run across machines

#*******************************************************************************
# Prepare the extra command line
#*******************************************************************************
cmdLineExtra=""
if [ ! -z "$NSLOTS" ]; then
    cmdLineExtra="$cmdLineExtra --nThreads=$NSLOTS"
fi

# Use a local directory on the computation node instead of doing temporary stuff on the network drive
if [ ! -z "$TMPDIR" ]; then
   cmdLineExtra="$cmdLineExtra --localdir=$TMPDIR"
fi

# toolsdir
if [ ! -z "$TOOLSDIR" ]; then
   cmdLineExtra="$cmdLineExtra --toolsdir=$TOOLSDIR"
fi

# Put your custom stuff in custom.sh, next to this file.
[ -x custom.sh ] && . custom.sh

#*******************************************************************************
# Now run the algorithm
#*******************************************************************************

# $@ is command line, writes cmd line to settings.txt
echo $@ $cmdLineExtra > settings.txt

EXECUTOR=$TOOLSDIR/all_tools/python2

# executes BINARY
$EXECUTOR -O $BINARY

exitCode=$?

DRMERR="errors.txt"
CEERR="logs/errors.txt"

if [ $exitCode -ne 0 -a -s "$CEERR" -a -f "$DRMERR" ]; then
        # append DRM errors to CE errors
        cat "$DRMERR" > "$CEERR"
fi

if [ -d "$JOBTMPDIR" ]; then
        # remove tempdir we have created
        rm -rf "$JOBTMPDIR"
fi
