#############################################################################
# 
# Set the environment variables here that are specific for this job
# 
# Algorithm: Genotyping
# Contact: mpatel5@cdc.gov
# 
# 
#############################################################################


# Set the environment variable so samtools can find the htslib
LD_LIBRARY_PATH="$TOOLSDIR/samtools_dir/lib:/opt/python/lib:/data/ceroot_dev:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH

# Make sure that we can find BWA and SAMTOOLS
PATH="$PATH:$TOOLSDIR"
export PATH