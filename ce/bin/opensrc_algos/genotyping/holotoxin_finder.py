###################################################################
#
# Reads based holotoxin finder
#
# Author: Milan Patel, with some direction from key mentors ;)
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import sys
import json
from functools import partial
from collections import namedtuple, defaultdict

from tools.environment import (
    log_message,
    log_error,
    log_progress,
    log_ephemeral,
    log_algo_version,
    write_results,
    full_path,
    valid_dir
)

from tools.dbinfo import (
    DbInfo,
    SequenceInfo
)

from tools.tools import (
    is_fasta,
    parse_fasta,
    popen
)