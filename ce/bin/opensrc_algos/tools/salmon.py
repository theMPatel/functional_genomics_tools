###################################################################
#
# Salmon runtime interface
#
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################


import os
import sys

from .tools import (
    popen,
    add_cmdline_args,
    add_cmdline_kwargs
)

from .environment import (
    log_message,
    log_error,
    log_warning,
    valid_dir,
    full_path
)

def run_salmon(env, *args, **kwargs):

    results_dir = os.path.join(env.localdir, 'salmon_results')

    salmon_path = os.path.join(
        env.toolsdir,
        'all_tools',
        'salmon'
    )

    valid_dir(results_dir)

    if not os.path.exists(salmon_path):
        raise RuntimeError('Missing salmon executable!')

    cmd_args = [
        salmon_path
    ]

    add_cmdline_args('salmon', cmd_args, args)
    add_cmdline_kwargs('salmon', cmd_args, kwargs)

    cmd_args.extend(['-o', results_dir])

    log_message('Running salmon args: {}'.format(
                ' '.join(cmd_args)))

    return_code, out, err = popen(cmd_args)

    if return_code:
        log_error(err.strip())
        raise RuntimeError('Error running salmon')
