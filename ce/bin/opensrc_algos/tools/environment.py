###################################################################
#
# Environment specifics for the genotyping algorithm
# You can use this for basically anything you want to
# implement into the Calculation Engine
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import json
import logging
import inspect
import collections
import logging.handlers
from copy import deepcopy
from datetime import datetime

__all__ = [
    'sanitize_path',
    'valid_dir',
    'check_dir',
    'log_progress',
    'log_message',
    'log_error',
    'log_warning',
    'write_results',
    'Environment',
    'ResultWriter'
    ]

def get_stack_len():
    return len(inspect.stack())-1

def get_message_depth(base_depth, extra=0):
    """
    Assumes that this is an internal logging function
    that the call stack looks like:

    worker -> log message -> check message depth
    """
    curr_depth = get_stack_len()-2

    if base_depth >= 0:
        base_depth = -curr_depth

    return {'depth':curr_depth+base_depth+extra}

def sanitize_path(path):
    # Required to remove quotes from bionumerics
    # caller script on Calculation Engine
    return path.replace('"', '').replace("'", "")

def valid_dir(path):
    # Check to make sure the directory exists
    if os.path.exists(path):
        return

    if not os.path.isdir(path):
        os.makedirs(path)

def check_dir(path):
    return os.path.isdir(path)

def full_path(path):
    return os.path.abspath(os.path.realpath(path))

base_depth = 0
def log_algo_version(algo_version=None, settings=None, env=None):

    database_version = '?'

    version_path = settings['version']

    if version_path is not None:

        version_path = env.get_sharedpath(version_path)

        with open(version_path, 'r') as f:

            version_info = json.load(f)

        algo_version = version_info.get('algorithm', algo_version)
        database_version = version_info.get('database', '?')

    # Create the version strings for output
    algo_version_str = 'Using algorithm version: {}'.format(
        algo_version)

    algo_database_str = 'Using database version: {}'.format(
        database_version)

    depth = get_message_depth(base_depth)

    # Log the poop out of them
    logging.info(algo_version_str, {'depth' : depth})
    logging.info(algo_database_str, {'depth' : depth})

def log_progress(msg):
    logging.log(logging.NOTSET, msg)

def log_message(msg, extra=0):
    depth = get_message_depth(base_depth, extra)
    logging.info(msg, depth)

def log_error(msg, extra=0):
    depth = get_message_depth(base_depth, extra)
    logging.error(msg, depth)
    logging.info('There was an error, check logs!', depth)

def log_warning(msg, extra=0):
    depth = get_message_depth(base_depth, extra)
    logging.warning(msg, depth)

def log_exception(msg, extra=0):
    depth = get_message_depth(base_depth, extra)
    logging.exception(msg, depth)

def write_results(name, content):
    if not ResultWriter.current:
        print('{} -> {}'.format(name, content))

    else:
        ResultWriter.current.add_result(name, content)

class Environment(object):

    def __init__(self, settings):

        if settings is None:
            raise RuntimeError('Missing settings for environment'
                ' creation.')

        # These are the settings we should be able to expect
        self._attrs = ['_resultsdir', '_localdir',
            '_shareddir', '_toolsdir', '_tempdir']

        # Run the setup tasks
        self.setup(settings)

    def get_sharedpath(self, path):
        # Make sure to replace the environment location
        path = path.replace('[SHAREDDIR]', os.path.dirname(self.shareddir))
        
        # normalize the path, just in case
        path = os.path.normpath(path)

        return path

    def copy(self):

        # Use this to dynamically create a new environment class
        # you can use this to change the paths for specific modules
        # like making a specific tmp directory for each
        # module
        return deepcopy(self)

    def setup(self, settings):
        # For each of the attr, set up the variables
        # to accept new values
        for attr in self._attrs:
            setattr(self, attr, None)

        # Get the paths for the environment
        for setting in settings:

            if hasattr(self, '_'+setting):

                # Get the path
                path = sanitize_path(settings[setting])
                
                # Make the directories
                valid_dir(path)

                # Set the attributes
                setattr(self, '_'+setting, path)

            elif hasattr(self, setting):

                # Get the path
                path = sanitize_path(settings[setting])

                # Make sure it exists
                valid_dir(path)

                # Set the attribute
                setattr(self, setting, path)


        # Make sure that all of the settings that were provided
        # are 'real'
        for attr in self._attrs:
            if getattr(self, attr) is None:
                raise RuntimeError('Missing environment variable: {}'.format(
                    attr))

        # Make sure we set ourselves up for logging
        self._logdir = os.path.join(self._resultsdir, 'logs')

        # Create the correct resultsdir
        self._resultsdir = os.path.join(self._resultsdir, 'results', 'raw')

        # Make sure it exists
        valid_dir(self._logdir)

        # Set the max thread out allocated to this task
        if 'nThreads' in settings:
            # The number of threads we can use.
            # You need to remember to subtract one
            # if you are going to let other
            # subprocesses take from this pool
            # to account for the genotyping algorithm itself
            # only when you are calling your subprocess though ;) 
            self._threads = int(settings['nThreads'])

            if self._threads <= 1:
                self._threads = 2

        else:
            self._threads = 2

    @property
    def localdir(self):
        return self._localdir

    @localdir.setter
    def localdir(self, path):
        valid_dir(path)
        self._localdir = path

    @property
    def shareddir(self):
        return self._shareddir

    @shareddir.setter
    def shareddir(self, path):
        check_dir(path)
        self._shareddir = path

    @property
    def toolsdir(self):
        return self._toolsdir

    @property
    def logdir(self):
        return self._logdir

    @property
    def resultsdir(self):
        return self._resultsdir

    @property
    def threads(self):
        return self._threads

    @property
    def tempdir(self):
        return self._tempdir

class SingleWriteFileHandler(logging.FileHandler):
    """
    This class is to be used in order to write to
    the progress file which is a single
    write file (only a single message can be in this file)
    """
    def emit(self, record):
        """
        Emit a record

        Because we want to overwrite the file everytime we 
        emit, we will always reopen the stream before calling
        StreamHandler.emit
        """

        self.stream = self._open()
        logging.StreamHandler.emit(self, record)

class ProgressFilter(object):
    """
    The progress file can only contain a number
    thus if we log anything besides a single number
    it should be thrown out for this handler
    """

    def filter(self, record):
        """
        *** FROM THE DOCS ***
        Determine if the specified record is to be logged.
        Is the specified record to be logged? Returns 0 for no, nonzero for
        yes. If deemed appropriate, the record may be modified in-place.
        """
        try:
            float(record.getMessage())
        except:
            return 0
        else:
            return 1

class TabbedModifier(object):
    """
    We can use this class to get the pretty
    depth information into the log files by modifying
    the log record at the filtering stage
    """

    def filter(self, record):
        """
        Determine if the the log record has a depth parameter in its
        argument list
        """

        if isinstance(record.args, collections.Mapping) and 'depth' in \
            record.args:

            new_msg = ''.join(['\t']*record.args['depth']) + \
                record.getMessage()
            
            record.msg = new_msg

        return 1

def bn_file_prep(filename, mode='a', **kwargs):
    """
    This function adds a hash to the start of the file

    Don't know the reasoning but bionumerics requires a hash
    symbol at the start of a log file to know that it needs
    to start parsing the messages 
    
    """
    with open(filename, mode) as f:
        f.write('#\n')

_base_formatter = logging.Formatter(
    fmt='%(asctime)s\t%(levelname)s\t%(message)s',
    datefmt='%Y-%m-%d %I:%M:%S %p'
)

_error_formatter = logging.Formatter(
    fmt='%(asctime)s\t%(levelname)s\n'
        '%(filename)s\t%(lineno)d\t%(funcname)s\n'
        '->%(message)s',

    datefmt='%Y-%m-%d %I:%M:%S %p'
)

base_logfiles = {

    'messages': {
        'handler'   : logging.FileHandler,
        'filename'  : 'messages.txt',
        'level'     : logging.INFO,
        'filter'    : TabbedModifier(),
        'format'    : None,
        'file_prep'  : bn_file_prep
    },
    'progress' : {
        'handler'   : SingleWriteFileHandler,
        'filename'  : '__progress__.txt',
        'level'     : logging.NOTSET,
        'filter'    : ProgressFilter(),
        'format'    : None,
        'file_prep'  : None
    },
    'single_message': {
        'handler'   : SingleWriteFileHandler,
        'filename'  : '__message__.txt',
        'level'     : logging.INFO,
        'format'    : None,
        'filter'    : None,
        'file_prep'  : None
    },
    'error': {
        'handler'   : logging.FileHandler,
        'filename'  : 'errors.txt',
        'level'     : logging.ERROR,
        'format'    : _error_formatter,
        'filter'    : None,
        'file_prep'  : None
    },
    'warnings' {
        'handler'   : logging.FileHandler,
        'filename'  : 'warnings.txt',
        'level'     : logging.WARNING,
        'format'    : None,
        'filter'    : None,
        'file_prep'  : None
    }
}

def initialize_logging(log_dir):

    # Make sure the dir exists
    valid_dir(log_dir)

    root = logging.getLogger()

    for handler, parameters in base_logfiles.iteritems():

        file_path = os.path.join(
            log_dir,
            parameters['filename']
        )

        if parameters['file_prep']:
            if callable(parameters['file_prep']):
                parameters['file_prep'](file_path)

        handler = parameters['handler'](file_path)
        
        handler.setLevel(parameters['level'])

        if parameters['format']:
            handler.setFormatter(parameters['format'])
        else:
            handler.setFormatter(_base_formatter)

        if parameters['filter']:
            handler.addFilter(parameters['filter'])


        root.addHandler(handler)

class ResultWriter(object):

    current = None

    def __init__(self, resultsdir):
        self._resultsdir = resultsdir

        valid_dir(self._resultsdir)

        ResultWriter.current = self

    def add_result(self, name, content):

        path = os.path.join(self._resultsdir, name)

        with open(path, 'w') as f:
            f.write(content)