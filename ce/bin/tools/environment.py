###################################################################
#
# Environment specifics for the genotyping algorithm
# You can use this for basically anything you want to
# implement into the Calculation Engine
# 
# Author: Milan Patel, with some direction from key mentors ;)
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import json
from copy import deepcopy
from datetime import datetime

__all__ = [
    'sanitize_path',
    'valid_dir',
    'check_dir',
    'time_now',
    'log_progress',
    'log_message',
    'log_error',
    'log_warning',
    'write_results',
    'Environment',
    'Logger',
    'ResultWriter'
    ]

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

def time_now():
    # Returns the locale specific datetime
    return datetime.now().strftime('%c')

def make_pretty(msg, depth):
    return ''.join(['\t']*depth) + msg + '\n'

def log_algo_version(algo_version=None, settings=None, env=None):

    if not all([algo_version, settings, env]):
        return

    database_version = '?'

    version_path = settings['version']

    if version_path is not None:

        version_path = env.get_sharedpath(version_path)

        with open(version_path, 'r') as f:

            version_info = json.load(f)

        database_version = version_info.get('version_info', '?')

    # Create the version strings for output
    algo_version_str = 'Using algorithm version: {}'.format(
        algo_version)

    algo_database_str = 'Using database version: {}'.format(
        database_version)

    # Log the poop out of them
    log_message(algo_version_str, 2)
    log_message(algo_database_str, 2)

def log_progress(msg, depth=0):
    if not Logger.current:
        print(msg)
    
    else:
        Logger.current.log_progress(msg)

def log_message(msg, depth=0):
    if not Logger.current:
        print('MSG\t' + make_pretty(msg, depth))

    else:
        Logger.current.log_msg(
            time_now() + '\tMSG\t' + make_pretty(msg, depth)
        )

def log_ephemeral(msg, depth=0):
    if not Logger.current:
        print('MSG\t' + make_pretty(msg, depth))    

    else:
        Logger.current.log_ephemeral(
            time_now() + '\tMSG\t' + make_pretty(msg, depth)
        )

def log_error(msg, depth=0):
    if not Logger.current:
        print('ERR\t' + make_pretty(msg, depth))

    else:
        Logger.current.log_error(
            time_now() + '\tERR\t' + make_pretty(msg, depth))

def log_warning(msg, depth=0):
    if not Logger.current:
        print('WRN\t' + make_pretty(msg, depth))

    else:
        Logger.current.log_warning(
            time_now() + '\tWRN\t' + make_pretty(msg, depth))

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

class Logger(object):
    # Provides the logging functionality of the program

    current = None

    def __init__(self, logdir):
        
        self._logdir = logdir

        # Run the setup
        self.setup()

        if Logger.current is None:
            Logger.current = self

    def setup(self):

        # Make sure that the path exists
        valid_dir(self._logdir)

        # The important files
        self._logfiles = {
            'messagesfile': 'messages.txt',
            'progressfile' : '__progress__.txt',
            'messagefile' : '__message__.txt',
            'errorfile' : 'errors.txt',
            'warningsfile': 'warnings.txt'
        }

        # Used to hold onto the files for logging
        self._filehandles = {}

        # Create the open file handles and save them
        for attr, filename in self._logfiles.iteritems():
            
            # Create the path of the file name
            path = os.path.join(self._logdir, filename)

            # Create the file handle
            if attr in ['progressfile', 'messagefile']:
                # These are write-over files, thus we close them
                self._filehandles[attr] = self.initialize_file(path, action='close')

            else:
                # These are append files, keep them open
                self._filehandles[attr] = self.initialize_file(path)

            # Set the attribute
            setattr(self, attr, self._filehandles[attr])


        # This needs to be at the start of the __messages__.txt
        # file
        self._log('#\n', ['messagesfile'])

    def initialize_file(self, path, action='open'):
        # Regardless of action, this function
        # hits the filesystem with a write of the requested
        # file
        if action == 'close':
            open(path, 'w').close()
            return path

        elif action == 'open':
            return open(path, 'w', 0)

    def _log(self, msg, files):
        # Writes to the file, make sure your out message
        # is ready to go before calling this method
        # make sure files is an iterable
        for file in files:

            attr = getattr(self, file)

            # __messages__ and progress files are
            # write from beginning of files
            # their attributes are the path to the file
            # rather than the open file handle

            if isinstance(attr, basestring):

                with open(attr, 'w') as f:
                    f.write(msg)

            else:
                attr.write(msg)

    @classmethod
    def close_files(cls):

        file_handles = cls.current._filehandles

        for file, handle in file_handles.iteritems():

            if not isinstance(handle, basestring):
                os.fsync(handle.fileno())
                handle.close()

    def log_progress(self, msg):
        # There can only be one thing in the progress file which is the
        # number relating to the progress
        if isinstance(msg, int):
            msg = str(msg)

        files = ['progressfile']

        self._log(msg, ['progressfile'])

    def log_error(self, msg):

        files = [
            'messagesfile',
            'messagefile',
            'errorfile'
        ]
        
        self._log(msg, files)

    def log_warning(self, msg):

        files = [
            'warningsfile'
        ]

        self._log(msg, files)

    def log_msg(self, msg):

        files = [
            'messagesfile',
            'messagefile'
        ]

        self._log(msg, files)

    def log_ephemeral(self, msg):

        files = [
            'messagefile'
        ]

        self._log(msg, files)

class ResultWriter(object):

    current = None

    def __init__(self, resultsdir):
        self._resultsdir = resultsdir

        valid_dir(self._resultsdir)

        if ResultWriter.current is None:
            ResultWriter.current = self

    def add_result(self, name, content):

        path = os.path.join(self._resultsdir, name)

        with open(path, 'w') as f:
            f.write(content)