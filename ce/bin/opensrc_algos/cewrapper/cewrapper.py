###################################################################
#
# Wrapper for the genotyping algorithm (or really anything)
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import json
import argparse 
import importlib
import os
import sys
import traceback

this_file = os.path.realpath(__file__)

base_path = os.path.dirname(
    os.path.dirname(this_file)
)

if base_path not in sys.path: 
    sys.path.append(base_path)

from tools.environment import (
    Environment,
    initialize_logging,
    graceful_shutdown_logging,
    ResultWriter,
    log_message,
    log_progress,
    log_error,
    log_exception,
    log_algo_params
)

from tools.custom_parser import CustomParser
from tools.config import Settings

if __debug__:
    import pdb

class MyArgumentParser(argparse.ArgumentParser):

    # Override to split space based cmdline args
    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()

# Parse the commandline arguments
def parse_cmdline():

    # create the parser
    parser = MyArgumentParser(fromfile_prefix_chars='@')

    # The environment variables
    parser.add_argument('--clientVersion', type=int)

    parser.add_argument('--nThreads',
        help='Number of threads', type=int, default=2)

    parser.add_argument('--localdir', default='./local',
        help='Local working directory', type=str)

    parser.add_argument('--tempdir',
        help='Temporary shared directory', type=str)

    parser.add_argument('--resultsdir',
        help='Results directory', type=str)

    parser.add_argument('--shareddir',
        help='Shared directory', type=str)

    parser.add_argument('--toolsdir',
        default='', help='Tools directory', type=str)

    parser.add_argument('--algorithm',
        default='', help='The module to be run', type=str)

    # We can now send the settings as a json file:
    genotyper_settings_path = os.path.join(
        os.getcwd(), 'genotyper_settings.json')

    genotyper_settings = []
    if os.path.exists(genotyper_settings_path):
        with open(genotyper_settings_path, 'r') as f:
            genotyper_settings = json.load(f)

    # TODO REPLACE THIS FUNCTIONALITY
    # 
    # This will store all of the args globally in
    # the CustomParser class
    CustomParser(genotyper_settings)

    # Get the known arguments of the command line
    try:
        args, remaining = parser.parse_known_args(['@settings.txt'])
    except:
        args, remaining = parser.parse_known_args()

    #return the arguments object
    return args, remaining

def main_throw_args(args, remaining, execable_modules):
    # Actually calls the genotyping algorithm
    
    # Initialize the environment
    env = Environment(vars(args))

    # Set up the logger
    initialize_logging(env.logdir)

    # Set up the results writer
    ResultWriter(env.resultsdir)

    # Log the called parameters
    log_message('Initializing cewrapper..')

    # Write the args out
    log_algo_params(vars(args))

    # The name of the genotyping module
    module_name = execable_modules.modules[args.algorithm]

    # Module specific configurations...
    if hasattr(execable_modules, args.algorithm):
        module_settings = getattr(execable_modules, args.algorithm)
    else:
        module_settings = {"version": None}

    # Import the module
    module = importlib.import_module(module_name, base_path)
    
    # Log the progress
    log_progress(0)

    #run the module's main function
    module.main(args, remaining, env, module_settings)
        
    #and we are done!
    log_progress(100)
    log_message('Done running algorithm: {}!'.format(
        args.algorithm))    

def main_throw():
        
    # Parse cmdline arguments
    args, remaining = parse_cmdline()

    # The config path that will load the module that
    # is requested
    config_path = os.path.join(base_path, 'config.json')

    # Load the settings
    execable_modules = Settings(settingspath=config_path)
    
    # Run the main program with arguments
    main_throw_args(args, remaining, execable_modules)

def main():

    return_code = 0

    try:
        main_throw()

    except Exception:
        log_exception('')
        return_code = 1

    finally:
        return return_code

if __name__=='__main__':
    # Get the return value of the algorithm
    retval = main()

    # Clear our logging contexts
    graceful_shutdown_logging()

    # Exit
    sys.exit(retval)
