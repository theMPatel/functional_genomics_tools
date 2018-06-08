#############################################################################
# 
# An installer for any algorithms you want to install into the CE
# Author: Milan Patel
# Date: 4/27/2018
# 
#############################################################################

import os
import pwd
import grp
import sys
import json
import shutil
import logging
import argparse
import traceback
import subprocess as sp
from functools import partial

_BASE_PATH = '/data'
_REPO = os.path.basename(os.getcwd())
_LOGGING_PATH = os.path.join(_BASE_PATH, 'deployment_logs')

# Define the log size in bytes
_MAX_LOG_SIZE = 10**8
_MAX_LOG_FILES = 10
_DEFAULT_LOGGING_LEVEL = logging.DEBUG
_LOGGER = None

if not os.path.exists(_LOGGING_PATH):
    os.makedirs(_LOGGING_PATH)

def parse_cmdline():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--env',
        type=basestring,
        required=True
    )

    parser.add_argument(
        '--root',
        type = basestring,
        required=True
    )

    args = parser.parse_args()

    return parser.root, parser.env

def setup_logging(
    fmt='%(asctime)s\t%(levelname)s\t%(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=_DEFAULT_LOGGING_LEVEL,
    environment='dev'):
    
    global _LOGGER
    _LOGGER = logging.getLogger()

    rotating_log_file_path = os.path.join(
        _LOGGING_PATH,
        _REPO,
        '{env}_{repo}_log.log'.format(
            env=environment,
            repo=_REPO
        )
    )

    handler = logging.handlers.RotatingFileHandler(
        rotating_log_file_path,
        maxBytes = _MAX_LOG_SIZE,
        backupCount = _MAX_LOG_FILES
    )
    handler_formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)
    handler.setFormatter(handler_formatter)

    _LOGGER.addHandler(handler)
    _LOGGER.setLevel(level)

def popen(args, stdout=None, stderr=None, cwd=None, shell=False):

    if not isinstance(args, list):
        raise RuntimeError('Provided arguments must be of type list')

    if not stderr:
        stderr = sp.PIPE

    if not stdout:
        stdout = sp.PIPE

    if not cwd:
        cwd=os.getcwd()

    child = sp.Popen(args, stdout=stdout, stderr=stderr, cwd=cwd, shell=shell)

    out, err = child.communicate()

    return child.returncode, out, err

def full_path(path):
    return os.path.abspath(os.path.realpath(path))

def check_py_version_range(low=(2,7), high=(2,9,9)):
    return low <= sys.version_info <= high

def sanitize_path(path):
    return path.replace('"', '').replace("'", "")

def change_perms(path, perm=0o755):
    if not os.path.exists(path):
        raise RuntimeError('Invalid file path: {}'.format(path))

    os.chmod(path, perm)

_USER = 'ce'
_GROUP = 'ce'
# _USER = 'nft9'
# _GROUP = 'users'
_UID = pwd.getpwnam(_USER).pw_uid
_GUID = grp.getgrnam(_GROUP).gr_gid

def change_owner(path, uid=_UID, guid=_GUID):
    os.chown(path, uid, guid)

def copy_func(src, dst):
    shutil.copy2(src, dst)

def dos2unix(filepath):
    dos2unix = ['dos2unix', filepath]
    if os.path.exists(filepath):
        results = popen(dos2unix)

    if results[0]:
        raise RuntimeError("Failed converting file: {}".format(filepath))

def git_tree(repo_directory):

    git_tree_args = [
        'git',
        'ls-tree',
        '-r',
        '--full-tree',
        '--name-only',
        'HEAD'
    ]

    results = popen(git_tree_args, cwd=repo_directory)

    if results[0]:
        raise RuntimeError('Failed to get repo tree')

    file_tree = results[1].split('\n')

    return file_tree

def get_version(repo_directory):

    git_version = ['git', 'describe']

    results = popen(git_version, repo_directory)

    if results[0]:
        raise RuntimeError('Error getting version tag: {}'.format(
            results[2]))

    return results[1]

def deploy(file_tree, environment, root, repo_directory):

    if environment == 'prod':
        environment = ''

    else:
        environment = '_' + environment

    root_path = root+environment

    ce_root_dir = os.path.join(
        _BASE_PATH,
        root_path
    )

    if not os.path.exists(ce_root_dir):
        raise RuntimeError('Invalid root directory for deployment: {}'.format(
            ce_root_dir
            )
        )

    for file in file_tree:

        # We only want to grab things that are from the
        # ce directory
        paths = file.split(os.sep)

        if paths[0] != 'ce' or paths[0] != 'shared':
            continue

        src = full_path(
            os.path.join(repo_directory, file)
        )

        paths[0] = ce_root_dir

        dst = os.path.join(*paths)

        copy_func(src, dst)
        dos2unix(dst)
        change_owner(dst)
        change_perms(dst)

    repo_version = get_version(repo_directory)
    version_dir = os.path.join(ce_root_dir, 'versions', _REPO)

    if not os.path.exists(version_dir):
        os.makedirs(version_dir)

    version_path = os.path.join(version_dir, 'version.txt')

    with open(version_path, 'w') as f:
        f.write(repo_version)

def main():

    root, env = parse_cmdline()
    setup_logging(environment=env)

    src_tree = git_tree(os.getcwd())

    deploy(src_tree, env, root, os.getcw())


if __name__ == '__main__':

    ret_val = 0
    try:
        main()

    except:
        ret_val = 1
        _LOGGER.exception('Error deploying code!')

    sys.exit(ret_val)