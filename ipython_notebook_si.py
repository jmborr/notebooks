#!/usr/bin/env python

import os
import argparse
import sys
import re
from tempfile import mkstemp

from pdb import set_trace as tr

def checkexists(file, doexit=True):
    '''Check if file or directory exists, possibly exit if not'''
    if not os.path.exists(file):
        if doexit:
            raise IOError , '{0} not found'.format(file)
            sys.exit(1)
        else:
            return False
    return True

def scan_notebook(name, extensions):
    '''Search the notebook file for pattern matchings
    '''
    handle,tmpname = mkstemp(dir='/tmp', prefix='junk')
    os.system("cat {0}|tr -d '\' > {1}".format(name,tmpname)) # remove backslash plague
    s = open(tmpname).read()
    files = []
    for extension in extensions:
        files += re.findall('\(files/([^)]*\.{0})\)'.format(extension),s)  # Check local link, with format (files/link)
        matches = re.findall('src=.*files/(.*.{0}.*width)'.format(extension), s) 
        files += [ match.split(extension)[0]+extension for match in matches] # clean up backspaces mess
    os.system('/bin/rm {0}'.format(tmpname))
    files = list(set(files)) # remove duplicates
    return files

def filterin(files, source):
    ''' Check files exist under source directory'''
    filtered = []
    for file in files:
        if checkexists(os.path.join(source,file),doexit=False):
            filtered.append(file)
    return filtered

def filterout(files, source, destination):
    ''' Check files newer or non-existing in destination'''
    filtered = []
    for file in files:
        if not checkexists(os.path.join(destination,file),doexit=False):
            filtered.append(file)
        else:
            t1 = os.path.getmtime(os.path.join(source, file))
            t2 = os.path.getmtime(os.path.join(destination, file))
            if t1 > t2:
                filtered.append(file)
    return filtered

def update_destination(files, source, destination):
    '''Copy the files from source to destination'''

    for file in files:
        sourcepath = os.path.join(source,file)
        destpath =  os.path.join(destination,file)
        destdir = os.path.dirname(destpath)
        if not checkexists(destdir, doexit=False):
            os.system( '/bin/mkdir -p {0}'.format(destdir) )
        os.system( '/bin/cp {0} {1}'.format(sourcepath,destpath) )
    if files:
        print 'Transferred files:\n', '\n'.join(files), "\nend of list of updated files"
    else:
        print 'No files transferred'

parser = argparse.ArgumentParser(description='Update supporting files for the ipython notebook')
parser.add_argument('name', help='filepath to the ipython notebook')
parser.add_argument('source', help='source directory containing the files.')
parser.add_argument('destination', help='destination directory.')
extensions='agr,conf,gif,in,jpeg,jpg,mpg,pdb,pdf,png,py,vmd,mol2,sdf,cif,pptx'
parser.add_argument('--extensions', help='only files with selected extensions will be fetched. Default: "{0}"'.format(extensions))
args=parser.parse_args()

for argument in (args.name, args.source, args.destination):
    checkexists(argument)

if not args.extensions:
    extensions = extensions.split(',')
else:
    extensions = args.extensions.strip().split(',')

files = scan_notebook(args.name, extensions) # Find files in source
files = filterin(files, args.source) # Discard non-existing files in source
files = filterout(files, args.source, args.destination) # Discard identical files in destination
update_destination(files, args.source, args.destination) # Transfer the files
