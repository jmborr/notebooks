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
    '''Search the notebook file for pattern matchings. Also include the notebook itself
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
    files.append(os.path.basename(name))  # include notebook itself
    return files

def filterin(files, source):
    ''' Check files exist under source directory'''
    filtered = []
    for file in files:
        if checkexists(os.path.join(source,file),doexit=False):
            filtered.append(file)
    return filtered

def fileTooBig(files, source, size_limit):
    "Units are bytes"
    filtered = []
    for file in files:
        if os.stat(os.path.join(source,file)).st_size < size_limit:
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
    '''Copy the files from source to destination
    '''

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
parser.add_argument('name', help='ipython notebook file name')
parser.add_argument('source', help='source directory containing the files.')
parser.add_argument('destination', help='destination directory.')
extensions='agr,cif,conf,cpptraj,doc,docx,f,gif,in,jpeg,jpg,mol2,mpg,pbs,pdb,pdf,png,PNG,pptx,ptraj,py,sdf,sh,vmd,xml'
parser.add_argument('--extensions', help='only files with selected extensions will be fetched. Default: "{0}"'.format(extensions))
size_limit=10000000L #10MB
parser.add_argument('--sizelimit', help='fetch only files with size under limit; Units are bytes; Default={0}'.format(size_limit))
args=parser.parse_args()

args.name=os.path.join(args.source,os.path.basename(args.name))
for argument in (args.name, args.source, args.destination):
    checkexists(argument)

if not args.extensions:
    extensions = extensions.split(',')
else:
    extensions = args.extensions.strip().split(',')

if not args.sizelimit:
    args.sizelimit=size_limit

files = scan_notebook(args.name, extensions) # Find files in source
files = filterin(files, args.source) # Discard non-existing files in source
files = fileTooBig(files, args.source, args.sizelimit)
files = filterout(files, args.source, args.destination) # Rules to discard found files
update_destination(files, args.source, args.destination) # Transfer the files
