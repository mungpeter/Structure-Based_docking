#!/usr/bin/env python3
import sys,glob,re,gzip,bz2

## Handle gzip and bzip2 file if the extension is right. otherwise, just open
## outuput: file handle
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'r')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'r')
  else:
    handle = open(file_name, 'r')

#  print "## Opening "+file_name
  return handle


##########################################################################
## remove "remarks" and "empty line" from read-line
## output: list of strings
def remove_remark(handle):
  New_File = []
  # Read each line as a list of characters
  with handle as fi:
    Lines = filter(None, (l.rstrip() for l in fi))
    for line in Lines:
      Line = []
      if re.search(r'^#', line): continue
      for char in list(line.rstrip()):
        if char != '#': Line.append(char)
        else: break
      if len(Line) > 0:
        New_File.append(''.join(Line))
  return New_File


##########################################################################
def cmp(x, y):
  """
    Replacement for built-in function cmp that was removed in Python 3

    Compare the two objects x and y and return an integer according to
    the outcome. The return value is negative if x < y, zero if x == y
    and strictly positive if x > y.
  """

  return (x > y) - (x < y)