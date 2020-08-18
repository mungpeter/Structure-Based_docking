#!/usr/bin/python

##########################################
##
##	Peter M.U. Ung @ MSSM
##
##	v1.0	- 13.10.30
##	v2.0	- 13.11.07 - added ability to ready multiple directories
##
##	v3.0	- 13.11.16 - this is a NEW script
##	v4.0	- 13.11.18 - if *.vina_score.txt is found, use it to check
##			     the redock result instead
##
##	Read in the individual vina-score files in each designated folder 
##	and generate a single vina-score file.
##	The finalized vina-score.txt will be put into the /result folder
##	This script will parse the vina-score.txt to look for molecules that
##	Vina has failed to dock the first run. 
##	
##	The individual vina-score files should have extension ".temp"
##	and the finalized vina-score file will have extension ".vina_score.txt"
##
##	## for versions before v3.0
##	## Preprocess the VINA docking results to generate a txt file with
##	## vina score and filename to speed up ranking process
##
##########################################

import sys,re,glob

if len(sys.argv) < 2:
  msg = """\n  ## Usage: x.py [directories with docked pdbqt files] ##
              e.g. /> x.py 21_p0.1 21_p0.2 21_p0.3
	      e.g. /> x.py "21_p0.*"
    # If previous result (*.vina_score.txt, *.vina_fail.txt) are found,
    # will check for the redock result instead (*.vina_fail.txt, *.redock.temp\n"""
  sys.exit(msg)


############################################################################
Dir_Names = []

sys.argv.pop(0)
for x in sys.argv: Dir_Names.append(glob.glob(x)[0])
print '  ## Found '+str(len(Dir_Names))+' directory ##'
print Dir_Names


def main ():
  for directory in Dir_Names:
    Files = []
 
  ## If found that the folder contains previously concatenated result file
  ## and Redock score file, then
    redock_file = glob.glob(directory+'/'+directory+'.redock.temp')
    if redock_file:
      print '  ## # Found ReDock result in: '+directory+' # ##'
      print redock_file

      fail_file  = glob.glob(directory+'/'+directory+'.vina_fail.txt')
      if not fail_file:
        print ' ##### Cannot find vina_fail in: '+directory+' Nothing is done ######'
      else:
        print '  ## # Found vina_fail in: '+directory+' # ##'
        print fail_file

      score_file = glob.glob(directory+'/'+directory+'.vina_score.txt')
      if not score_file:
        print ' ##### Cannot find vina_score in: '+directory+' Nothing is done ######'
      else:
        print '  ## # Found vina_score in: '+directory+' # ##'
        print score_file
        CheckRedock(directory, score_file[0], redock_file[0], fail_file[0])

  ## else, just concatenate all result file and tell how many docking failed
    else:
      Files = glob.glob(directory+'/*.temp')
      print '  ## Found '+str(len(Files))+' Vina Score files in: '+directory+' ##'
      CheckResult(directory, Files)    


############################################################################
def CheckResult(directory, Files):
  SCORE = open(directory+'/'+directory+'.vina_score.txt', 'w')
  FAIL  = open(directory+'/'+directory+'.vina_fail.txt',  'w')	# Vina didnt dock
  Failed = []
  print "    ## Compiling Vina Scores in "+directory+" ##"

  for temp_file in Files:
      with open(temp_file, 'rh') as f:
        for line in f:
	  if not re.search(r'REMARK VINA RESULT', line):
	    Failed.append(line.split('::')[0]+"\n")
          else:
            SCORE.write(line)
  print "      ## Generated "+directory+".vina_score.txt ##"
  print "      ## "+directory+" has "+str(len(Failed))+" ligands need to redock ##"
  for x in sorted(Failed): FAIL.write(x)
  SCORE.close()
  FAIL.close()

  print "\n### Remember to check the $folder.vina_fail.txt for missing ligands ###\n"


###########################################################################

## Compare the *.redock.temp result to the *.vina_fail.txt result
## If the result are not the same (different number of entry), exit 
def CheckRedock(directory, score_file, redock_file, fail_file):
  Scores = []
  Redock = []
  Failed = []
  Fail_Names = []

  with open(score_file, 'rh') as f:
    for line in f: Scores.append(line)

  with open(fail_file, 'rh') as f:
    for line in f:
      Fail_Names.append(line)

  with open(redock_file, 'rh') as f:
    for line in f:
      if re.search(r'REMARK VINA RESULT', line):
        Scores.append(line)
        Redock.append(line)
      else:
        Failed.append(line.split('::')[0]+"\n")

  if len(Fail_Names) != len(Redock):
    sys.exit("  ### Only "+str(len(Redock))+" Out of "+str(len(Fail_Names))+" are redocked. Redo. ###")
  if len(Fail_Names) > len(Redock):
    sys.exit("  ### Only "+str(len(Redock))+" Out of "+str(len(Fail_Names))+" are redocked. Redo. ###")
#  if len(Fail_Names) < len(Redock):
#    sys.exit("  ### Number of redock: "+str(len(Redock))+" is larger than number of failed in "+fail_file+": "+str(len(Fail_Names)+". Assumed redock.temp has already been incorporated into the vina_score.txt. Check to confirm. ###")

  SCORE = open(directory+'/'+directory+'.vina_score.txt', 'w')
  FAIL  = open(directory+'/'+directory+'.vina_fail.txt',  'w')

  for x in sorted(Scores): SCORE.write(x)
  for x in sorted(Failed): FAIL.write(x)
  SCORE.close()
  FAIL.close()
  print "\n  ## # Regenerated NEW Score file: "+score_file+" # ##\n"
  print "  ## # "+directory+" has "+str(len(Failed))+" Failed docking # ##\n"


#####################################################
if __name__ == '__main__':
    main()
