#!/usr/bin/perl -w

$HOME  = $ARGV[0];	# home directory
$FOLD  = $ARGV[1];	# ligand folder directory (1 level up)
$RSLT  = $ARGV[2];	# home directory of result
$LIB   = $ARGV[3];      # ligand folder name
$TOTAL = $ARGV[4];	# total of ligand found in folder
$STEP  = $ARGV[5];	# how many ligand to run on each CPU
$PBS   = $ARGV[6];	# template pbs for vina run
$CONF  = $ARGV[7];
$QUEUE = $ARGV[8];
$WALL  = $ARGV[9];
$NODES = $ARGV[10];

$submit = $ARGV[11];

print "  ## $TOTAL ligands found in ".$FOLD."/".$LIB." ##\n";

$NUM   = 1+int( ($TOTAL/$STEP) + ($TOTAL/$STEP)/abs($TOTAL*2/$STEP) );
$START = 0;
$LAST  = 0;  
print "  ## Number of 'step' PBS will be run: ".$NUM." ##\n";

for ($i = 0; $i < $NUM; $i++) {

  $out = $LIB.".".$i.".pbs";
  open OUT, "> $out";
  open INP, "< $PBS";

  $LAST = $START + $STEP - 1;

  while (<INP>) {

    s/RUN/$i/;
    s/FOLD/$FOLD/;
    s/HOME/$HOME/;
    s/RSLT/$RSLT/;
    s/LIB/$LIB/;
    s/START/$START/;
    s/LAST/$LAST/;
    s/CONF/$CONF/;
    s/QUEUE/$QUEUE/;
    s/WALL/$WALL/;
    s/NODES/$NODES/;
    print OUT;
  }

  $START = $LAST + 1;

  close OUT;
  close INP;
  
  system("qsub $out") if $submit;
}
