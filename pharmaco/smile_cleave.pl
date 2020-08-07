#!/usr/bin/perl -w
use strict;

# 1 August 2007
# Peter Ung @ U of Michigan, Medicinal Chemistry Department
# cleave the ZINC SMILES string file .smi into small pieces, with user defines
# the number of SMILES string in each file.

die "\nUsage: x.pl <number> <filename> <output index>
--> eg: x.pl 10,000 zinc.smi zinc6
    Default number: 25,000\n\n" if @ARGV != 3;

@_ = ();
my $num = 1;
my $fileNum = 0;
my $limit = 25000;
$limit = shift @ARGV;
$limit =~ s/,//g;
my $outName = pop @ARGV;

while (<>)
{
  if ($num > $limit)
  {
    print "Read ", $limit*($fileNum+1), "lines.\n";
    $num = 1;
    $fileNum++;
    close OUT;
  }
  my $file = $outName."_"."$fileNum".".smi";
  die "Cannot open $file" unless open OUT, ">> $file";
  if (/ZINC/)
  {
    $num++;
    print OUT;
  }
}
print "-->> ", ($fileNum+1)*$limit+$num, " lines read";
