#!/usr/bin/perl -w

# 5 August 2007 / 7 Aug 2007
# Writer: Peter Ung @ U of Michigan, Medicinal Chemistry Department
#
# Read in ZINC database SMILES strings file .smi and calculate the MW of
# each entry. Also recognize molecule with Aromatic carbon. Molecules with
# aromatic carbons, heteroatoms (O, N, S, P)  and fall in the selected MW 
# range are selected and output.
# WARNING: For S(+4) and P(+3) states the calculator cannot correctly 
#          calculate the MW and always off by + 2au.
#          Also, if simpified SMILES format is used for substituted 5-membered 
#          ring nitrogen, eg n1cccc1 instead of N1C=CC=C1, the calculated MW
#          will be off by -1au for each such atom.

die "\nUsage: x.pl <min M> <max M> <SMILE in> <out prefix> <optional out size>
-> Output molecules with Aromatic carbon in the selected MW range
   WARNING!! Cannot process S(+4) and P(+3) states (-2au off)
   WARNING!! Cannot process substituted 5-membered ring Nitrogen (-1au off)\n\n"
     if (@ARGV <=3 || @ARGV >= 6);

  # atomic mass
  $Cm = 12.011; $Om = 15.9994; $Nm = 14.00674;
  $Sm = 32.065; $Pm = 30.97376; $Hm = 1.0078;
  $Fm = 18.998403; $Clm = 35.453; $Brm = 79.904; $Im = 126.9045;

$num = 0;
$fileNum = 0;
$limit = 50000;

$limit = pop @ARGV if (@ARGV == 5);
$minM = shift @ARGV;
$maxM = shift @ARGV;
$prefix = pop @ARGV;

sub massCal
{
  # n, s, p represent the alternative states of N(-3), S(-2), P(-3)
  $C = $O = $N = $S = $P = $F = $Cl = $Br = $I = 0;
  $n = $s = $p = 0;
  $arom = $aroN = $doubleBond = $cycl = $pos = $neg = 0;

  @smile = @_;
  $explicitH = pop @smile;
  $end = $#smile;

  for ( $i = 0; $i <= $end; $i++ )
  {
    # C
    $C++ if ($smile[$i] =~ /C/);
    $C++ if ($smile[$i] =~ /c/);
    $O++ if ($smile[$i] =~ /O|o/);

    # N -- assume -3 state, find +5 state
    if ($smile[$i] =~ /N/)
    { 
      # if N in +5 state N+(=O)O-
      if (($i+1) <= $end && $smile[$i+1] eq "+")
      {
        if (($i+2) <= $end && $smile[$i+2] =~ /\(/)
        {
          if ($smile[$i+3] =~ /=/ && $smile[$i+4] =~ /O/)
          { $n++; }
          else
          { $N++; }
        }
        else
        { $N++; }
      }
      else
      { $N++; }
    }
    if ($smile[$i] =~ /n/)
    {
      $N++;
      if (($i+1) <= $end && $smile[$i+1] eq /\(/)
      { $aroN++; }
    }

    # S -- assume -2 state, find +6 state;
    # WARNING -- +4 state will be counted as +6 state
    if ($smile[$i] =~/S|s/)
    {
      if (($i+1) <= $end && $smile[$i+1] =~ /\(/) 
      {
        # if S in +6 state
        if ($smile[$i+2] =~ /=/ && $smile[$i+3] =~ /O|N/)
        { $s++; }
        else
        { $S++; }
      }
      # if S is S2(=O)
      elsif (($i+1) <= $end && $smile[$i+1] =~ /\d/)
      {
        if (($i+2) <= $end && $smile[$i+2] =~ /\(/)
        {
          if ($smile[$i+3] =~ /=/ && $smile[$i+4] =~ /O|N/)
          { $s++; }
          else
          { $S++; }
        }
        else
        { $S++; }
      }
      else
      { $S++; }
    }

    # P -- assume -3 state, find +5 state;
    # WARNING -- +3 state will be counted as +5 state
    if ($smile[$i] =~ /P|p/)
    {
      # if P is in +5 state, P(=O)
      if (($i+1) <= $end && $smile[$i+1] =~ /\(/)
      {
        if ($smile[$i+2] =~ /=/ and $smile[$i+3] =~ /O|S|N/)
        { $p++; }
        else
        { $P++; }
      }
      # if P2(=O)
      elsif (($i+1) <= $end && $smile[$i+1] =~ /\d/)
      {
        if (($i+2) <= $end && $smile[$i+2] =~ /\(/)
        { 
          if ($smile[$i+3] =~ /=/ && $smile[$i+4] =~ /O|N|S/)
          { $p++; }
          else
          { $P++; }
        }
        else
        { $P++; }
      }
      else
      { $P++; }
    }

    # halogens
    $F++ if ($smile[$i] =~ /F/);
    $I++ if ($smile[$i] =~ /I/);
    $Br++ if ($smile[$i] =~ /B/ && ($i+1) <= $end && $smile[$i+1] =~ /r/);
    if ($smile[$i] eq "C")
    { 
      if (($i+1) <= $end && $smile[$i+1] =~ /l/)
      { $Cl++; $C--; }
    }
    $C-- if $smile[$i] =~ /C/ && ($i+1) <= $end && $smile[$i+1] =~ /a|d|e|m|s|u|o/;

    # degree of unsaturation
    $doubleBond++ if ($smile[$i] =~ /=/);
    $doubleBond += 2 if ($smile[$i] eq "#");
    $arom++ if ($smile[$i] =~ /c|n|p/);
    $cycl++ if ($smile[$i] =~ /\d/);

    # charges
    $pos++ if ($smile[$i] eq "+");
    $neg-- if ($smile[$i] eq "-");

  }  #end of loop

  $H = $C*2 + 2 -($F + $Cl + $Br + $I + $cycl + $arom + ($doubleBond*2)) +
       $pos + $neg + $explicitH +
       ($N)*2 + ($n-$N) + $aroN +
       ($s*4) +
       ($P+$p)*2 + ($p-$P);
  $mass = $C*$Cm + $O*$Om + ($N+$n)*$Nm + ($S+$s)*$Sm + ($P+$p)*$Pm +
          $F*$Fm + $Cl*$Clm + $Br*$Brm + $I*$Im +
          $H*$Hm;

#  $form = "C".$C."H".$H;
#  $form .= "N".($N+$n) if $N > 0;
#  $form .= "O".$O if $O > 0;
#  $form .= "S".($S+$s) if $S > 0;
#  $form .= "P".($P+$p) if $P > 0;
#  $form .= "F".$F if $F > 0;
#  $form .= "Cl".$Cl if $Cl > 0;
#  $form .= "Br".$Br if $Br > 0;
#  $form .= "I".$I if $I > 0;

  return $mass;
} #end of subroute &massCal

sub explicitH
{
  # count explicit H of aromatic N and P
  $countH = 0;

  for ($i = 0; $i <= $#_; $i++)
  {
    if ($_[$i] =~ /n|p/)
    {
      $countH++ if (($i+1) <= $#_ && $_[$i+1] =~ /H/)
    }
    $_[$i] = Z if ($_[$i] =~ /H/);
  }
  push @_, $countH;
  return @_;
}

$file = $prefix."_".$fileNum.".smi" if ($num == 0);
@_ = ();
while (<>)
{ 
  if ($num >= $limit)
  {
    $num = 0;
    $fileNum++;
    close OUT;
    $file = $prefix."_".$fileNum.".smi";
  }
  die "Cannot open $file" unless open OUT, ">> $file";
  if (/ZINC/)
  {
    chomp;

    if ($_ =~ /c/)  #select those with Aromatic C
    {
      @line = split;
      # select for compd with heteroatoms -- at least soluble
      if ($line[0] =~ /N|n|O|o|S|s|P|p/)
      {
        # remove charge number in [] eg [NH2+]
        $line[0] =~ s/(\[\w*)[0-9](.\])/$1$2/g;

        # remove auxilary characters, [] @ / \
        $line[0] =~ s/\[|\]|@|\/|\\//g;

        # preserve H of aromatic N and S and P
        @temp = split //, $line[0];
        @char = &explicitH(@temp);

        $mass = &massCal(@char);
#        print OUT "$_ ($mass) $explicitH\n";   # debugging check
        unless ($mass < $minM || $mass > $maxM)
        {
          print OUT "$_\n";        
          $num++;
        }
      }
    }
  }
}
