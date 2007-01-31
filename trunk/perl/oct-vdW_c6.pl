#!<PERL>

# Physical constants
my $P_BOHR    = 0.529177;
my $P_HARTREE = 2.0*13.60580;
my $P_C       = 137.036;
my $M_PI      = 4 * atan2(1,1);

if($#ARGV != 1){
  print "usage: $PROGRAM_NAME cross_section1 cross_section2\n";
  exit 1;
}

# read files
open(IN1, "<$ARGV[0]");
open(IN2, "<$ARGV[1]");

# skip headers (this also skips 0 frequency)
while(($_=<IN1>) && ($_=~/^#/)){};
while(($_=<IN2>) && ($_=~/^#/)){};

# read data
my @omega  = ();
my @sigma1 = ();
my @sigma2 = ();

my $jj = 0;
while($_=<IN1>){
  /^\s*(\S+)\s*(\S+)/;
  $omega [$jj] = $1/$P_HARTREE;         # convert to atomic units
  $sigma1[$jj] = $2/($P_BOHR*$P_BOHR);

  $_=<IN2>;
  /^\s*(\S+)\s*(\S+)/;
  die "cross section files are not compatible" if($1/$P_HARTREE != $omega [$jj]);
  $sigma2[$jj] = $2/($P_BOHR*$P_BOHR);

  $jj++;
}

close(IN1);
close(IN2);

# calculate c6 integral
my $sum = 0.0;
for($ii = 0; $ii<$jj; $ii++){
  $alpha1 = ($sigma1[$ii]/$omega[$ii]) * $P_C/(4.0*$M_PI);
  $alpha2 = ($sigma2[$ii]/$omega[$ii]) * $P_C/(4.0*$M_PI);

  $sum += $alpha1*$alpha2;
}

$sum *= $omega[0]*3.0/$M_PI;

print "C6 = $sum [a.u.]\n";
$sum *= $P_BOHR**6/$P_HARTREE;
print "C6 = $sum [eV A^6]\n";
