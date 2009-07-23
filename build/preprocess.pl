#!/usr/bin/env perl

my $fin         = $ARGV[0];
my $ndebug      = $ARGV[1] eq "no";
my $nlong_lines = $ARGV[2] eq "no";
my $nline_num   = $ARGV[3] eq "no";
my $nf90_forall = $ARGV[4] eq "no";

my $tmpfile = "/tmp/oct_preprocess.$$";

#regular expression to match opening and closing parentesis
my $NestedGuts = qr{
  (?{ local $d=0 }) # Set depth to 0
  (?>
    (?:                  
      [^()]+        # you could see some non-parenthesis text
    | \( (?{$d++})  #  ...increment the depth
    | \) (?(?{ $d!=0 }) (?{ $d-- }) | (?!))     #  ...in which case decrement the depth...
    )*
  )
  (?(?{$d!=0})(?!))
}x;

open(IN,  "<$fin");
open(OUT, ">$tmpfile");

my $ncount_forall;

while($_ = <IN>)
{
  if($ndebug){
      next if /push_sub/;
      next if /pop_sub/;
  }
  if($nlong_lines){
      s/\\newline/\n/g;
      s/\\cardinal/#/g;
  }
  if($nline_num){
      next if /^#/;
  }
  if($nf90_forall){
      if(m/^\s*forall\s*\(($NestedGuts)\)\s*(?<body>.*)$/xi){
	  # get rid of commas that are inside parenthesis
	  my $newstr = "";
	  my $dd = 0;
	  my $str = $1;
	  for($i=0; $i<length($str); $i++){
	      my $c = substr($str, $i, 1);
	      $dd++ if($c eq '(');
	      $dd-- if($c eq ')');
	      $newstr .= ($c eq ',' && $dd > 0) ? '#comma#' : $c;
	  }

	  my @loops = split(/,/x, $newstr);
	  $ncount_forall = $#loops + 1;
	  my $i;
	  for($i=0; $i<$ncount_forall; $i++){
	      $loops[$i] =~ s/:/,/;
	      $loops[$i] =~ s/#comma#/,/;
	      print OUT "do ", $loops[$i], "\n";
	  }
	  if($+{body}){
	      print OUT "$+{body}\n";
	      for($i=0; $i<$ncount_forall; $i++){
		  print OUT "end do\n";
	      }
	  }
	  next;
      }
      if(m/^\s*end\s+forall/i){
	  my $i;
	  for($i=0; $i<$ncount_forall; $i++){
	      print OUT "end do\n";
	  }
	  next;
      }
  }
  print OUT;
}

close(IN);
close(OUT);
`mv -f $tmpfile $fin`;
