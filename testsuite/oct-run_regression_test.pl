#!/usr/bin/env perl
#
# Copyright (C) 2005-2014 H. Appel, M. Marques, X. Andrade, D. Strubbe
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# $Id$

use warnings;
use Getopt::Std;
use File::Basename;
use Fcntl ':mode';
use Time::HiRes qw(gettimeofday tv_interval);
use Scalar::Util qw(looks_like_number);

sub usage {

  print <<EndOfUsage;

 Copyright (C) 2005-2014 H. Appel, M. Marques, X. Andrade, D. Strubbe

Usage: oct-run_regression_test.pl [options]

    -n        dry-run
    -v        verbose
    -h        this usage
    -D        name of the directory where to look for the executables   
    -s        exec suffix for the executables
    -f        filename of testsuite [required]
    -p        preserve working directories
    -l        copy output log to current directory
    -m        run matches only (assumes there are work directories)

Exit codes:
    0         all tests passed
    1..253    number of test failures
    254       test skipped
    255       internal error

Report bugs to <octopus-devel\@tddft.org>
EndOfUsage

  exit 0;
}


sub set_precision{
  my $p = $_[0];
  if($p ne "default"){
    $precnum = 1.0*$p;
  }elsif($_[1] =~ m/_single/){
    $precnum = 0.001
  } else {
    $precnum = 0.0001
  }
}

# Check whether STDOUT is a terminal. If not, no ANSI sequences are
# emitted.
if(-t STDOUT) {
    $color_start{blue}="\033[34m";
    $color_end{blue}="\033[0m";
    $color_start{red}="\033[31m";
    $color_end{red}="\033[0m";
    $color_start{green}="\033[32m";
    $color_end{green}="\033[0m";
} else {
    $color_start{blue}="";
    $color_end{blue}="";
    $color_start{red}="";
    $color_end{red}="";
    $color_start{green}="";
    $color_end{green}="";
}

if (not @ARGV) { usage; }

getopts("nlvhD:c:f:s:pm");

# avoid warnings 'used only once: possible typo'
$useless = $opt_h;
$useless = $opt_l;

# Default values
use File::Temp qw/tempdir/;

# Handle options
$opt_h && usage;

my $exec_directory;
if($opt_D) {
 $exec_directory = $opt_D;
 if($exec_directory !~ /^\//){
  $exec_directory = get_env("PWD")."/$exec_directory";
 }
} else {
 $exec_directory = "/usr/bin";
}

if(length($opt_f) == 0) {
    die255("ERROR: You must supply the name of a test file with the -f option.\n");
}

# Find out which executables are available.
opendir(EXEC_DIRECTORY, $exec_directory) || 
 die255("ERROR: Could not open the directory $exec_directory to look for executables");
@octopus_execs = grep { /^oct/ } readdir(EXEC_DIRECTORY);
closedir(EXEC_DIRECTORY);

# determine exec suffix
$exec_suffix = "";
if($opt_s)  { $exec_suffix = $opt_s; }

$aexec = get_env("EXEC");
$global_np = get_env("OCT_TEST_MPI_NPROCS");

$mpiexec = get_env("MPIEXEC");
$machinelist = get_env("MACHINELIST");
if ("$mpiexec" eq "") { $mpiexec = `which mpiexec 2> /dev/null`; }
chomp($mpiexec);

# mpiexec without arguments (to check if it is available)
$mpiexec_raw = $mpiexec;
$mpiexec_raw =~ s/\ (.*)//;

if ("$mpiexec_raw" ne "") {
    if(!( -e "$mpiexec_raw")) {
	print "mpiexec ($mpiexec_raw) does not exist\n";
    } elsif(!( -x "$mpiexec_raw")) {
	print "mpiexec ($mpiexec_raw) is not executable\n";
    }
}

# default dummy value, for dry run so MPIEXEC need not be set properly
if("$mpiexec" eq "") {
    $mpiexec = "mpiexec";
}

# default number of processors for MPI runs is 2
$np = 2;


# Figure out which are the executables to test
my @executables;
if($opt_n) {
    @executables = "octopus";
} else {
    find_executables();
}

# This variable counts the number of failed testcases.
$failures = 0;

$tempdirpath = get_env("TEMPDIRPATH");
if ("$tempdirpath" eq "") { $tempdirpath = '/tmp'; }
if (! -d $tempdirpath) { mkdir $tempdirpath; }

# Loop over all the executables.
foreach my $octopus_exe (@executables){

  set_precision("default", $octopus_exe);
  $test_succeeded = 1;

  $pwd = get_env("PWD");
  if (!$opt_m) {
      $workdir = tempdir("$tempdirpath/octopus.XXXXXX");
      chomp($workdir);

      system ("rm -rf $workdir");
      mkdir $workdir;
      
      $scriptname = "$workdir/matches.sh";
      open(SCRIPT, ">$scriptname") or die255("ERROR: could not create '$scriptname'.\n");
      print SCRIPT "#\!/usr/bin/env bash\n\n";
      print SCRIPT "perl $pwd/$0 -m -D $exec_directory -f $pwd/$opt_f\n";
      close(SCRIPT);
      chmod 0755, $scriptname;
      
      $matchdir = $workdir;
  } else {
      $workdir = $pwd;
  }

  # testsuite
  open(TESTSUITE, "<".$opt_f ) or die255("ERROR: cannot open testsuite file '$opt_f'.\n");

  $command = $octopus_exe;

  while ($_ = <TESTSUITE>) {

    # remove trailing newline 
    chomp; 
    # remove trailing whitespace 
    $_ =~ s/^\s+//; 
    # skip blank lines
    next if (length($_) == 0);

    # skip comments
    next if /^#/;

    if ( $_ =~ /^Test\s*:\s*(.*)\s*$/) {
      $test{"name"} = $1;
      print "$color_start{blue} ***** $test{\"name\"} ***** $color_end{blue} \n\n";
      print "Using workdir    : $workdir\n";
      if($opt_p) {
	  print "Workdir will be saved.\n";
      }
      print "Using executable : $octopus_exe\n";
      print "Using test file  : $opt_f \n";
    } elsif ( $_ =~ /^Enabled\s*:\s*(.*)\s*$/) {
      %test = ();
      $enabled = $1;
      $enabled =~ s/^\s*//;
      $enabled =~ s/\s*$//;
      $test{"enabled"} = $enabled;

      if ( $enabled eq "No") {
          print STDERR "Test disabled: skipping test\n\n";
	  if (!$opt_p && !$opt_m) { system ("rm -rf $workdir"); }
	  exit 254;
      } elsif ( $enabled ne "Yes") {
	  die255("ERROR: Unknown option 'Enabled = $enabled' in testsuite file.\n\n");
	  if (!$opt_p && !$opt_m) { system ("rm -rf $workdir"); }
      }
    } elsif ( $_ =~ /^Programs/) {
        # handled earlier
    } elsif ( $_ =~ /^Options/) {
        # handled earlier
    } elsif ( $_ =~ /^TestGroups/) {
        # handled by oct-run_testsuite.sh
    } else {
      if ( $enabled eq "") {
	die255("ERROR: Testsuite file must set Enabled tag before another (except Test, Programs, Options, TestGroups).\n\n");
      }

      if ( $_ =~ /^Util\s*:\s*(.*)\s*$/) {
	$command = "$exec_directory/$1";
	if( ! -x "$command") {
	  $command = "$exec_directory/../utils/$1";
	}

	if( ! -x "$command") {
	  print "\nCannot find utility : $1 . Skipping utilities test \n";
	  if (!$opt_p && !$opt_m) { system ("rm -rf $workdir"); }
	  exit $failures;
	}
      }

      elsif ( $_ =~ /^Not_Util/) {
        $command = $octopus_exe;
      }

      elsif ( $_ =~ /^Processors\s*:\s*(.*)\s*$/) {
	$np = $1;
      }

      elsif ( $_ =~ /^Input\s*:\s*(.*)\s*$/) {

        $input_base = $1;
        $input_file = dirname($opt_f) . "/" . $input_base;
      
	if ( $opt_m ) {
	    print "\n\nFor input file : $input_file\n\n";
	    $return_value = 0;
	    $matchdir = "$workdir/$input_base";
	} else {
          if( -f $input_file ) {
            print "\n\nUsing input file : $input_file \n";
            system("cp $input_file $workdir/inp");
            # Ensure that the input file is writable so that it can
            # be overwritten by the next test.
            $mode = (stat "$workdir/inp")[2];
            chmod $mode|S_IWUSR, "$workdir/inp";
          } else {
            die255("ERROR: could not find input file: $input_file\n");
          }
      
	  print "\nStarting test run ...\n";

	  $command_suffix = $command;

	  # serial or MPI run?
	  if ( $octopus_exe =~ /mpi$/) {
            if("$global_np" ne "") {
		$np = $global_np;
            }
	    # utility runs should be on only one processor
	    if ( $command_suffix !~ /mpi$/) {
		$np = 1;
	    }
	    # we do not need to care if mpiexec works if we are just doing a dry run
	    if( -x "$mpiexec_raw" || $opt_n) {
	      if ("$mpiexec" =~ /ibrun/) { # used by SGE parallel environment
		  $specify_np = "";
		  $my_nslots = "MY_NSLOTS=$np";
	      } elsif ("$mpiexec" =~ /runjob/) { # used by BlueGene
		  $specify_np = "--np $np --exe";
		  $my_nslots = "";
	      } elsif ("$mpiexec" =~ /poe/) { # used by IBM PE 
                  $specify_np = ""; 
                  $my_nslots = "MP_PROCS=$np"; 
	      } else { # for mpirun and Cray's aprun
		  $specify_np = "-n $np";
		  $my_nslots = "";
	      }
	      $command_line = "cd $workdir; $my_nslots $mpiexec $specify_np $machinelist $aexec $command_suffix > out";
	    } else {
	      print "No mpiexec found: Skipping parallel test \n";
	      if (!$opt_p && !$opt_m) { system ("rm -rf $workdir"); }
	      exit 254;
	    }
	  } else {
	      $command_line = "cd $workdir; $aexec $command_suffix > out ";
	  }

# MPI implementations generally permit using more tasks than actual cores, and running tests this way makes it likely for developers to find race conditions.
	  if($np > 4) {
	      print "Note: this run calls for more than the standard maximum of 4 MPI tasks.\n";
	  }

	  print "Executing: " . $command_line . "\n";

	  if ( !$opt_n ) {
	    $test_start = [gettimeofday];
	    $return_value = system("$command_line");
	    $test_end   = [gettimeofday];

	    $elapsed = tv_interval($test_start, $test_end);
	    printf("\tElapsed time: %8.1f s\n\n", $elapsed);

	    if($return_value == 0) {
	      print "Finished test run.\n\n";
	      printf "%-40s%s", " Execution", ": \t [ $color_start{green}  OK  $color_end{green} ] \n";
	      
	    } else {
	      print "\n\nTest run failed with exit code $return_value.\n";
	      print "These are the last lines of output:\n\n";
	      print "----------------------------------------\n";
	      system("tail -20 $workdir/out");
	      print "----------------------------------------\n\n";

	      printf "%-40s%s", " Execution", ": \t [ $color_start{red} FAIL $color_end{red} ] \n\n";

	      $failures++;
	      $test_succeeded = 0;      
	    }
	    $test{"run"} = 1;
	  }

	  # copy all files of this run to archive directory with the name of the
	  # current input file
	  mkdir "$workdir/$input_base";
	  @wfiles = `ls -d $workdir/* | grep -v inp`;
	  $workfiles = join("",@wfiles);
	  $workfiles =~ s/\n/ /g;
	  system("cp -r $workfiles $workdir/inp $workdir/$input_base");
	}
      }

      elsif ( $_ =~ /^Precision\s*:\s*(.*)\s*$/) {
	set_precision($1, $command) ;
      }

      elsif ( $_ =~ /^match/ ) {
	  if (!$opt_n && $return_value == 0) {
	      if(run_match_new($_)){
		  printf "%-40s%s", "$name", ":\t [ $color_start{green}  OK  $color_end{green} ] \t (Calculated value = $value) \n";
		  if ($opt_v) { print_hline(); }
	      } else {
		  printf "%-40s%s", "$name", ":\t [ $color_start{red} FAIL $color_end{red} ] \n";
		  print_hline();
		  $test_succeeded = 0;
		  $failures++;
	      }
	  }
      } else {
	  die255("ERROR: Unknown command '$_'\n");
      }
    }

  }

  if ($opt_l && !$opt_m && !$opt_n)  { system ("cat $workdir/out >> out.log"); }
  if (!$opt_p && !$opt_m && $test_succeeded) { system ("rm -rf $workdir"); }

  print "\n";
  close(TESTSUITE)
}

exit $failures;


sub find_executables {
  my $name = "";
  $options = ""; # initialize in case no options specified

  open(TESTSUITE, "<".$opt_f ) or die255("ERROR: cannot open testsuite file '$opt_f'.\n");
  while ($_ = <TESTSUITE>) {

    if ( $_ =~ /^Test\s*:\s*(.*)\s*$/) {
      $name = $1;
    }

    if ( $_ =~ /^Options\s*:\s*(.*)\s*$/) {
      $options = $1;
    }

    if ( $_ =~ /^Programs\s*:\s*(.*)\s*$/) {
      my $i = 0;
      foreach my $program (split(/;/, $1)) {
	$program =  "$program$exec_suffix";
	$program =~ s/^\s+//;
	foreach my $x (@octopus_execs) {
	  $valid = $program cmp $x;
	  if(!$valid) {
	    # check if the executable was compiled with the required options
	    $has_options = 1;
	    foreach my $y (split(/;/, $options)){
	      $command_line = "$exec_directory/$x -c | grep -q $y";
	      $rv = system($command_line);
	      $has_options = $has_options && ($rv == 0)
	    }

	    if($has_options) {
	      $executables[$i] = "$exec_directory/$x";
	      $i = $i+1;
	    }
	  }
	}
      }
    }

  }
  close(TESTSUITE);

  if($name eq "") {
      print STDERR "ERROR: No name was provided with Test tag.\n";
      exit 254;
  }

  # Exit if no suitable executable was found.
  if( @executables == 0 ){
    print STDERR "$color_start{blue} ***** $name ***** $color_end{blue} \n\n";
    print STDERR "$color_start{red}No valid executable$color_end{red} found for $opt_f\n";
    print STDERR "Skipping ... \n\n";
    exit 254;
  }
}

sub run_match_new {
  die255("ERROR: Have to run before matching\n") if !$test{"run"} && !opt_m;

  # parse match line
  my ($line, $match, $pre_command, $ref_value, $off);
  $line = $_[0];
  $line =~ s/\\;/_COLUMN_/g;
  ($match, $name, $pre_command, $ref_value) = split(/;/, $line);
  $pre_command =~ s/_COLUMN_/;/g;
  $ref_value =~ s/^\s*//;
  $ref_value =~ s/\s*$//;

  # parse command
  $pre_command =~ /\s*(\w+)\s*\((.*)\)/;

  my $func = $1;
  my $params = $2;

  # parse parameters
  $params =~ s/\\,/_COMMA_/g;
  my @par = split(/,/, $params);
  for($params=0; $params <= $#par; $params++){
    $par[$params] =~ s/_COMMA_/,/g;
    $par[$params] =~ s/^\s*//;
    $par[$params] =~ s/\s*$//;
  }

  if($func eq "SHELL"){ # function SHELL(shell code)
    check_num_args(1, 1, $#par, $func);
    $pre_command = $par[0];

  }elsif($func eq "LINE") { # function LINE(filename, line, column)
    check_num_args(3, 3, $#par, $func);
    if($par[1] < 0) { # negative number means from end of file
      $line_num = "`wc -l $par[0] | awk '{print \$1}'`";
      $pre_command = "awk -v n=$line_num '(NR==n+$par[1]+1)' $par[0]";
    } else {
      $pre_command = "awk '(NR==$par[1])' $par[0]";
    }
    $pre_command .= " | cut -b $par[2]-";

  }elsif($func eq "LINEFIELD") { # function LINE(filename, line, field)
    check_num_args(3, 3, $#par, $func);
    if($par[1] < 0) { # negative number means from end of file
      $line_num = "`wc -l $par[0] | awk '{print \$1}'`";
      $pre_command = "awk -v n=$line_num '(NR==n+$par[1]+1) {printf \$$par[2]}' $par[0]";
    } else {
      $pre_command = "awk '(NR==$par[1]) {printf \$$par[2]}' $par[0]";
    }

  }elsif($func eq "GREP") { # function GREP(filename, 're', column <, [offset>])
    check_num_args(3, 4, $#par, $func);
    if($#par == 3) {
	$off = $par[3];
    } else {
	$off = 0;
    }
    # -a means even if the file is considered binary due to a stray funny character, it will work
    $pre_command = "grep -a -A$off $par[1] $par[0] | awk '(NR==$off+1)'";
    $pre_command .= " | cut -b $par[2]-";

  }elsif($func eq "GREPFIELD") { # function GREPFIELD(filename, 're', field <, [offset>])
    check_num_args(3, 4, $#par, $func);
    if($#par == 3) {
	$off = $par[3];
    } else {
	$off = 0;
    }
    # -a means even if the file is considered binary due to a stray funny character, it will work
    $pre_command = "grep -a -A$off $par[1] $par[0]";
    $pre_command .= " | awk '(NR==$off+1) {printf \$$par[2]}'";
    # if there are multiple occurrences found by grep, we will only be taking the first one via awk

  }elsif($func eq "SIZE") { # function SIZE(filename)
    check_num_args(1, 1, $#par, $func);
    $pre_command = "ls -lt $par[0] | awk '{printf \$5}'";

  }else{ # error
    printf STDERR "Unknown command '$func'\n";
    return 0;
  }

  # 'set -e; set -o pipefail' (bash 3 only) would make the whole pipe series give an error if any step does;
  # otherwise the error comes only if the last step failed.
  $value = qx(cd $matchdir && $pre_command);
  # Perl gives error code shifted, for some reason.
  $exit_code = $? >> 8;
  if($exit_code) {
      print STDERR "Match command failed: $pre_command\n";
      return 0;
  }

  # extract numeric string (including possibility of NaN)
  if($value =~ /([0-9\-+.eEdDnNaA]+)/) {
      $value = $1;
      chomp $value;
  } else {
      $value = "";
  }

  if(length($value) == 0) {
      print STDERR "Match command returned nothing: $pre_command\n";
      return 0;
  }

  if(!looks_like_number($value)) {
      print STDERR "Match command returned non-numeric value '$value': $pre_command\n";
      return 0;
  }

  # at this point, we know that the command was successful, and returned a number.
  $success = (abs(($value)-($ref_value)) <= $precnum);

  if(!$success || $opt_v) {
    print_hline();
    print "Match".$name.":\n\n";
    print "   Calculated value : ".$value."\n";
    print "   Reference value  : ".$ref_value."\n";
    print "   Difference       : ".abs($ref_value - $value)."\n";
    print "   Tolerance        : ".$precnum."\n\n";
  }

  return $success;
}

sub print_hline {
  print "\n-----------------------------------------\n\n";
}

# return value of environment variable (specified by string argument), or "" if not set
sub get_env {
    if(exists($ENV{$_[0]})) {
	return $ENV{$_[0]};
    } else {
	return "";
    }
}

# args: min num args, max num args, args given, function name
sub check_num_args {
    my $min_num_args   = $_[0];
    my $max_num_args   = $_[1];
    my $given_num_args = $_[2]+1;
    my $func_name      = $_[3];

    if($given_num_args < $min_num_args) {
	die255("$func_name given $given_num_args argument(s) but needs at least $min_num_args.\n");
    }
    if($given_num_args > $max_num_args) {
	die255("$func_name given $given_num_args argument(s) but can take no more than $max_num_args.\n");
    }
}

sub die255 {
    print STDERR $_[0];
    exit 255;
}
