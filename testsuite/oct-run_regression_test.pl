#!/usr/bin/perl
#
# Copyright (C) 2005 Heiko Appel
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
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# $Id$

use Getopt::Std;
use File::Basename;
use Fcntl ':mode';


sub usage {

  print <<EndOfUsage;

 Copyright (C) 2005 by Heiko Appel

Usage: oct-run_regression_test.pl [options]

    -n        dry-run
    -v        verbose
    -h        this usage
    -D        name of the directory where to look for the executables   
    -s        exec suffix for the executables
    -c        create template
    -f        filename of testsuite
    -i        print inputfile
    -p        preserve working directories
    -l        copy output log to current directory
    -m        run matches only (assumes there are work directories)

Exit codes:
    0         all tests passed
    255       test skipped
    1..254    number of test failures

Report bugs to <octopus-devel\@tddft.org>
EndOfUsage

  # Option -d is ignored for the moment.
  #    -d        working directory for the tests

  exit 0;
}


sub create_template {
  $date = `date +"%d.%m.%y"`;
  chomp($date);
  $arch = `uname -a`;
  chomp($arch);
  $author = `whoami`;
  chomp($author);
  $author =~ s/^(\w)(.*)/\u$1$2/;
  $cvs_id = "\$Id: oct-run_regression_test.pl 2423 2006-09-24 21:25:52Z acastro \$";

  open(TEMPLATE, ">".$opt_c );

  print TEMPLATE <<EndOfTemplate;
# -*- coding: utf-8 mode: shell-script -*-
# $cvs_id

Test       : $opt_c
Author     : $author
Date       : $date
Arch       : $arch
Release    : 3.1.0pre1
Programs   : octopus
TestGroups : short-run
Enabled    : Yes

Input: 01-template.01-ground_state.inp

# add your own matches
#
# Example (of course you have to uncomment the lines :)
# match ; Total energy       ; GREP(static/info, 'Total       =', 20) ; -146.81378481
# match ; Eigenvalue   [1up] ; GREP(static/info, '1   up', 13) ; -14.466667
# match ; Forces [step  1] ; LINE(td.general/coordinates, -21, 270) ; 1.700891538586e-01

EndOfTemplate

  close(TEMPLATE);
  print "Template written to: $opt_c \n";
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

# Check, if STDOUT is a terminal. If not, not ANSI sequences are
# emitted.
if(-t STDOUT) {
    $color_start{blue}="\033[34m";
    $color_end{blue}="\033[0m";
    $color_start{red}="\033[31m";
    $color_end{red}="\033[0m";
    $color_start{green}="\033[32m";
    $color_end{green}="\033[0m";
}

if (not @ARGV) { usage; }

# Option -d is ignored for the moment.
#getopts("nlvhD:c:f:d:s:ipm");
getopts("nlvhD:c:f:s:ipm");

# Default values
use File::Temp qw/tempdir/;

# Handle options
$opt_h && usage;
$opt_c && create_template;
####NEW This has to be handled in a differenty way.
#if($opt_d)  { $workdir = $opt_d; }
####

my $exec_directory;
if($opt_D) {
 $exec_directory = $opt_D;
 if($exec_directory !~ /^\//){
  $exec_directory = $ENV{PWD}."/$exec_directory";
 }
} else {
 $exec_directory = "/usr/bin";
}

# Find out which executables are available.
opendir(EXEC_DIRECTORY, $exec_directory) || 
 die "Could not open the directory $exec_directory to look for executables";
@octopus_execs = grep { /^oct/ } readdir(EXEC_DIRECTORY);
closedir(EXEC_DIRECTORY);

# determine exec suffix
$exec_suffix = "";
if($opt_s)  { $exec_suffix = $opt_s; }

# MPI stuff
$mpirun = $ENV{MPIRUN};
if ("$mpirun" eq "") { $mpirun = `which mpirun`; }
chomp($mpirun);

# mpirun without arguments (to check if it is available)
$mpirun_raw = $mpirun;
$mpirun_raw =~ s/\ (.*)//;

# default number of processors for MPI runs is 2
$np = 2;


# Figure out which are the executables to test
my @executables;
find_executables();

# This variable counts the number of failed testcases.
$failures = 0;

# Loop over all the executables.
foreach my $octopus_exe (@executables){

  set_precision("default", $octopus_exe);

  $workdir = tempdir('/tmp/octopus.XXXXXX');
  chomp($workdir);

  if (!$opt_m) {
    system ("rm -rf $workdir");
    mkdir $workdir;
  }

  # create script for cleanups in the current workdir
  $mscript = "$workdir/clean.sh";
  open(SCRIPT, ">$mscript") or die "could not create script file\n";
  print SCRIPT "#\!/bin/bash\n\n";
  print SCRIPT "rm -rf tmp static exec *_tmp *_static out.oct out ds* td.* \n";
  close(SCRIPT);
  chmod 0755, $mscript;

  # testsuite
  open(TESTSUITE, "<".$opt_f ) or die "cannot open testsuite file\n";

  while ($_ = <TESTSUITE>) {

    # skip comments
    next if /^#/;

    if ( $_ =~ /^Test\s*:\s*(.*)\s*$/) {
      $test{"name"} = $1;
      if(!$opt_i) {
	print "$color_start{blue} ***** $test{\"name\"} ***** $color_end{blue} \n\n";
	print "Using workdir    : $workdir \n";
	print "Using executable : $octopus_exe\n";
	print "Using test file  : $opt_f \n";
      }
    }

    if ( $_ =~ /^Enabled\s*:\s*(.*)\s*$/) {
      %test = ();
      $enabled = $1;
      $enabled =~ s/^\s*//;
      $enabled =~ s/\s*$//;
      $test{"enabled"} = $enabled;
    }

    # Running this regression test if it is enabled
    if ( $enabled eq "Yes" ) {

      if ( $_ =~ /^Processors\s*:\s*(.*)\s*$/) {
	$np = $1;
      }

      if ( $_ =~ /^Input\s*:\s*(.*)\s*$/) {
	$input_base = $1;
	$input_file = dirname($opt_f) . "/" . $input_base;

	if( -f $input_file ) {
	  print "\n\nUsing input file : $input_file \n";
	  system("cp $input_file $workdir/inp");
	  # Ensure, that the input file is writable so that it can
	  # be overwritten by the next test.
	  $mode = (stat "$workdir/inp")[2];
	  chmod $mode|S_IWUSR, "$workdir/inp";
	} else {
	  die "could not find input file: $input_file\n";
	}


	if ( !$opt_m ) {
	  if ( !$opt_n ) {
	    print "\nStarting test run ...\n";

	    $octopus_exe_suffix = $octopus_exe;

	    # serial or MPI run?
	    if ( $octopus_exe_suffix =~ /mpi$/) {
	      if( -x "$mpirun_raw") {
		print "Executing: cd $workdir; $mpirun -np $np $octopus_exe_suffix > out \n";
		system("cd $workdir; $mpirun -np $np $octopus_exe_suffix > out ");
	      } else {
		print "No mpirun found: Skipping parallel test \n";
		exit 255;
	      }
	    } else {
	      print "Executing: cd $workdir; $octopus_exe_suffix > out \n";
	      system("cd $workdir; $octopus_exe_suffix > out ");
	    }
	    system("sed -n '/Running octopus/{N;N;N;N;N;N;p;}' $workdir/out > $workdir/build-stamp");
	    print "Finished test run.\n\n"; }
	  else {
	    if(!$opt_i) { print "cd $workdir; $octopus_exe_suffix > out \n"; }
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

	# file for shell script with matches
	$mscript = "$workdir/$input_base/matches.sh";
	open(SCRIPT, ">$mscript") or die "could not create script file\n";
	# write skeleton for script
	print SCRIPT "#\!/bin/bash\n\n";
	close(SCRIPT);
	chmod 0755, $mscript;
      }

      if ( $_ =~ /^Precision\s*:\s*(.*)\s*$/) {
	set_precision($1, $octopus_exe) ;
      }

      if ( $_ =~ /^match/ && !$opt_n) {
	if(run_match_new($_)){
	  print "$name: \t [ $color_start{green}  OK  $color_end{green} ] \n";
	  $test_succeded = 1;
	} else {
	  print "$name: \t [ $color_start{red} FAIL $color_end{red} ] \n";
	  $test_succeded = 0;
	  $failures++;
	}
      }

    } else {
      if ( $_ =~ /^RUN/) { print " skipping test\n"; }
    }

  }

  if ($opt_l)  { system ("cat $workdir/out >> out.log"); }
  if (!$opt_p && !$opt_m && $test_succeded) { system ("rm -rf $workdir"); }

  print "\n";
  close(TESTSUITE)
}

exit $failures;


sub find_executables(){
  my $name;

  open(TESTSUITE, "<".$opt_f ) or die "cannot open testsuite file\n";
  while ($_ = <TESTSUITE>) {

    if ( $_ =~ /^Test\s*:\s*(.*)\s*$/) {
      $name = $1;
    }

    if ( $_ =~ /^Programs\s*:\s*(.*)\s*$/) {
      my $i = 0;
      foreach my $program (split(/;/, $1)) {
	$program =  "$program$exec_suffix";
	$program =~ s/^\s+//;
	foreach my $x (@octopus_execs) {
	  $valid = $program cmp $x;
	  if(!$valid) {
	    $executables[$i] = "$exec_directory/$x";
	    $i = $i+1;
	  }
	}
      }
    }

  }
  close(TESTSUITE);

  # Die if not suitable executable was found.
  if( @executables == 0 ){
    print stderr "$color_start{blue} ***** $test{\"name\"} ***** $color_end{blue} \n\n";
    print stderr "$color_start{red}No valid executable$color_end{red} found for $opt_f\n";
    print stderr "Skipping ... \n\n";
    exit 255;
  }
}

sub run_match_new(){
  die "Have to run before matching" if !$test{"run"} && !opt_m;

  # parse match line
  my $line, $match, $name, $pre_command, $ref_value;
  $line = @_[0];
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

  if    ($func eq "SHELL"){ # function SHELL(shell code)
    $pre_command = $par[0];

  }elsif($func eq "LINE") { # function LINE(filename, line, column)
    $pre_command = "cat $par[0]";
    if($par[1] > 0){
      $pre_command .= " | tail -n +$par[1] | head -n 1";
    } else {
      $pre_command .= " | tail -n  $par[1] | head -n 1";
    }
    $pre_command .= " | cut -b $par[2]- | perl -ne '/\\s*([0-9\\-+.eEdD]*)/; print \$1'";

  }elsif($func eq "GREP") { # function GREP(filename, 're', column <, offset>)
    my $off = 1*$par[3];
    $pre_command = "n=\$(cat -n $par[0] | grep $par[1] | head -n 1 | awk '{print \$1;}');";
    $pre_command .= " cat $par[0] | tail -n +\$((n+$off)) | head -n 1";
    $pre_command .= " | cut -b $par[2]- | perl -ne '/\\s*([0-9\\-+.eEdD]*)/; print \$1'";

  }else{ # error
    printf stderr "Unknown command '$func'\n";
    return 0;
  }

  printf("%s\n", $pre_command) if ($opt_v);

  my $value = `cd $workdir; $pre_command`;

  # append the command and the regexp also to the shell script matches.sh in the
  # current archive directory
  open(SCRIPT, ">>$mscript");
  print SCRIPT "
echo '", "="x4, " [ $name - pre command ]'
$pre_command
echo
echo '", "-"x4, " [ $name - ref value   ]'
echo $ref_value
export LINE=`$pre_command`
perl -e 'print \"Match: \".(abs(\$ENV{LINE}-($ref_value)) <= $precnum ? \"OK\" : \"FAILED\");'
echo
echo";
  close(SCRIPT);

  return (abs(($value)-($ref_value)) <= $precnum);
}

