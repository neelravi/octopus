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

Report bugs to <appel\@physik.fu-berlin.de>.
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
#  $cvs_id = qw($Id$);
  $cvs_id = "\$Id: oct-run_regression_test.pl 2423 2006-09-24 21:25:52Z acastro \$";

  open(TEMPLATE, ">".$opt_c );

  print TEMPLATE <<EndOfTemplate;
# -*- coding: utf-8 mode: shell-script -*-
# $cvs_id

Test       : $opt_c
Author     : $author
Date       : $date
Arch       : $arch
Release    : 2.2.0pre1
Programs   : octopus
TestGroups : short-run
Enabled    : Yes

Input: 01-template.01-ground_state.inp

# add your own matches
#
# Examples (of course you have to uncomment the lines :)
# Match ; TotalEnergy ; grep -A 6 '^Energy:' static/info ; \\s*Total\\s*=\\s*-0.33802679
# Match ; Eigenvalues ; grep -A 2 'Eigenvalues \\[' static/info | tail -1 ; \\s*-0.13735\\d*\\s*2.0000\\d*
EndOfTemplate

  close(TEMPLATE);
  print "Template written to: $opt_c \n";
  exit 0;
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
# default number of processors for MPI runs is 2
$np = 2;


# Figure out which are the executables to test
open(TESTSUITE, "<".$opt_f ) or die "cannot open testsuite file\n";
my @executables;
while ($_ = <TESTSUITE>) {

 if ( $_ =~ /^Test\s*:\s*(.*)\s*$/) {
  $test{"name"} = $1;
 }

 if ( $_ =~ /^Programs\s*:\s*(.*)\s*$/) {
  $i = 0;
  foreach my $program (split(/;/,$1)) {
   $program =  "$program$exec_suffix";
   $program =~ s/^\s+//;
   foreach my $x (@octopus_execs) {
    $valid = $program cmp $x;
    if(!$valid) {
     $executables[$i] = "$exec_directory/$x";
     $i = $i+1;
    };
   }
  } 
 }
}
close(TESTSUITE);

# Die if not suitable executable was found.
if( @executables == 0 ){
 print "$color_start{blue} ***** $test{\"name\"} ***** $color_end{blue} \n\n";
 print "$color_start{red}No valid executable$color_end{red} found for $opt_f\n";
 print "Skipping ... \n\n";
 exit 255;
}


# This variable counts the number of failed testcases.
$failures = 0;

# Loop over all the executables.
foreach my $octopus_exe (@executables){

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
 print SCRIPT "rm -rf tmp static status *_tmp *_static out.oct out ds* td.* \n";
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
	   if( -x "$mpirun") {
	     print "Executing: cd $workdir; $mpirun -np $np $octopus_exe_suffix > out 2>&1 \n";
	     system("cd $workdir; $mpirun -np $np $octopus_exe_suffix > out 2>&1");
	   } else {
	     print "No mpirun found: Skipping parallel test \n";
	     exit 255;
	   }
	 } else {
	   print "Executing: cd $workdir; $octopus_exe_suffix > out 2>&1 \n";
	   system("cd $workdir; $octopus_exe_suffix > out 2>&1");
	 }
	 system("grep -B2 -A5 'Running octopus' $workdir/out > $workdir/build-stamp");
	 print "Finished test run.\n\n"; }
       else {
	 if(!$opt_i) { print "cd $workdir; $octopus_exe_suffix < inp > out 2>&1 \n"; }
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

   if ( $_ =~ /^Match/ && !$opt_n) {
     $- = s/\\;/_COLUMN_/g;
     ($match, $name, $pre_command, $regexp) = split(/;/, $_);
     $name        =~ s/_COLUMN_/;/g;
     $pre_command =~ s/_COLUMN_/;/g;
     $regexp      =~ s/_COLUMN_/;/g;

     if(!opt_m) {
       die "have to run before matching" if !$test{"run"};
     }

     if ($opt_v) { print "$pre_command \n"; }
     if ($opt_v) { print "$regexp \n"; }
     $regexp =~ s/^\s*//;
     $regexp =~ s/\s*$//;
     $lineout = `cd $workdir; $pre_command`;

     # append the command and the regexp also to the shell script matches.sh in the
     # current archive directory
     open(SCRIPT, ">>$mscript");
     print SCRIPT "echo ", "="x60, "[ $name - pre command ] \n";
     print SCRIPT "$pre_command\n";
     print SCRIPT "echo ", "-"x60, "[ $name - regular expression ] \n";
     print SCRIPT "echo $regexp\n";
     print SCRIPT "echo;echo\n";
     close(SCRIPT);

     if ( $lineout =~ /$regexp/ ) {
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

 close(TESTSUITE)

}


if (!$opt_i) { print "\n\n\n"; }

exit $failures;
