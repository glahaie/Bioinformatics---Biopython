#!/usr/bin/env perl
##
## Usage: mavid.pl seqs.fasta
##
## A script to coordinate the actions of mavid, which include:
## - iterative phylogenetic tree building (using CLUSTALW by default or
##   fastDNAml)
## - progressive alignment (using mavid)
## - iterative refinement  
## Two arguments are required: the name of the sequence file to process
## and the directory where temporary files will be stored.
## 
## Options:
##  -treeiter   Number of iterations on treebuilding step (default 2)
##  -treefile   User-specified guide tree in Newick format
##  -h          Prints this help screen
##  -ml         Uses fastDNAml for tree building (may be slower than CLUSTALW)
##  -r          Uses iterative refinement on the last iteration of the alignment.

#
# Change the directories below to the directories where you have the 
# programs installed.
#

my $MAVID_DIRECTORY = "";
my $CLUSTALW_DIRECTORY = "";
my $FASTDNAML_DIRECTORY = "";


use warnings;
use strict;
use English qw( -no_match_vars );
use Getopt::Long;
use FileHandle;
use File::Basename;
use File::stat;
use Time::localtime;
use Cwd;


## Locations of the required executable files - make sure they're in your path
my %exe = (
  fastDNAml => "../utils/fasta2phylip/fasta2phylip",
  clustalw  => "clustalw2",
  mavid     => "./mavid",
  randtree  => "../utils/randtree/randtree",
  root_tree      => "../utils/root_tree/root_tree"
);
warn "$0: WARNING: can't run all required executables\n" if(checkExes(\%exe));


## The main options structure - options used in this program.
my %option = (
  h          => 0,
  ml         => 0,
  r          => 0,
  treeiter   => 2,
  treefile   => "",
);
GetOptions(
  "h"            => \$option{h},
  "ml"           => \$option{ml},
  "r"            => \$option{r},
  "treeiter=i"   => \$option{treeiter},
  "treefile=s"   => \$option{treefile},
) || die("Problem processing command-line options: $!\n");


usage(0) if($option{h});
usage(1) if(@ARGV != 1);


my $seq = shift(@ARGV);
my $dir = dirname($seq); $dir =~ s@\/*$@\/@;
my $mfafile = "${dir}mavid.mfa";
my $phyfile = "${dir}mavid.phy";
my $tree = "${dir}mavid.ph";

if( !$option{treefile} ){
    runrandtree($seq,\%option,\%exe);
    rename( "treefile", $tree );
}
else {
    `cp $option{treefile} $tree`;
}

for(my $i=0;$i < $option{treeiter}; $i++)
{
  run_mavid($tree,$seq,\%option,\%exe);
  if( $option{ml} ){
      runfastDNAml($phyfile,$tree,\%option,\%exe);
  }
  else {
      runclustalw($mfafile,$tree,\%option,\%exe);
  }
  runroot_tree($tree,\%option,\%exe);
}
exit(0);

#
# Usage function - will print out commented lines at the top of this 
# file until it encounters a non-commented line.  (Skips the first.)
#
sub usage {
  my ($die) = @_;
  open(MYSELF,$0) || die "$0: can't open $0 to print usage!\n";
  while(<MYSELF>) {
    if(/^\#!/) {
      next;
    } elsif (/^\s*\#(.*)$/) {
      eval "print STDERR \"$1\n\"";
    } else {
      last;
    }
  }
  close MYSELF;
  die if $die;
}

##
## Verify we have execute authority on all files specified
## in the supplied hash (passed by reference)
##
sub checkExes {
  my($exe) = @_;
  my($problem) = 0;

  foreach my $type (keys %$exe) {
    my $result = `bash -c "type -p $exe->{$type}"`;
    chomp $result;
    if($result =~ /^\s*$/) {
      warn "Can't execute ",$exe->{$type},"\n";
      $problem = 1;
    }
  }
  
  return($problem);
}


# Run the mavid executable (assuming that the tree file has been created)
sub run_mavid {
  my($tree,$seq,$opt,$exe) = @_;
  my $command = $exe->{mavid};
  if( $option{r} ){
      $command .= " -r";
  }
  $command .= " $tree $seq";
  print "$command\n";
  warn("WARNING: non-zero return from mavid.  Call was:\n$command\n") if(system($command));
}


# Run the fastDNAml executable (assuming that the tree file has been created)
sub runfastDNAml {
  my($phyfile,$tree,$opt,$exe) = @_;
  my $command = $exe->{fastDNAml} . " < $phyfile";
  print "$command\n";
  my ($pid) = (`$command | grep "Tree also written to"` =~ /treefile\.(.*)$/);
  if( !$pid ){
      die "Error while running fastDNAml. Aborting.\n";
  }
  `mv treefile.$pid $tree`;
  `rm checkpoint.$pid`;
}


# Run the clustalw executable (assuming that the tree file has been created)
sub runclustalw {
  my($mfafile,$tree,$opt,$exe) = @_;
  my $command = $exe->{clustalw} . " $mfafile -tree";
  print "$command\n";
  warn("WARNING: non-zero return from clustalw.  Call was:\n$command\n") if(system($command));
  my $ofile = $mfafile;
  $ofile =~ s/\.[^.]*$//;
  $ofile .= ".ph";
  if( !($ofile eq $tree) ){
      `mv $ofile $tree`;
  }
}


# Run the randtree executable
sub runrandtree {
  my($seq,$opt,$exe) = @_;
  my $output;
  my $command = $exe->{randtree} . " $seq";
  print "$command\n";
  $output = `$command`;
  if( $output !~ /;/ ){
      die "Error while running randtree. Aborting.\n";
  }
  open OUT, ">$tree";
  print OUT $output;
  close OUT;
}


sub runroot_tree {
  my($tree,$opt,$exe) = @_;
  my ($command,$output);

  $command = $exe->{root_tree} . " $tree";
  print "$command\n";
  $output = `$command`;
  if( $output !~ /;/ ){
      die "Error while running root_tree. Aborting.\n";
  }
  open OUT, ">$tree";
  print OUT $output;
  close OUT;
}
