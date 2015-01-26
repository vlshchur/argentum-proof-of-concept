#!/usr/bin/perl

use strict;
use warnings;

my($fn, $lines);
my @m;
my $tmp;
my $sub;
my $mode;
my $ll = -1;
my $counter=0;

print "Hello! Here are some instruments to edit MACS file.\n";
if (!defined $ARGV[0]){
	die("Script called without file name.\n");
}
$fn = $ARGV[0];
if (!defined $ARGV[1] && $ARGV[1] != "-t" && $ARGV[1] != "-d" && $ARGV[1] != "-s"){
	die("Choose mode:\n-t to leave only topology\n-d to delete trees\n-s to show symbol structure (replace tabs by TAB and spaces by SPA)\n-bin to leave only binary.\n");
}
$mode = $ARGV[1];

if (defined $ARGV[2]){
	$ll = $ARGV[2];
}
$mode = $ARGV[1];

open INP, "<$fn" or die $!;
open OUT, ">out-macs.txt" or die $!;

if ($mode eq "-t"){
	while (<INP>){
		$tmp = $_;
		$tmp =~ s/\:[0-9]*\.[0-9]*,/,/g;
		$tmp =~ s/:[0-9]*\.[0-9]*\)/\)/g;
		print OUT $tmp;
	}
}
if ($mode eq "-d"){
	while (<INP>){
		$tmp = $_;
		if ($tmp !~ /^NEWICK_TREE/){
			print OUT $tmp;
		}
	}
}
if ($mode eq "-bin"){
	print OUT "#BINARY file created by macs-edit.pl script from $fn";
	while (<INP>){
		$tmp = $_;
		if ($tmp !~ /^SITE/){
			next;
		}
		@m = $tmp =~ /([01]*)$/;
#		print join(", ", @m) . "\n";
		print OUT "\n" . $1;
		$counter++;
		if ($counter == $ll){
			last;
		}
	}
}
if ($mode eq "-s"){
	while (<INP>){
		$tmp = $_;
		$tmp =~ s/\t/ TAB /g;
		$tmp =~ s/ / SPA /g;
		print $tmp;
	}
}
print "\n\n";
close INP;
close OUT;

