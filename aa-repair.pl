#!/usr/bin/perl

use strict;
use warnings;

my $counter=0;
my $counter1=0;
my($fn, $lines);
my @m;
my $tmp;
my $tmp1;
my $sub;
my $repStr=0;

print "Hello! Let us repair ancestral alleles.\n";
if (!defined $ARGV[0]){
	die("Script called without file name.");
}
if (defined $ARGV[1]){
	$lines = $ARGV[1];
}
else{
	$lines = -1;
}
$fn = $ARGV[0];
open INP, "<$fn" or die $!;
open OUT, ">out.vcf" or die $!;
open OUTBIN, ">out-bin.txt" or die $!;
print OUTBIN "#BINARY file created by aa-repair.pl script from $fn\n";

while (<INP>){
	$tmp = $_;
	@m = $tmp =~ m/([ATCG])(.)([ATCG])(.*AA=)([ATCG])/;
#	if ($counter == 19875){
#		print "$counter1 !!!\n";
#	}
	$counter1++;
	if ($tmp !~ /^#/ && !@m){
		next;
	}
	if ($tmp =~ /^#/){
		print OUT $tmp;
		next;
	}
#	if ($tmp =~ /\|[234]/ || $tmp =~ /[234]\|/){
#		next;
#	}
	#print "$m[0] $m[1] $m[2]\n";
#	if ($m[0] ne $m[4] && $m[2] ne $m[4] ){
#		next;
#	}
	if ($m[2] eq $m[4]){
#		next;
		#$tmp =~ s/$m[0]$m[1]$m[2]$m[3]$m[4]/$m[2]$m[1]$m[0]$m[3]$m[4]/;
		$repStr++;
		$tmp =~ s/0\|0/2\|2/g;
		$tmp =~ s/1\|1/0\|0/g;
		$tmp =~ s/2\|2/1\|1/g;
		$tmp =~ s/0\|1/2\|2/g;
		$tmp =~ s/1\|0/0\|1/g;
		$tmp =~ s/2\|2/1\|0/g;
		if ($repStr%5000 == 1){
			print "$repStr strings already repaired.\nline number = $counter.\n";
		}
	}
#	print OUT $tmp;
	@m = $tmp =~ m/[01]\|[01]/g;
	$tmp1 = "";
	foreach (@m){
		$_ =~ tr/\|//d;
		$tmp1 = $tmp1.$_;
	}
	print OUTBIN $tmp1."\n";
	$counter ++;
	if ($counter == $lines){
		last;
	}
}
print "\n$repStr strings repaired.\n";
close INP;
close OUT;
close OUTBIN;

