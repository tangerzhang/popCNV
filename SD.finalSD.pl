#!/usr/bin/perl -w

my %countdb;
open(IN, "global.align.out") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	$countdb{$data[0]}++;
	}
close IN;

open(OUT, "> final.SD.bed") or die"";
open(INN, "ext.SD.bed") or die"";
while(<INN>){
	chomp;
	my ($name,$a,$b) = split(/\s+/,$_);
	my $key          = $name."|".$a."..".$b;
	next if(!exists($countdb{$key}) or $countdb{$key}<2);
	my @tmp          = split(/\_/,$name);
	my $ctg          = $tmp[0];
	my $sa           = $tmp[1];
	my $ctg_a        = $sa + $a - 1;
	my $ctg_b        = $sa + $b - 1;
	my $l            = abs($ctg_a-$ctg_b) + 1;
	next if($l<1000);
	print OUT "$ctg	$ctg_a	$ctg_b\n";
	}
close INN;
close OUT;


