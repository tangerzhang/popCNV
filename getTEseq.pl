#!/usr/bin/perl -w

die "perl $0 ref.fasta repeatmask.out \n" if(!defined($ARGV[0]) or !defined($ARGV[1]));

my %refdb;
my $scaf;
open(IN, $ARGV[0]) or die"";
while(<IN>){
	chomp;
	if(/>/){
		$scaf = $_;
		$scaf =~ s/\s+.*//g;
		$scaf =~ s/>//g;
	}else{
		$refdb{$scaf} .= $_;
		}
	}
close IN;

foreach my $scaf(sort keys %refdb){
	$refdb{$scaf} =~ s/\s+//g;
	}

open(OUT, ">TEseq.fasta") or die"";
open(IN, $ARGV[1]) or die"";
while(<IN>){
	chomp;
	next if(/position/);
	next if(/score/);
	$_ =~ s/^\s+//g;
	my @data = split(/\s+/,$_);
	next if(@data<2);
	my $name = $data[4];
	my $a    = $data[5];
	my $b    = $data[6];
	my $l    = abs($b-$a)+1;
	my $seq  = substr($refdb{$name},$a,$l);
	my $newn = $name."|".$a."|".$b;
	print OUT ">$newn\n$seq\n";
	}
close IN;
close OUT;


