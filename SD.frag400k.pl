#!/usr/bin/perl -w

die "Usage: perl $0 ref.fasta\n" if(!defined $ARGV[0]);

my $refSeq = $ARGV[0];

open(OUT, "> frg400k.fasta") or die"";
open(IN, $refSeq) or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($name,$seq) = split(/\n/,$_,2);
	$name           =~ s/\s+.*//g;
	$seq            =~ s/\s+//g;
	my $str         = uc $seq;
	my $len         = length $seq;
	my $gapN        = ($str =~ tr/N//);
#	next if($len - $gapN < 400000);
	my ($a,$b,$l,$frgseq,$frgname);
	for(my $i=1;$i<=$len;$i=$i+400000){
		$a            = $i;
		$t            = $a + 400000 - 1;
		$b            = ($t<$len)?$t:$len;
		$l            = abs($b - $a) + 1;
		$frgseq       = substr($seq,$a-1,$l);
		$frgname      = $name."_".$a."_".$b;
		print OUT ">$frgname\n$frgseq\n";
		}
	}
close IN;
close OUT;

open(OUT, "> frg400k.rr.fasta") or die"";
open(TMP, "> posidb.txt") or die"";
open(IN, "frg400k.fasta") or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($frg,$seq) = split(/\n/,$_,2);
	my @strdb      = split('',$seq);
	my $mp         = 1;
	foreach my $i(0..$#strdb){
		my $op       = $i+1; ###original position
		my $mp       = ($strdb[$i] eq "N")?"N":$mp++;
		print TMP "$frg	$op	$mp\n";
		}
	my $noN_seq   = $seq;
	   $noN_seq   =~ s/N//g;
	my $len_noN   = length $noN_seq;
	next if($len_noN<100);
	print OUT ">$frg\n$noN_seq\n";
	}
close IN;
close OUT;
close TMP;

