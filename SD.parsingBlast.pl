#!/usr/bin/perl -w


print "Reading input fasta: frg400k.fasta ...\n";
my %fadb;
my $gene;
open(DB, "frg400k.fasta") or die"";
while(<DB>){
	chomp;
	if(/>/){
		$gene = $_;
		$gene =~ s/>//g;
	}else{
		$fadb{$gene} .= $_;
		}
	}
close DB;

print "Reading original position ...\n";
my %posidb; ###record original position and modified position
open(IN, "posidb.txt") or die"";
while(<IN>){
	chomp;
	my @tmp = split(/\s+/,$_);
	next if($tmp[2] eq 'N');
	$posidb{$tmp[0]}->{$tmp[2]} = $tmp[1];
	}
close IN;

print "Identify potential segmental duplication ...\n";
###filter and identify segmental duplication
my %storedb; ###store original position
open(OUT, "> SD.ori_posi.tmp") or die""; ###output SD original position
open(IN, "self.aln.out") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $f1   = $data[0];
	my $a1   = $data[6];
	my $b1   = $data[7];
	my $f2   = $data[1];
	my $a2   = $data[8];
	my $b2   = $data[9];
	my $e    = $data[10];
	next if($e>1e-10);
	my $key1 = $f1.",".$a1.",".$b1;
	my $key2 = $f2.",".$a2.",".$b2;
	next if($key1 eq $key2);
###original position before removing N
	my ($ao,$bo) = sort {$a<=>$b} ($posidb{$f1}->{$a1}, $posidb{$f1}->{$b1});
	my $name = $f1."|o".$ao."_".$bo;
	my $l    = abs($ao-$bo) + 1;
	next if($l<1000);
	print OUT "$f1	$ao	$bo\n";
	foreach my $i ($ao..$bo){
		$storedb{$f1}->{$i}++;
		}
	}
close IN;
close OUT;

print "Extending SD regions ...\n";
%posidb = ();
#####extend SD regions
open(OUT, "> ext.SD.bed") or die"";
foreach my $ctg (keys %storedb){
	my @tmp   = grep {$_} sort {$a<=>$b} keys %{$storedb{$ctg}};
	my $str   = &extend (@tmp);
  my @strdb = split(/,/,$str);
  foreach my $ele(@strdb){
  	$ele    =~ s/\.\./	/g;
  	print OUT "$ctg	$ele\n";
  	}
	}
close OUT;

print "Output potential SD sequences ...\n";
open(OUT, "> SD.potential.fasta") or die"";
open(INN, "ext.SD.bed") or die"";
while(<INN>){
	chomp;
	my ($ctg,$a,$b) = split(/\s+/,$_);
	my $l           = $b - $a + 1;
	my $seq         = substr($fadb{$ctg},$a,$l);
	print OUT ">$ctg|$a..$b\n$seq\n";
	}
close INN;
close OUT;

sub extend {
	my @ext   = @_;
	my $line  = $ext[0];
	for (my $i=1;$i<=$#ext;$i++){
		if(defined($ext[$i+1]) and $ext[$i]+1!=$ext[$i+1]){
			$line  .= "..".$ext[$i].",".$ext[$i+1]     if($ext[$i+1]+1 == $ext[$i+2]);
			$line  .= "..".$ext[$i].",".$ext[$i+1].",".$ext[$i+2] if($ext[$i+1]+1 != $ext[$i+2]);
			$i++;
		}elsif(!defined($ext[$i+1])){
			$line  .= "..".$ext[$i];
			last;
			}               
		}
	return $line;
	}
