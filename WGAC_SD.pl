#!/usr/bin/perl -w

use Getopt::Std;
getopts "i:l:t:d:s:";


if ((!defined $opt_i)|| (!defined $opt_l)) {
    die "************************************************************************
    Usage: WGAC_SD.pl -i ref.fasta -l repeat.lib -t cpu -s stepN
      -h : help and usage.
      -i : ref.fasta
      -l : customed repeat lib, fasta format
      -t : threads, default 12
      -d : debug use only, default -d 1
           clean intermediate files
      -s : int, start from step N, default 1
************************************************************************\n";
}else{
  print "************************************************************************\n";
  print "Version demo\n";
  print "Copyright to Tanger, tanger.zhang\@gmail.com\n";
  print "RUNNING...\n";
  print "************************************************************************\n";    
        }

my $threads      = (defined $opt_t)?$opt_t:12;
my $refSeq       = $opt_i;
my $debug        = (defined $opt_d)?$opt_d:1;
my $step         = (defined $opt_s)?$opt_s:1;
my $fileExist;
##########0. check data ...
print "checking data ...\n";
my $refname = `grep '>' $refSeq`;
die "symbol \"_\" is not allowed in reference header, please remove \n" if ($refname =~ /\_/);

##########1. running RepeatMasker to identify repeats ...

if($step==1){
print "1. running RepeatMasker to identify repeats ...\n";
print "RepeatMasker -pa $threads -lib $opt_l $refSeq\n";
system("RepeatMasker -pa $threads -lib $opt_l $refSeq 2&>log.txt ");
print "Get TE sequences in TEseq.fasta\n";
my $rpout = $refSeq.".out";
print "getTEseq.pl $refSeq $rpout\n";
system("getTEseq.pl $refSeq $rpout");
system("rm $rpout");
print "...Done...\n";
$step++;
}
##########2. split genome to 400-kb trackable fragment and remove repeats
if($step==2){
print "\n2. split genome to 400-kb trackable fragment and remove repeats ...\n";
my $mask_refSeq = $refSeq.".masked";
system("SD.frag400k.pl $mask_refSeq");
print "...Done...\n";
$step++;
}
##########3. self-alignment
if($step==3){
print "\n3. self-alignment ...\n"; 
system("lastdb -uNEAR -R01 dblast frg400k.rr.fasta");
print "multiT_runLastal.pl frg400k.rr.fasta frg400.rr.fasta $threads 88 self.aln.out\n";
system("multiT_runLastal.pl frg400k.rr.fasta frg400.rr.fasta $threads 88 self.aln.out");

print "parsing blast results ...\n";
system("SD.parsingBlast.pl");
print "...filtering fragments ...\n";
print "...self-align finished...\n";
$step++;
}
###########4. global alignment
if($step==4){
print "\n4. global alignment ...\n";
system("lastdb -uNEAR -R01 dblast $refSeq");
print "multiT_runLastal.pl SD.potential.fasta $refSeq $threads 90 global.align.out\n";
system("multiT_runLastal.pl SD.potential.fasta $refSeq $threads 90 global.align.out");
print "parsing global alignment results..\n";
print "Get final Segmental duplication ...\n";
system("SD.finalSD.pl");
print "Number of SD: ";
system("wc -l final.SD.bed");
print "\nTotal Size of SD: ";
system("perl \-e 'while\(\<\>\)\{\$a=(split)[1];\$b=(split)[2];\$sum+=\$b-\$a+1}print \"\$sum\n\"' final.SD.bed");
print "\n";
print "...Global Alignment finished...\n";
$step++;
}
###########clean intermediate files 

#die "All processes finished ...\n" if($debug==0);
#print "\ncleaning intermediate files ...\n" if($debug==1);
#system("rm ext.SD.bed frg400k* SD.ori_posi.tmp SD.potential.fasta posidb.txt dblast*");
#system("rm self.aln.out global.align.out");
#system("rm -rf faByChr");





