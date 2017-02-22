#!/usr/bin/perl -w

use threads;
use threads::shared; 
use Getopt::Std;
getopts "i:r:t:s:";


if ((!defined $opt_i)|| (!defined $opt_r)) {
    die "************************************************************************
    Usage: thorough_mask.pl -i ref.fasta -r TEseq.fasta -t cpu -s SD_interval.txt
      -h : help and usage.
      -i : ref.fasta
      -r : TE sequences
      -t : threads, default 12
      -s : segental duplication, optional
************************************************************************\n";
}else{
  print "************************************************************************\n";
  print "Version demo\n";
  print "Copyright to Tanger, tanger.zhang\@gmail.com\n";
  print "RUNNING...\n";
  print "************************************************************************\n";    
        }
my @threads;
my $refSeq  = $opt_i;
my $TEseq   = $opt_r;
my $maxThreads = (defined $opt_t)?$opt_t:12;
my $SDfile  = (defined $opt_s)?$opt_s:"none";

print "1. RepeatMask reference genome with -s \n";
system("ln -s $refSeq ./ref.fasta");
print "RepeatMasker -pa $maxThreads -s -e ncbi -gff -lib $TEseq ref.fasta \n";
system("RepeatMasker -pa $maxThreads -s -e ncbi -gff -lib $TEseq ref.fasta 2&>>log.txt");
print "Done RepeatMask ...\n";

print "2. identify tandom repeats\n";

system("rm -rf faByChr");
system("mkdir faByChr");

my %infordb; 
my $count = 0;
open(IN, "ref.fasta") or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($name,$seq) = split(/\n/,$_,2);
	$name  =~ s/\s+.*//g;
	next if($name eq "" or !defined($name));
	$seq   =~ s/\s+//g;
	my $outfile = $name.".fasta";
	open(my $out, ">faByChr/$outfile") or die"";
	print $out ">$name\n$seq\n";
	close $out;
	$count++;
	$infordb{$count}->{'name'} = $name;
	$infordb{$count}->{'seq'}  = $seq;
	}
close IN;

my $num_of_job = keys %infordb;
my $usedThreads;
$usedThreads   = $num_of_job - 1 if($maxThreads>$num_of_job);
$usedThreads   = $maxThreads     if($maxThreads<$num_of_job);

my $i=1;
while($i<=$num_of_job){
	while(scalar(threads->list())<=$usedThreads){
		my $scaf_n = $infordb{$i}->{'name'};
		next if(!(defined $scaf_n));
		threads->new(\&run_trf,$scaf_n);
		$i++;
		}
  foreach $thread(threads->list(threads::all)){
     if($thread->is_joinable()){
       $thread->join();
     }
   }	
	}

foreach $thread(threads->list(threads::all)){
    $thread->join();
   }

open(OUT, "> tandem.all.bed") or die"";
while(my $dat = glob "*.dat"){
	my $cont = $dat;
	   $cont =~ s/.fasta.1.1.2.80.5.200.2000.dat//g;	
	open(my $fh, $dat) or die"";
	my $content = <$fh>;
	my @linedb  = split(/\n/,$content);
	foreach my $line(@linedb){
		my @data = split(/\s+/,$line);
		next if(@data < 10);
		my $str  = $data[14];
		my $a    = $data[0];
		my $b    = $data[1];
 		my $len  = abs($data[0]-$data[1]);
 		print OUT "$cont	$a	$b	$len	$str\n";
		}
	close $fh;
	}
close OUT;
system("rm *.dat");
print "\nDone tandem duplication ...\n";

print "3. run windowmasker ... \n";
print "windowmasker -mk_counts -in ref.fasta -out mk_counts.txt\n";
print "windowmasker -ustat mk_counts.txt -in ref.fasta -out phase2.out -dust true\n";

system("windowmasker -mk_counts -in ref.fasta -out mk_counts.txt ");
system("windowmasker -ustat mk_counts.txt -in ref.fasta -out phase2.out -dust true ");
print "\nDone window masker ...\n";

print "4. generate unique genome ...\n\n";
my %rpdb;
print "a. reading phase2.out \n\n";
open(IN, "phase2.out") or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($scaf,$info) = split(/\n/,$_,2);
	$scaf            =~ s/\s.*//g;
	my @data         = split(/\n/,$info);
	foreach my $i(0..$#data){
		my ($a,$b) = split(/\s+/,$data[$i]);
		$a = $a - 36; $b = $b + 36;
		$a = 0 if($a<0); 
		foreach my $j($a..$b){
			$rpdb{$scaf}->{$j} += 1;
			}
		}
	}
close IN;

print "b. reading repeatmasker result \n\n";
open(IN, "grep -v '#' ref.fasta.out.gff |") or die"";
$content = <IN>;
@linedb  = split(/\n/,$content);
foreach my $line(@linedb){
	my ($scaf,$a,$b) = (split/\s+/,$line)[0,3,4];
	$a = $a - 36; $b = $b + 36;
	$a = 0 if($a<0); 	
	foreach my $j ($a..$b){
		$rpdb{$scaf}->{$j} += 1;
		}
	}
close IN;

print "c. reading tandem.all.bed \n\n";
open(IN, "tandem.all.bed") or die"";
$content = <IN>;
@linedb  = split(/\n/,$content);
foreach my $line(@linedb){
	my ($scaf,$a,$b) = (split/\s+/,$line)[0,1,2];
	$a = $a - 36; $b = $b + 36;
	$a = 0 if($a<0); 		
	foreach my $j($a..$b){
		$rpdb{$scaf}->{$j} += 1;
		}
	}
close IN;

if($SDfile ne "none"){
	print "d. reading segmental duplication \n\n";
	open(INN, $SDfile) or die"";
  $content = <INN>;
  @linedb  = split(/\n/,$content);
  foreach my $line(@linedb){
  	my ($scaf,$a,$b) = (split/\s+/,$line)[0,1,2];
	  $a = $a - 36; $b = $b + 36;
	  $a = 0 if($a<0); 	  	
  	foreach my $j($a..$b){
  		$rpdb{$scaf}->{$j} += 1;
  		}
  	}
  close INN;	
	
}
my $num_of_N = 0;

open(OUT, "> unique.genome.fasta ") or die"";
foreach my $i(1..$num_of_job){
	my @basedb = split('',$infordb{$i}->{'seq'});
	my $scaf   = $infordb{$i}->{'name'};
	foreach my $j(0..$#basedb){
		$basedb[$j] = "N" if(exists($rpdb{$scaf}->{$j}));
		$num_of_N++       if(exists($rpdb{$scaf}->{$j}));
		next if(!exists($rpdb{$scaf}->{$j}));
		}
	my $str = join('',@basedb);
	print OUT ">$scaf\n$str\n";
	}
close OUT;

system("rm log.txt");

print "Total masked size: $num_of_N\n";

sub run_trf{
 my $scaf = shift;
 my $cmd  = "trf faByChr/".$scaf.".fasta 1 1 2 80 5 200 2000 -d -h ";
 system($cmd);
  }



