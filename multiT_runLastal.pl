#!/usr/bin/perl -w

use threads;
use threads::shared; 

die "Usage: perl $0 query.fasta database.fasta threads identity out.aln\n" 
if(!defined($ARGV[0]) || !defined($ARGV[1]) || !defined($ARGV[2]) || !defined($ARGV[4]));

my $query      = $ARGV[0];
my $dataB      = $ARGV[1];
my $maxThreads = int $ARGV[2]/2;
my $identity   = (defined $ARGV[3])?$ARGV[3]:60;
my $outaln     = $ARGV[4];
system("rm -rf faByChr");
system("mkdir faByChr");
system("rm -rf tmp");
system("mkdir tmp");
#system("rm dblast*");
#system("lastdb -uNEAR -R01 dblast $dataB");

my %infordb;
my %namedb;
my $count = 0;
open(IN, $query) or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($name,$seq) = split(/\n/,$_,2);
	$name  =~ s/\s+.*//g;
	$seq   =~ s/\s+//g;
	$count++;
	$infordb{$count}->{'name'}  = $name;
	$infordb{$count}->{'seq'}   = $seq;
	my $rname       = $count;	
  $namedb{$rname} = $name;
  $infordb{$count}->{'rname'} = $rname;
  my $outfile     = $rname.".fasta";
	open(my $out, ">faByChr/$outfile") or die"";
	print $out ">$rname\n$seq\n";
	close $out;	

	}
close IN;

my $num_of_job = keys %infordb;
my $i=1;
while($i<=$num_of_job){
	while(scalar(threads->list())<=$maxThreads){
		my $scaf_n = $infordb{$i}->{'rname'};
		  next if(!defined($scaf_n));
		threads->new(\&run_lastal,$scaf_n);
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

open(OUT, "> $outaln") or die"";
open(IN, "cat tmp/*.out |") or die"";
my $content = <IN>;
@linedb = split(/\n/,$content);
foreach my $line (@linedb){
	my @data = split(/\s+/,$line);
	$data[0] = $namedb{$data[0]} if(exists($namedb{$data[0]}));
	foreach my $i(0..$#data){
		print OUT "$data[$i]	";
		}
	print OUT "\n";
	}
close IN;
close OUT;
system("rm -rf tmp/");

sub run_lastal{
  my $scaf = shift;
  my $cmd  = "lastal -P 2 dblast faByChr/".$scaf.".fasta -a 180 -b 1 -q 80 -r 30 |maf-convert blasttab |awk '\$3>$identity' > tmp/".$scaf.".out";
  system($cmd);	
	}
