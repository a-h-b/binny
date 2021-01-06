#!/usr/bin/env perl

use strict;
use Getopt::Long;

my ($hmmscan_result_file,$gf_file,$namearg,$out_file,$cutoff,$geneno,$keepAll,$trimlast,$trimall);
GetOptions('i|hmm=s' => \$hmmscan_result_file,
	   'a|annotationgff=s' => \$gf_file,
	   'n|name=s' => \$namearg,
	   'o|out=s' => \$out_file,
	   's|bitscore=f' => \$cutoff,
	   'g|gennumber=i' => \$geneno,
	   'k|keepall' => \$keepAll,
	   'l|trimlast' => \$trimlast,
	   't|trimall' => \$trimall,);

if ($gf_file eq $out_file) {
	die "Input and ouput file should not be the same!\n";
}
if (not defined $namearg) {
	die "a name for the field must be given using the -n argument!\n"
}

if (defined $geneno) {
	$cutoff = &log2($geneno);
}elsif (not defined $cutoff){
	$cutoff = 0;
}


my %hit_hash=();
my @all = ();
my $line = ();

### ESS

print STDOUT "# hmm result: $hmmscan_result_file\n";

close(FH);
open(FH, "$hmmscan_result_file") or die "Failed to open $hmmscan_result_file\n";
@all = <FH>;	
shift @all; # to get rid of header
while (@all) {
	$line = shift @all;
	chomp($line);
	next if $line=~/^#/;
	my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,@desc_of_target) = split(/ +/, $line);
	#print STDERR "query_name : $query_name line: $line\n";		
	#push @{$ess_hmmscan_result_hash{$hit_name}}, $line;
	if ($seqscore >= $cutoff) {
		if ($trimlast) {
			$query_name =~ s/(.+)_.+/$1/; #s/_.+//g;
		}elsif ($trimall) {
			$query_name =~ s/_.+//g;
		}
		if (!$hit_hash{$hit_name}) {
			$hit_hash{$hit_name} = [$query_name,$seqscore];
		}else{
			if ($seqscore > $hit_hash{$hit_name}[1] && !$keepAll) {
				$hit_hash{$hit_name} = [$query_name,$seqscore];
			}elsif ($seqscore == $hit_hash{$hit_name}[1] || $keepAll){
				my $tmp_qn = $hit_hash{$hit_name}[0].";".$query_name;
				$hit_hash{$hit_name} = [$tmp_qn,$seqscore];
			}
		}
	}
} 
close(FH);

### GFF
print STDOUT "# gff: $gf_file\n";
my $outfile=$out_file;
open(FH, "$gf_file") or die "Failed to open $gf_file\n";
open RESULT,">$outfile";
@all = <FH>;	
while (@all) {
	$line = shift @all;
	chomp($line);
	next if $line=~/^#/;
	my ($contig, $tool, $type, $start,$end,$point,$sense,$zero,$description) = split(/\t/, $line);
	my @descriptors = split(/;/, $description);
	my $product = $descriptors[0];
	$product =~ s/ID=//;
	if ($hit_hash{$product}) {
		$line =~ s/;$//;
		my @hmmgenes = split(/;/, $hit_hash{$product}[0]);
		my $addarg = "";
		foreach my $hmmgene (@hmmgenes){
			$addarg .= ";".$namearg."="."$hmmgene"
		}
		print RESULT "$line".$addarg."\n";
	}#else{
	#	print RESULT "$line\n";
	#}
} 
close(FH);
close(RESULT);

%hit_hash = ();


sub log2 {
	my $n = shift;
        return log($n)/log(2);
}

