#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Fasta;
use Getopt::Long;

my ($fastaFile,$gffFile,$outFile,$cutFile);
my $cutoff = 1000;
GetOptions('f|fasta=s' => \$fastaFile,
	   'g|gff=s' => \$gffFile,
	   'o|output=s' => \$outFile,
	   'l|log=s' => \$cutFile,
	   'c|cutoff=i' => \$cutoff);

my %rRNAs = ();

open (IN, $gffFile) or die "Failed to open $gffFile\n";
while (my $line = <IN>){
    chomp $line;
	my @att = split("\t",$line);
	if ($att[2] eq "rRNA"){
		push @{$rRNAs{$att[0]}[0]}, $att[3];
		push @{$rRNAs{$att[0]}[1]}, $att[4];
	}
}
close(IN);

open(OUT, '>', $outFile) or die "could not open file '$outFile' \n";
open(CUT, '>', $cutFile) or die "could not open file '$cutFile' \n";
my $db = Bio::DB::Fasta->new( $fastaFile );
my @ids = $db->get_all_ids;
foreach my $contig (@ids){
    my $sequence = $db->seq($contig);
    if  (!defined( $sequence )) {
        print STDERR "Sequence $contig not found. \n";
        next;
    }
    my $lenSeq = length $sequence;
    unless ($lenSeq < $cutoff) {
		if (exists $rRNAs{$contig}) {
			my @useStarts = ();
			my @useEnds = ();
			if ($#{$rRNAs{$contig}[0]} > 0) {
				my @starts = @{$rRNAs{$contig}[0]};
				my @ends = @{$rRNAs{$contig}[1]};
				my @idx = sort { $starts[$a]+0 <=> 0+$starts[$b] } 0 .. $#starts;
				@starts = @starts[@idx];
				@ends = @ends[@idx];
				my $currentStart = $starts[0];
				my $currentEnd = $ends[0];
				for (my $i = 0; $i <= $#starts; $i++) { 
					if ($starts[$i] <= $currentEnd){
						$currentEnd = $ends[$i] if ($ends[$i] >= $currentEnd)
					} else {
						unshift @useStarts, $currentStart;
						unshift @useEnds, $currentEnd;
						$currentStart = $starts[$i];
						$currentEnd = $ends[$i];
					}
				}
				unshift @useStarts, $currentStart;
				unshift @useEnds, $currentEnd;
			} else {
				push @useStarts, $rRNAs{$contig}[0][0];
				push @useEnds, $rRNAs{$contig}[1][0];
			}
			my $lenCutSeq = length $sequence;
			for (my $j = 0 ; $j <= $#useStarts ; $j++ ) {
				my $startcut = $useStarts[$j] -1; 
				my $length = $useEnds[$j] - $startcut;
				$lenCutSeq = (length $sequence) - $length;
				if ($lenCutSeq >= $cutoff) {
					substr($sequence,$startcut,$length) = "";
				} else {
					$length = (length $sequence) - $cutoff;
					substr($sequence,$startcut,$length) = "";
					$lenCutSeq = length $sequence;
					last;
				}
			}
			print CUT "$contig\t", "$lenSeq\t", "$lenCutSeq\n";
		}
		print OUT ">$contig\n", "$sequence\n";
    }
}
close(OUT);
close(CUT);


sub get_all_ids {

 grep {!/^__/} keys %{shift->{offsets}}

}

