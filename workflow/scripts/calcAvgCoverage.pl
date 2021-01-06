#!/usr/bin/perl
#
use strict;

my $input=$ARGV[0];
my $contigs=$ARGV[1] if $ARGV[1];

my %ids=();
my @ids=();

if ($contigs){
    open(CONTIGS,$contigs) or die $!;
    while (my $line=<CONTIGS>){
        chomp($line);
        next if $line !~/\S/;
        if ($line=~/^>(.+)$/){
            my $id=$1;
            push @ids, $id;
        } 
    }
    close(CONTIGS);
}

open(FILE,$input) or die $!;
while(my $line=<FILE>){
	chomp($line);
    next if $line !~/\S/;
	my @F=split(/\t/,$line);
	if(exists($ids{$F[0]})){
		$ids{$F[0]}+=$F[1]*$F[4];
	}else{
		$ids{$F[0]}=$F[1]*$F[4];
        push @ids,$F[0] if !$contigs;
	}
}
close(FILE);

foreach my $id (@ids){
    if(exists($ids{$id})){
        printf("%s\t%.5f\n",$id,$ids{$id});
    }else{
        printf("%s\t%.5f\n",$id,0);    
    }
}


