#!/usr/bin/env perl 
use strict;

my $outfile = shift or die $!;

my %data = ();

my @warn = ();

foreach my $infile (@ARGV){
	print STDOUT "Reading ".$infile."\n";
	open (INFILE, "$infile") or die $!;
	my $crisprct=0;
	while(my $line=<INFILE>){
		chomp($line);
		next if $line =~/#/;
		my @tabs = split("\t",$line);
		my $attr = pop @tabs;
		my $tabstr = join(";",@tabs);
		my @attrs = split(";",$attr);
		# note=CRISPR with 19 repeat units;rpt_family=CRISPR;rpt_type=direct
		if ($attrs[0] =~ /^note=CRISPR/){
			$crisprct++;
			$attr="ID=CRISPR_".$crisprct.";locus_tag=CRISPR_".$crisprct.";".join(";",@attrs);
			@attrs = split(";",$attr);
		}
		my $id = "";
		my @tmpattr = ();

		foreach my $anno (@attrs){
			if($anno =~ /^ID=/){
				($id = $anno) =~ s/^ID=//;
			}else{
				#print STDERR $anno."\n";
				if($anno ne ""){
					my @annname=split("=",$anno);
					if($annname[0] =~ /-/){
						my $old_annname=$annname[0];
						$annname[0] =~ s/-/_/g;
						$anno=join("=",@annname);
						push @warn,"$old_annname was replaced by $annname[0].";
					}
					push @tmpattr,$anno;
				}
			}
		}
		if(exists($data{$tabstr})){
			if($data{$tabstr}{'id'} eq $id){
				push @{$data{$tabstr}{'att'}},@tmpattr;
			}else{
				print STDERR "# WARNING: two different IDs for the same locus.\n";
			}
		}else{
			$data{$tabstr}{'id'} = $id;
			@{$data{$tabstr}{'att'}} = @tmpattr;
		}
	}
	close(INFILE);
}

my @uniqueWarn = do { my %seen; grep { !$seen{$_}++ } @warn };
foreach my $warnout (@uniqueWarn){
	print STDOUT $warnout."\n";
}	
	
open (my $OUT, '>', $outfile);
foreach my $outline (keys(%data)){
	my @tabs = split(";",$outline);
	my @tmpattr = "ID=".($data{$outline}{'id'});
	my @unique = do { my %seen; grep { !$seen{$_}++ } @{$data{$outline}{'att'}} };
	push @tmpattr,@unique;
	print $OUT join("\t", @tabs)."\t".join(";", @tmpattr)."\n";
}

close($OUT)



