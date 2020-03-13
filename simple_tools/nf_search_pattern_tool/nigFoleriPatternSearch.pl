#!/usr/bin/env perl

## 12/03.2020
## David Couvin
## Moco Vincent


################################################
## search pattern in nigleria Foleri species
################################################

use strict;
use warnings;

#variables
my $sequencesFile;  # file for sequences

my $patternsFile = "pattern.txt"; # file for patterns

my $summary = "Nigleria_Foleri_summary.xls";


if (@ARGV<1) {
	user_help($0);
	exit 1;
}


for (my $i=0; $i<=$#ARGV; $i++) {

	if ($ARGV[$i]=~/-input/ or $ARGV[$i]=~/-i/) {
		$sequencesFile = $ARGV[$i+1];
    }
	
}

print "Working....\n";

compare_sequences_patterns($sequencesFile, $patternsFile, $summary);

print "End of process \n";

#  compare sequences patterns
sub compare_sequences_patterns {
	my ($sequencesFile, $patternsFile, $summary) = @_;
	
	my %hashSequences = extract_informations($sequencesFile, ">");
	my %hashPatterns = extract_informations($patternsFile, "\\d");

	my %revHashPatterns = reverse %hashPatterns;

	my %hashPatternSeq;
	my %hashSeqPattern;
	
	my $length = values %hashSequences;
	
	print "number of sequences found : $length\n";

	my @ListPattern ;

for my $seq (keys %hashSequences) {
	@ListPattern = match_pattern($hashSequences{$seq}, %hashPatterns);
	@ListPattern = control_pattern(@ListPattern);
	@{$hashSeqPattern{$seq}} = @ListPattern;
}

for my $patt (keys %hashPatterns) {
	@{$hashPatternSeq{$patt}} = ();
}

for my $seq (keys %hashSeqPattern) {

	@ListPattern = @{$hashSeqPattern{$seq}};
	
	foreach my $patt (keys %revHashPatterns) {
		foreach my $seqPatt (@ListPattern) {
			if ($patt eq $seqPatt) { 
				push @{$hashPatternSeq{$revHashPatterns{$patt}}}, $seq;
			}
		}
	}
}
write_summary($summary, %hashPatternSeq);	
}
#------------------------------------------------------------------------------
#  get pattern present in sequence
sub match_pattern {
	my ($seq, %hashPatterns) = @_;
	
	my @matchPattern;
	
	for my $pattern  (values %hashPatterns) {
		if ($seq =~ /$pattern/ ) { push @matchPattern, $pattern; } #  && not grep (/$pattern/, @matchPattern)
	}
	
	return @matchPattern;
}
#------------------------------------------------------------------------------
#  extract sequences and patterns from file and return them as list
sub control_pattern {
	my (@matchPattern) = @_;
	
	my @patternsToDelete;
	my @newMatchPatterns;
	
	foreach my $chekPattern (@matchPattern) {
		foreach my $pattern (@matchPattern) {
			if ( $pattern ne $chekPattern && $pattern =~ /$chekPattern/) {
				push @patternsToDelete, $chekPattern;
			}
		}	
	}
	
	foreach my $pattern (@matchPattern) {
		my $found =  grep {/^$pattern$/} @patternsToDelete;
		if ($found == 0) { push @newMatchPatterns, $pattern; }
	}
	
	return @newMatchPatterns;
}
#------------------------------------------------------------------------------
#  extract sequences and patterns from file and return them as list
sub extract_informations {
	my ($infoFile, $header) = @_;
	
	my %hashInformations;
	my $key;
	
	open(INF, "<", $infoFile) or die "open file error $!:";
	while(<INF>) {
		chomp;
		if ($_ =~ /^$header/) { $key = $_; }
		if ($_ =~ /^[ATGC]{1,}$/) { $hashInformations{$key} = $_; }
	}
	close (INF) or die "close file error $!:";
	
	return %hashInformations;
}
#------------------------------------------------------------------------------
#  write summary
sub write_summary {
	my ($summary, %hashPatternSeq) = @_;
	
	open(SUM, ">", $summary) or die "error open file $!:";
	for my $pattern (keys %hashPatternSeq) {
	
		my @nameSeq = @{$hashPatternSeq{$pattern}};  
	
		print SUM "motif $pattern : \n";
		
		foreach  my $name (@nameSeq) {
			print SUM "$name\n";
		}
		print SUM "\n";
	}
	close(SUM) or die "error close file $!:";
}
#------------------------------------------------------------------------------
#  user help
sub user_help {
	print <<HEREDOC;
	
	Name : 
		$0.
	
	Synopsis :
		a Perl script to get pattern in nigleria Foleri sequences
		
	Usage :
	
	  example : 
	     perl $0 -input sequence.txt 
		 
HEREDOC
}