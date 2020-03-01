package Calculation;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT = qw(nucleotid_percent  log2 entropy gc_percent);

# compute percentage of nucleotid
sub nucleotid_percent {
	my ($fastaFile) = @_;
	my $A = 0;
	my $T = 0;
	my $G = 0;
	my $C = 0;
	my $length = 0;
	open (FASTA, "<", $fastaFile) or die ("Could not open $!");
	while (<FASTA>) {
		chomp;
		my $line = $_;
		if ($line !~ />/) {
			my @seq = split(//, uc($line));
			for my $NC (@seq) {
				$length++;
				if ($NC eq "A") {$A++;}
				if ($NC eq "T") {$T++;}
				if ($NC eq "G") {$G++;}
				if ($NC eq "C") {$C++;}
			}
		}
	}
	my $A_PERCENT = ($A/$length) * 100;
	my $T_PERCENT = ($T/$length) * 100;
	my $G_PERCENT = ($G/$length) * 100;
	my $C_PERCENT = ($C/$length) * 100;
	return ($A_PERCENT, $T_PERCENT, $G_PERCENT, $C_PERCENT)
}
#------------------------------------------------------------------------------
sub log2 {
  my $n = shift;
  return (log($n)/ log(2));
}
# Function allowing to calculate conservation of DRs based on entropy
sub entropy {
  my($file) = @_; #file given as parameter is an alignment Fasta file
  
  open F, $file or die ("an error occurred while opening $file\n");
  my @lines = <F>;
  close F;

  my (%words, $total, @text);

  my $seqLength=0;
  my @tableRows=();
  my @tableCols=();

  my @table=();

  my $countLine = 0;

  foreach my $line (@lines) {
	chomp $line;

	my @words = split /[^a-zA-Z]+/, $line;
	#Replace to treat Fasta sequence
	#Retrieve each letter by column
	if ($line !~ />/)
	{
		$seqLength=length($line);
		#print "Line is : $line\n";
		for (my $i = 0; $i < $seqLength; $i++) {
			my $value = substr($line, $i, 1);
			$table[$countLine][$i] = $value;

			#Instanciate $tableCols with values from columns
			$tableCols[$i][$countLine] = $value;
		}		
		$countLine++;
	}
  }

  # entropy
  my $sumEntropy =0;

  for (my $j = 0; $j < $seqLength; $j++) {	
	my $sizeTable = @{$tableCols[$j]};
	my %element =();

	foreach my $elem (@{$tableCols[$j]}) {
		if ($elem eq "-") { 
		  $element{$elem} = 0; 
		} 
             	else{
		  $element{$elem}++;
		}
	}
	my $entropy;
	foreach my $word (@{$tableCols[$j]}) {
		if ($word ne "-"){
		  my $prob = $element{$word} / $sizeTable;
		  $entropy += log2($prob); 
		}
	}

	$entropy *= -1;
	$entropy = $entropy/$sizeTable;

	#Calculate Sum of entropies and divide by $seqLength ... Adding (1- $entropy) to get percentage-like results
	$sumEntropy = $sumEntropy + (1 - $entropy); 
  }

  #PRINT FINAL Entropy (Sum + "Mean")
  my $finalResult = ($sumEntropy / $seqLength) * 100; 
  
  if($finalResult<0){ $finalResult = 0; }
  
  return $finalResult;
}
#------------------------------------------------------------------------------
# compute GC pourcent
sub gc_percent {
	my ($seq) = @_;
	
	my @charSeq = split(//, uc($seq));
	my %hashFlank = ();

	foreach my $v (@charSeq) {
		$hashFlank{$v} += 1;  
	}
  
	if (! $hashFlank{'G'}) { hashFlank{'G'} = 0;}
	if (! $hashFlank{'C'}) { $hashFlank{'C'} = 0;}

	if(length($seq) == 0) {
		return 0;
	}
	else {
		return (($hashFlank{'G'} + $hashFlank{'C'}) / (length($seq))) * 100;
	}
	
}

1;