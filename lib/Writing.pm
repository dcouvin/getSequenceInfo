package Writing;
use strict;
use warnings;

use Exporter qw(import);
use Cwd  qw(abs_path);
use Calculation;
use File::Basename qw(dirname);
use File::Copy qw(cp move);
use File::Path qw(rmtree);
use LWP::Simple qw(get);
use lib dirname(abs_path $0);
use POSIX qw(floor);
use Treatment;

our @EXPORT = qw(write_assembly get_taxonomic_ranks write_assembly_component
									get_fasta_and_report_sequence_ena_assembly 
										sequence_ena get_fasta_and_report_sequence_ena_other
											add_to_file write_html_summary);



# write general assembly file
sub write_assembly {
	my ($assemblyReport, $genomicFile, $genbankFile, $summary, $repositoryAssembly,
		    $chromosomesSummary, $plasmidsSummary, $scaffoldsSummary, 
				$contigsSummary, $specificSummary, $components, $kingdom, $actualOS, @header) = @_;
				
	my %hashInformations = ();
	my $seq = "";
	my $genomeName = "";
	my $country = "na";
	my $GCpercent = -1;
	my $entropyLevel = "na";
	my $taxId = "na";
	my $assemblyLine;
	my $pubmedId = "na";
	my $host = "na";
	my $isoSource = "na";
	my $species = "na";
	my $genus = "na";
	my $family = "na";
	my $order = "na";
	my $class = "na";
	my $phylum = "na";
	
	open(FIC, "<", $assemblyReport) or die "error open file $!:";
	while (<FIC>) {
		chomp;
		$_ =~ s/^#*//;
		if ($_ =~ /:/) {
			my @ligne = split(':', $_);
			if (defined $ligne[1]) {
				$ligne[1] = trim($ligne[1]);
				$hashInformations{$ligne[0]} = $ligne[1];
			}
		}
		if ($_  =~ /assembled-molecule/) { $assemblyLine = $_; }
	}
	close(FIC) or die "close file error $!:";
	
	my @header_report = keys(%hashInformations);
	
	open(FILE_SUMMARY, ">>", $summary) or die "error open file $!:";
	foreach my $k(@header) {
		if (grep $_ eq $k, @header_report) {
			my $information = $hashInformations{$k};
			
			if ($k =~ /Assembly name/) { $genomeName = $information; }
			
			if (($information =~ /^\s*$/) || ($information eq "")) {
				print FILE_SUMMARY "na\t";
			} else { 
				print FILE_SUMMARY $information . "\t";
			}
		} else {
			print FILE_SUMMARY "na\t";
		}
	}
	
	open(FIC2, "<", $genomicFile) or die "Could not open $!:";
	while (<FIC2>) {
		chomp;
		if ($_ !~ /^>/) { $seq .= $_; }
	}
	close(FIC2)  or die "Close file error $!:";
	
	if ($hashInformations{' Taxid'} !~ /\s+/) { $taxId = $hashInformations{' Taxid'} };

	$GCpercent = gc_percent($seq);
	
	if ($actualOS eq "linux") {
		($species, $genus, $family, $order, $class, $phylum) =  get_taxonomic_rank($taxId, "rankedlineage.dmp");
	}
	elsif ($actualOS eq "MSWin32") { 
		($species, $genus, $family, $order, $class, $phylum) = get_taxonomic_rank_genbank($genbankFile);
	}
	
	my ($a_percent, $t_percent, $g_percent, $c_percent) = nucleotid_percent($genomicFile);
	
	open(FIC3, "<", $genbankFile) or die "Could not open $!:";
	while(<FIC3>) {
		chomp;
		if ($_ =~ /\/country="(.*)"/) { $country = trim($1); }
		if ($_ =~ /PUBMED(.*)/) {  $pubmedId = trim($1); }
		if ($_ =~ /\/host="(.*)"/) {  $host = trim($1); }		
		if ($_ =~ /\/isolation_source="(.*)"/) {  $isoSource = trim($1); }		
	}
	close(FIC3) or die "Close file error $!:";	
	
	print FILE_SUMMARY $pubmedId . "\t" . $GCpercent . "\t" . $entropyLevel . "\t" . $species . "\t" . $genus . "\t" . $family ."\t" . 
												$order . "\t" . $class . "\t" . $phylum . "\t" . $kingdom . "\t" . $country . "\t" . $host . "\t" . $isoSource  . "\t" .
													$a_percent . "\t" . $t_percent . "\t" . $g_percent . "\t" . $c_percent  ."\n" ; 
									
	close(FILE_SUMMARY) or die "close file error $!:";
	
	write_assembly_component($genomicFile, $genomeName, $chromosomesSummary,
	$plasmidsSummary, $scaffoldsSummary, $contigsSummary, $specificSummary, $components);
}
#------------------------------------------------------------------------------
# get assembly component
sub write_assembly_component {
	my($multi_fasta, $assembly_name, $chromosomes_summary, $plasmids_summary,
			$scaffolds_summary, $contigs_summary, $specific_summary, $components) = @_;
			
	my $status = "na";
	my $level = "na";
	my $gcpercent;
	my $info;
	my $tmp_fasta_file = "sequence.txt";
	my @desc = ();
	
	# check each sequence from (multi-)fasta file
	my ($seq, $inputfile);
	#my @tabNcounts = ();

	# extract sequence details
	my $seqIO = Bio::SeqIO->new(-format=>'Fasta', -file=>$multi_fasta);
	
	while ($seq = $seqIO->next_seq()) {
		my $seqID = $seq->id; # ID of sequence
		my $seqDesc = $seq->desc; # Description of sequence
		my $globalSeq = $seq->seq;
		my $seqLength = $seq->length();
		
		open(TSEQ, ">", $tmp_fasta_file) or die("Could not open $!");
		print TSEQ $globalSeq;
		close(TSEQ);
		
		my($a_percent, $t_percent, $g_percent, $c_percent) = nucleotid_percent($tmp_fasta_file);
		$gcpercent = gc_percent($globalSeq);
		
		@desc = split(',', $seqDesc);
		
		if ($components) {
			if ($desc[0] =~ /$components/) {
				$info = $seqID . "\t" . $assembly_name ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t"
							. $gcpercent."\t". $a_percent ."\t". $t_percent ."\t". $g_percent ."\t". $c_percent . "\n";
				add_to_file($specific_summary, $info);
			}
		}	
		else {
			if ($desc[1]) {$level = $desc[1];}
			
			if ($desc[0] =~ /chromosome/) {
				$status = "chromosome";
				$info = $seqID . "\t" . $assembly_name ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t"
							. $gcpercent."\t". $a_percent ."\t". $t_percent ."\t". $g_percent ."\t". $c_percent . "\n";
				add_to_file($chromosomes_summary, $info);
			}
			elsif ($desc[0] =~ /plasmid/) {
				$status = "plasmid";
				$info = $seqID . "\t" . $assembly_name ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t"
							. $gcpercent."\t". $a_percent ."\t". $t_percent ."\t". $g_percent ."\t". $c_percent . "\n";
				add_to_file($plasmids_summary, $info);				
			} 
			elsif ($desc[0] =~ /scaffold/) {
				$status = "scaffold";
				$info = $seqID . "\t" . $assembly_name ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t"
							. $gcpercent."\t". $a_percent ."\t". $t_percent ."\t". $g_percent ."\t". $c_percent . "\n";
				add_to_file($scaffolds_summary, $info);					
			} 
			elsif ($desc[0] =~ /contig/) {
				$status = "contig";
				$info = $seqID . "\t" . $assembly_name ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t"
							. $gcpercent."\t". $a_percent ."\t". $t_percent ."\t". $g_percent ."\t". $c_percent . "\n";
				add_to_file($contigs_summary, $info);	
			}
		}	
	}
}
#------------------------------------------------------------------------------
# download fasta sequence and report on ENA with assembly ID
sub get_fasta_and_report_sequence_ena_assembly {
	my($sequenceID) = @_;
	my $tmp_file = "fichier.txt";
	my @id_list = ();
	my $id_chain = "";
	my $fasta_file = "";
	my $report_file = "";
	my $url = "https://www.ebi.ac.uk/ena/data/view/$sequenceID&display=xml";
	my $output = get($url);
	
	open(TMP, ">", $tmp_file) or die("could not open $!");
	print TMP $output;
	close(TMP) or die("could not close $!");
	
	open(TMP, "<", $tmp_file) or die("could not open $!");
	while(<TMP>){
		if($_ =~ /<CHROMOSOME accession="(.*)">/){
		push(@id_list, $1)
		}
	}
	close(TMP) or die("could not close $!");
	
	$id_chain = join(",", @id_list);
	$url = "https://www.ebi.ac.uk/ena/data/view/$id_chain&display=fasta";
	$output = get($url);
	$fasta_file = $sequenceID . ".fasta";
	open(FILE, ">", $fasta_file) or die("could not open $!");
	print FILE $output;
	close(FILE) or die("could not close $!");
	
	
	$report_file = $sequenceID . "_report.txt";
	for my $id (@id_list) {
		$url = "https://www.ebi.ac.uk/ena/data/view/$id&display=text&header=true";
		$output = get($url);
		open(FILE, ">>", $report_file) or die("could not open $!");
		print FILE $output;
		print FILE "##########################################################################\n\n";
		close(FILE) or die("could not close $!");
	}
	
	unlink "fichier.txt" or die " $!: error delete file fichier.txt";
	
	return ($fasta_file, $report_file);
}
#------------------------------------------------------------------------------
# download ENA
sub sequence_ena {
	my($sequenceID) = @_;
	my $assemblyRep = $sequenceID . "_folder";
	my $fastaFile;
	my $reportFile;

	if(-d $assemblyRep) { rmtree($assemblyRep); }
	mkdir $assemblyRep;
	
	if($sequenceID =~ /^GCA_.*/){
		($fastaFile, $reportFile) = get_fasta_and_report_sequence_ena_assembly($sequenceID);
	}
	else {
		($fastaFile, $reportFile) = get_fasta_and_report_sequence_ena_other($sequenceID);
	}
	move($fastaFile, $assemblyRep) or die "move failed: $!";
	move($reportFile, $assemblyRep) or die "move failed: $!";
} 
#------------------------------------------------------------------------------
# download fasta sequence and report on ENA with ENA ID 
sub get_fasta_and_report_sequence_ena_other {
	my($sequenceID) = @_;
	my $fasta_file = "";
	my $report_file = "";
	my $url;
	my $output;
	
	$url = "https://www.ebi.ac.uk/ena/data/view/$sequenceID&display=fasta";
	$output = get($url);
	$fasta_file = $sequenceID . ".fasta";
	open(FILE, ">", $fasta_file) or die("could not open $!");
	print FILE $output;
	close(FILE) or die "could not close $!";
	
	$url = "https://www.ebi.ac.uk/ena/data/view/$sequenceID&display=text&header=true";
	$output = get($url);
	$report_file = $sequenceID . "_report.txt";
	open(FILE, ">>", $report_file) or die("could not open $!");
	print FILE $output;
	close(FILE) or die "could not close $!";
	
	return ($fasta_file, $report_file);
}
#------------------------------------------------------------------------------
# add information to file
sub add_to_file {
	my ($file, $info) = @_;
	open(FILE, ">>",  $file) or die ("Could not open $!");
	print FILE $info;
	close(FILE);
}
#------------------------------------------------------------------------------
#   return taxonomic rank of species by tax id
sub get_taxonomic_rank {
	my($tax_id, $taxonomic_file ) = @_;
	
	my ($species,$genus,$family,$order,$class,$phylum);
	my @tmp_array = ($species, $genus, $family, $order, $class, $phylum);
	
	open(TFILE, "<", $taxonomic_file) or 
		die("Could not open $taxonomic_file: $!");
		
	while(<TFILE>) {
		chomp;
		my @tax_info = split(/\|/, $_);
		
		if ($tax_info[0] == $tax_id) {
			@tax_info  = trim_array(@tax_info);
		
			$tmp_array[0] = $tax_info[1];
			splice(@tax_info, 0, 3);
			
			for(my $i = 1; $i < $#tmp_array + 1; $i++) {
				if (length($tax_info[$i-1]) > 0) { $tmp_array[$i] = $tax_info[$i-1]; }
			}
			close(TFILE) or die "error close $taxonomic_file  $!:";
			return @tmp_array;
		}
	}
	close(TFILE) or die "error close $taxonomic_file  $!:";
}
#------------------------------------------------------------------------------
#   write html summary file
sub write_html_summary {
	my($summary) = @_;
	my $htmlFile = "summary.html";
	my $header = "";
	my @fileToRead = ();
	
	open(HTML, ">", $htmlFile) or die "error open HTML summary $!";
	print HTML "<!DOCTYPE html>\n";
	print HTML "<html>\n";
	print HTML " <head>\n";
	print HTML "  <title>Assembly summary</title>\n";
	print HTML " </head>\n";
	print HTML " <body>\n";
	print HTML "  <h2>Assembly Summary</h2>\n";
	close(HTML) or die "error close HTML summary $!";
	
	open(SUM, "<", $summary) or die "error open tsv summary $!";
	@fileToRead = <SUM>;
	close(SUM) or die "error close tsv summary $!";
	
	$header = splice(@fileToRead, 0, 1);
	
	for my $line (@fileToRead) {
		write_html_table($line, $htmlFile, $header);
	}
	
	open(HTML, ">>", $htmlFile) or die "error open HTML summary $!";
	print HTML " </body>\n";
	print HTML "</html>\n";
	close(HTML) or die "error close HTML summary $!";
}
#------------------------------------------------------------------------------
#   write html table for summary
sub write_html_table {
	my ($line, $htmlFile, $header) = @_;
	
	open(HTML, ">>", $htmlFile) or die "error open HTML summary $!";
	print HTML "  <table  border=\"1\" style=\"margin-bottom: 20px;\">\n";
	close(HTML) or die "error close HTML summary $!";
	add_table_content($line, $htmlFile, $header);
}
#------------------------------------------------------------------------------
#   add information to table
sub add_table_content {
	my ($line, $htmlFile, $header) = @_;
	
	my @assemblyHeader = split(/\t/, $header);
	my @assemblyInfo = split(/\t/, $line);
	my $length = $#assemblyHeader + 1;
	my $lineNb = 0;
	my $fullLine = floor($length / 7);
	my $title = "";
	my $info = "";
	
	open(HTML, ">>", $htmlFile) or die "error open HTML summary $!";
	print HTML "   <tr>\n";
	
	for (my $i = 1; $i < $length + 1; $i++) {
	
		$title = trim($assemblyHeader[$i - 1]);
		print HTML "    <th>$title</th>\n";
		
		if ($i % 7 == 0) {
			my $itr = 0;
			print HTML "   </tr>\n";
			print HTML "   <tr>\n";
			
			while($itr < 7) {
				$info = trim(splice(@assemblyInfo, 0, 1));
				print HTML "    <td>$info</td>\n";
				$itr ++;
			}
			print HTML "   </tr>\n";
			print HTML "   <tr>\n";
			$lineNb++;
		}
		
		if ($length % 7 != 0 && $lineNb == $fullLine) {
		
			$i++;
			
			while ($i < $length + 1) {
				$title = trim($assemblyHeader[$i - 1]);
				print HTML "    <th>$title</th>\n";
				$i++;
			} 
			
			print HTML "   </tr>\n";
			print HTML "   <tr>\n";
			
			my $subLength = $#assemblyInfo + 1;
			
			for(my $j = 0; $j < $subLength; $j++) {
				$info = trim(splice(@assemblyInfo, 0, 1));
				print HTML "    <td>$info</td>\n";
			}
			print HTML "   </tr>\n";
		}
	}
	print HTML "  </table>\n";
	close(HTML) or die "error close HTML summary $!";
}
#------------------------------------------------------------------------------
#getTaxonomicRanks (function allowing to get taxonomic ranks from Genbank file)
sub get_taxonomic_rank_genbank {
	my ($genbank) = @_;
	my ($species,$genus,$family,$order,$class,$phylum,$kingdomGB);

	my $seqio_object = Bio::SeqIO->new(-file => $genbank);
	my $seq_object = $seqio_object->next_seq;

	# legible and long
	my $species_object = $seq_object->species;
	my $species_string = $species_object->node_name;

	# get all taxa from the ORGANISM section in an array
	my @classification = $seq_object->species->classification;
	my $arraySize = @classification;

	if($arraySize == 7){
		($species,$genus,$family,$order,$class,$phylum,$kingdomGB) = @classification;
	}
	elsif($arraySize == 4){
		($species,$class,$phylum,$kingdomGB) = @classification;
	}
  	
	return ($species,$genus,$family,$order,$class,$phylum); 
}

1;