package Writing;
use strict;
use warnings;

use Exporter qw(import);
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0);
use Calculation;
use Treatment;

our @EXPORT = qw(write_assembly get_taxonomic_ranks write_assembly_component
									get_fasta_and_report_sequence_ena_assembly 
										sequence_ena get_fasta_and_report_sequence_ena_other
											add_to_file);



# write general assembly file
sub write_assembly {
	my ($assembly_report, $genomic_file, $genbank_file, $summary, $repositoryAssembly,
		    $chromosomes_summary, $plasmids_summary, $scaffolds_summary, 
				$contigs_summary, $specific_summary, $components, @header) = @_;
				
	my %hash_informations = ();
	my $seq = "";
	my $fasta_details = "";
	my $genome_name = "";
	my $country = "ND";
	my $GCpercent = -1;
	my $entropy_level = "ND";
	my $tax_id = "ND";
	my $assembly_line;
	my $pubmed_id = "ND";
	
	open(FIC, "<", $assembly_report) or die ("Could not open $!\n");
	while (<FIC>) {
		chomp;
		$_ =~ s/^#*//;
		if ($_ =~ /:/) {
			my @ligne = split(':', $_);
			$ligne[1] =~ s/^\s+//;
			$ligne[1] =~ s/\s+$//;
			$hash_informations{$ligne[0]} = $ligne[1];
		}
		if ($_  =~ /assembled-molecule/) { $assembly_line = $_; }
	}
	close(FIC);
	
	my @header_report = keys(%hash_informations);
	
	open(FILE_SUMMARY, ">>", $summary) or die ("Could not open $!\n");
	foreach my $k(@header) {
		if (grep $_ eq $k, @header_report) {
			my $information = $hash_informations{$k};
			if ($k =~ /Assembly name/) {
				$genome_name = $information;
			}
			if (($information =~ /^\s*$/) || ($information eq "")) {
				print FILE_SUMMARY "ND\t";
			}
			else {
				print FILE_SUMMARY $information . "\t";
			}
		}
		else {
			print FILE_SUMMARY "ND\t";
		}
		
	}
	
	open(FIC2, "<", $genomic_file) or die ("Could not open $!\n");
	while (<FIC2>) {
		chomp;
		if ($_ !~ /^>/) { $seq .= $_; }
	}
	close(FIC2);
	
	if ($hash_informations{' Taxid'} !~ /\s+/) { $tax_id = $hash_informations{' Taxid'} };

	$GCpercent = gc_percent($seq);
	my ($species,$genus,$family,$order,$class,$phylum,$kingdomGB) = get_taxonomic_rank_taxid($tax_id, "rankedlineage.dmp");
	my ($a_percent, $t_percent, $g_percent, $c_percent) = nucleotid_percent($genomic_file);
	
	open(FIC3, "<", $genbank_file) or die ("Could not open $!");
	while(<FIC3>) {
		chomp;
		if ($_ =~ /\/country="(.*)"/) { $country = $1; }
		if ($_ =~ /PUBMED(.*)/) {  $pubmed_id = trim($1); }		
	}
	close(FIC3);	
	
	print FILE_SUMMARY $pubmed_id."\t".$GCpercent."\t".$entropy_level. "\t". $species."\t".$genus."\t".$family."\t".
											$order."\t".$class."\t".$phylum."\t".$kingdomGB."\t". $country. "\t". $a_percent . "\t" .
											$t_percent . "\t" . $g_percent . "\t" .$c_percent  ."\n" ; 
	close(FILE_SUMMARY);
	
	write_assembly_component($genomic_file, $genome_name, $chromosomes_summary,
	$plasmids_summary, $scaffolds_summary, $contigs_summary, $specific_summary, $components);
}
#------------------------------------------------------------------------------
#function allowing to get taxonomic ranks from Genbank file
sub get_taxonomic_rank {
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

	if($arraySize == 7) {
		($species,$genus,$family,$order,$class,$phylum,$kingdomGB) = @classification;
	}
	elsif($arraySize == 4) {
		($species,$class,$phylum,$kingdomGB) = @classification;
	}
  	
	return ($species,$genus,$family,$order,$class,$phylum,$kingdomGB); 
}
#------------------------------------------------------------------------------
# get assembly component
sub write_assembly_component {
	my($multi_fasta, $assembly_name, $chromosomes_summary, $plasmids_summary,
			$scaffolds_summary, $contigs_summary, $specific_summary, $components) = @_;
			
	my $status = "ND";
	my $level = "ND";
	my $gcpercent;
	my $info;
	my $tmp_fasta_file = "/tmp/sequence.txt";
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
			## a factoriser
			if ($desc[0] =~ /chromosome/) {
				$status = "chromosome";
				$info = $seqID . "\t" . $assembly_name ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t"
							. $gcpercent."\t". $a_percent ."\t". $t_percent ."\t". $g_percent ."\t". $c_percent . "\n";
				add_to_file($chromosomes_summary, $info);
			} elsif ($desc[0] =~ /plasmid/) {
				$status = "plasmid";
				$info = $seqID . "\t" . $assembly_name ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t"
							. $gcpercent."\t". $a_percent ."\t". $t_percent ."\t". $g_percent ."\t". $c_percent . "\n";
				add_to_file($plasmids_summary, $info);				
			} elsif ($desc[0] =~ /scaffold/) {
				$status = "scaffold";
				$info = $seqID . "\t" . $assembly_name ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t"
							. $gcpercent."\t". $a_percent ."\t". $t_percent ."\t". $g_percent ."\t". $c_percent . "\n";
				add_to_file($scaffolds_summary, $info);					
			} elsif ($desc[0] =~ /contig/) {
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
	my $tmp_file = "/tmp/fichier.txt";
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
	for my $id(@id_list){
		$url = "https://www.ebi.ac.uk/ena/data/view/$id&display=text&header=true";
		$output = get($url);
		open(FILE, ">>", $report_file) or die("could not open $!");
		print FILE $output;
		print FILE "##########################################################################\n\n";
		close(FILE) or die("could not close $!");
	}
	return ($fasta_file, $report_file);
}
#------------------------------------------------------------------------------
# download ENA
sub sequence_ena {
	my($sequenceID) = @_;
	my $assemblyRep = $sequenceID . "_folder";
	my $fasta_file;
	my $report_file;

	if(-d $assemblyRep){system("rm -rf $assemblyRep");}
	mkdir $assemblyRep;
	
	if($sequenceID =~ /^GCA_.*/){
		($fasta_file, $report_file) = get_fasta_and_report_sequence_ena_assembly($sequenceID);
	}
	else{
		($fasta_file, $report_file) = get_fasta_and_report_sequence_ena_other($sequenceID);
	}
	
	system("mv $fasta_file $assemblyRep");
	system("mv $report_file $assemblyRep");
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
	close(FILE) or die("could not close $!");
	
	$url = "https://www.ebi.ac.uk/ena/data/view/$sequenceID&display=text&header=true";
	$output = get($url);
	$report_file = $sequenceID . "_report.txt";
	open(FILE, ">>", $report_file) or die("could not open $!");
	print FILE $output;
	close(FILE) or die("could not close $!");
	
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
sub get_taxonomic_rank_taxid {
	my($tax_id, $taxonomic_file) = @_;
	my $species = "ND";
	my $genus = "ND";
	my $family = "ND";
	my $order = "ND";
	my $class = "ND";
	my $phylum = "ND";
	my $kingdom = "ND";
	open(TFILE, "<", $taxonomic_file) or die("Could not open $!");
		while(<TFILE>) {
			chomp;
			my @tax_info = split(/\|/, $_);
			if ($tax_info[0] == $tax_id) {
				@tax_info  = trim_array(@tax_info);
				splice(@tax_info, 0, 2);
				my @tmp_array = ($species, $genus, $family, $order, $class, $phylum, $kingdom);
				for(my $i = 0; $i < $#tax_info + 1; $i++) {
					if (length($tax_info[$i]) > 0) { $tmp_array[$i] = $tax_info[$i] }
				}
				close(TFILE);
				return @tmp_array;
			}
		}
}


1;