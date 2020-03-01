#!/usr/bin/perl

use strict;
use warnings;

use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0) . '/lib';
use Bio::SeqIO;
use Bio::Species;
use Date::Calc qw(:all);
use File::Copy qw(cp move);
use LWP::Simple qw(get);
use Net::FTP;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
 use Archive::Tar;
use Writing;
use Download;

##########################################################################################
##  tool to find genomic informations on NCBI and ENA
##  version : 1.0
##########################################################################################

### main program
# Date and time of the current day (Beginning)
my ($start_year,$start_month,$start_day, $start_hour,$start_min,$start_sec) = Today_and_Now();

print "##################################################################\n";
print "## Welcome to program: $0!\n";
print "## Date: $start_year-$start_month-$start_day, $start_hour:$start_min:$start_sec\n";
print "## version : 1.0\n";
print "## « Copyright 2019 David Couvin, Moco Vincent »\n";
print "## licence GPL-3.0-or-later\n";
print "## This program is free software: you can redistribute it and/or modify\n";
print "##  it under the terms of the GNU General Public License as published by\n";
print "##  the Free Software Foundation, either version 3 of the License, or\n";
print "##  (at your option) any later version.\n";
print "##\n";
print "##  This program is distributed in the hope that it will be useful,\n";
print "##  but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
print "##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
print "##  GNU General Public License for more details.\n";
print "##\n";
print "##  You should have received a copy of the GNU General Public License\n";
print "##  along with this program.  If not, see <https://www.gnu.org/licenses/>. 1\n";
print "##################################################################\n\n";

# parameters
my $version = "1.0";

my $directory = "refseq";
	
my $kingdoms = ""; # kingdom of organism

my $releaseDate = ""; # sequence are download from this release date

my $components; # components of the assembly

my $species = ""; # species name

my $getSummary; # indicates if a new assembly report is required

my $representation = ""; # assembly level of the genome

my $quantity; # number of sequences to download

my $sequenceID;

my $ftpServor = "ftp.ncbi.nlm.nih.gov"; 

if (@ARGV<1) {
	help_user_simple($0);
	exit 1;
}

if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
	help_user_advance($0);
	exit 0;
}

if ($ARGV[0] eq "-version" || $ARGV[0] eq "-v") {
	program_version($0);
	exit 0;	
}

# requirements
for (my $i=0; $i<=$#ARGV; $i++) {
    	if ($ARGV[$i]=~/-kingdom/ or $ARGV[$i]=~/-k/) {
			$kingdoms = $ARGV[$i+1];
    	}
    	elsif ($ARGV[$i]=~/-directory/ or $ARGV[$i]=~/-dir/) {
      		$directory = $ARGV[$i+1];
    	}
    	elsif ($ARGV[$i]=~/-date/) {
      		$releaseDate = $ARGV[$i+1];
    	}
		elsif ($ARGV[$i]=~/-getSummary/ or $ARGV[$i]=~/-get/) {
		if (-e "assembly_summary.txt") { system("rm -f assembly_summary.txt");}
      		$getSummary = 1;
    	}
		elsif ($ARGV[$i]=~/-species/ or $ARGV[$i]=~/-s/ or $ARGV[$i]=~/-S/) {
			$species = $ARGV[$i+1];
			$species =~ tr/ /_/;
		}
		elsif ($ARGV[$i]=~/-representation/ or $ARGV[$i]=~/-r/ or $ARGV[$i]=~/-R/) {
			$representation = $ARGV[$i+1];
		}
		elsif ($ARGV[$i]=~/-components/ or $ARGV[$i]=~/-c/ or $ARGV[$i]=~/-C/) {
			$components = $ARGV[$i+1];
		}
		elsif ($ARGV[$i]=~/-quantity/ or $ARGV[$i]=~/-q/ or $ARGV[$i]=~/-Q/) {
			$quantity = int($ARGV[$i+1]);
		}
		elsif ($ARGV[$i]=~/-ena/ or $ARGV[$i]=~/-ENA/) {
			$sequenceID = $ARGV[$i+1];
		}
}

# Table containing kingdoms individually
my @kingdomTab = split(/,/, $kingdoms); 

foreach my $kingdom (@kingdomTab) {
	$kingdom = lc($kingdom);
	if ($kingdom eq "viruses") {$kingdom = "viral";}
	get_assembly_summary_species($releaseDate, $directory, $kingdom, $species, $representation);	
}

if ($sequenceID) {sequence_ena($sequenceID);}

my ($end_year,$end_month,$end_day, $end_hour,$end_min,$end_sec) = Today_and_Now();

my ($D_y,$D_m,$D_d, $Dh,$Dm,$Ds) =
      Delta_YMDHMS($start_year,$start_month,$start_day, $start_hour, $start_min, $start_sec,
                   $end_year, $end_month, $end_day, $end_hour,$end_min,$end_sec);



### subroutine 
# display global help document
sub help_user_simple {
	my $programme = shift @_;
	print STDERR  "Usage : perl $programme -k XXX -s \"XXX\"  -r \"XXX\" -date yyyy-mm-dd -get  \n";
	print "type -version or -v to get actual version\n";
	print "type -help or -h to get full help\n";
}
#------------------------------------------------------------------------------
# display full help document
sub help_user_advance {
	print <<HEREDOC;
	
	Name : 
		$0.
	
	Synopsis :
		a Perl script to get sequences informations
		
	Usage :
	
	  example : 
	     perl $0 -k bacteria -s "Helicobacter pylori"  -r "Complete Genome" -date 2019-06-01 -get 
	     perl $0 -k bacteria  -q 1  -r "Complete Genome" -date 2019-06-01 -get
						 
		
	Kingdom :
		archaea
		bacteria
		fungi
		invertebrate
		plant
		protozoa
		vertebrate_mammalian
		vertebrate_other
		viral
		
	General :
		-help or -h  
		-version or -v display actual program version
			
	Options :
		-get to obtain a new assembly summary
		perl $0 -k "XXX" -s "XXX"  -r "XXX" -date yyyy-mm-dd -get
		
		-k kingdom of the organism
		perl $0 -k "bacteria" -s"XXX" -r "XXX" -date  yyyy-mm-dd -get 
		
		-s specific species must be combin with -k option
		perl $0 -k "bacteria"  -s "Helicobacter pylori" -r "XXX" -date  yyyy-mm-dd -get
		
		-date  sequences are search from this date
		perl $0 -k "XXX"  -s "XXX" -r "XXX" -date  2019-06-01 -get
		
		-r allow to display sequences by a specific assembly level
		perl $0 -k "XXX"  -s "XXX" -r "Complete Genome" -date yyyy-mm-dd -get
			
		-r example Complete Genome, Chromosome, Scaffold, Contig 
		
		-q quantity of assembly if kingdom is specified but not species
		perl $0 -k "XXX"  -q 10  -r "XXX" -date  yyyy-mm-dd -get
		
		-c specific component of the assembly avoid to download all type of
		component of the assembly
		perl $0 -k "XXX"  -q 00  -r "XXX" -c chromosome -date  yyyy-mm-dd -get
		
		it is possible to use the -c option when only kingdom is mention
		perl getSequenceInfo2.pl -k XXX -q 00 -r "XXX" -c XXX  -date yyyy-mm-dd
		
HEREDOC
}
#------------------------------------------------------------------------------
# display program version 
sub program_version {
	my $programme = shift @_;
	print "\n $programme, version : $version\n";
	print "\n A perl script to get sequences informations\n";
}
#------------------------------------------------------------------------------
sub get_assembly_summary_species {
	my($releaseDate, $directory, $kingdom, $species, $representation) = @_;

	# assembly_summary.txt file from NCBI FTP site
	my $assembly_summary = "ftp://ftp.ncbi.nlm.nih.gov/genomes/$directory/$kingdom/assembly_summary.txt"; 
	
	# lineage folder
	my $lineage_folder = "ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz";

	# new assembly_summary information after modification
  	my $newAssemblySummary = "assembly_summary_$kingdom". $releaseDate.".tsv";

	# Assembly semicolon kingdom file
	my $newAssemblySemicolon = "assembly_semicolon_$kingdom".$releaseDate.".txt";

	# get assembly_summary.txt
	if ($getSummary) { system("wget $assembly_summary"); }
	
	# create $kingdomRep
	my $holdRep = ""; # to check if Repository already exist
	
	
	my $kingdomRep = $kingdom."_".$start_year."_".$start_month."_".$start_day;
	mkdir $kingdomRep unless -d $kingdomRep;
	
	
	# Repository for Assembly
	my $repositoryAssembly = "./assembly_repository_".$species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
	mkdir $repositoryAssembly unless -d $repositoryAssembly;
	
	$holdRep = "./" . $kingdomRep . "/" . $repositoryAssembly;
	if (-d $holdRep) { system("rm -rf $holdRep") }
	
	my $specificRep; # Repository for required component
	my $plasmidsRep; # Repository for plasmids
	my $chromosomesRep; # Repository for chromosomes
	my $scaffoldsRep; # Repository for scaffolds
	my $contigsRep; # Repository for contigs
	
	if ($components) {
		$specificRep = "./" . $components."_".$species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $specificRep unless -d $specificRep;
	}
	else {
		$plasmidsRep = "./plasmids_". $species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $plasmidsRep unless -d $plasmidsRep;
	
		$chromosomesRep = "./chromosomes_" . $species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $chromosomesRep unless -d $chromosomesRep;
		
		$scaffoldsRep = "./scaffolds_" . $species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $scaffoldsRep unless -d $scaffoldsRep;
		
		$contigsRep = "./contigs_" . $species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $contigsRep unless -d $contigsRep;
	}	
	
	# Count how many assemblies and sequences were found with given criteria
	# my $countAssemblies = 0;
   	# my $countSequences = 0;
	my %assemblyReportList;
	
	if (-e "assembly_summary.txt") {
		
		# initialiaze tar manipulation
		my $tar = Archive::Tar->new;
		
		# download taxdump folder
		system("wget $lineage_folder");
		$tar->read("new_taxdump.tar.gz");
		$tar->extract_file("rankedlineage.dmp");
		
		# Read file
		open (SUM, "assembly_summary.txt") or die "open : $!";
		while(<SUM>) {
			chomp;
			my @tab = split('\t', $_);	
			if ($_ !~  m/^#/ && ($tab[11] eq $representation) && ($tab[13] =~  m/Full/)) {
			
				$species =~ tr/_/ /;    #  delete the underscore in the species name
				
				if (($tab[7] =~ /$species/ && $kingdom ne "") or ($kingdom ne "" && $species eq "")) { 
					my @gcfInfo = split(/\//, $tab[19]);  
					my $gcfName = pop(@gcfInfo);
					my $realDate = $tab[14];
					$realDate =~ s/\//-/g;
					
					my $genbankFile = $tab[19] . "/" . $gcfName . "_genomic.gbff.gz";
					my $dnaFile = $tab[19] . "/" . $gcfName . "_genomic.fna.gz";
					my $assemblyReport = $tab[19] . "/" . $gcfName . "_assembly_report.txt";
					
					if ($releaseDate) {
						if ($realDate gt $releaseDate) {
						
							if ($quantity) { $quantity -= 1; }
							
							#$genbankFile = obtain_file($ftpServor, $genbankFile);
							#$dnaFile = obtain_file($ftpServor, $dnaFile);
							#$assemblyReport = obtain_file($ftpServor, $assemblyReport);
							
							# download_file($ftpServor, $dnaFile);
							# download_file($ftpServor, $genbankFile);
							#download_file($ftpServor, $assemblyReport);
							
							# download sequences and check number of "N" characters
							system("wget $dnaFile");
							my $fileFasta = $gcfName."_genomic.fna.gz";
							if (-e $fileFasta) {
								system("gunzip $fileFasta");
								$fileFasta = $gcfName."_genomic.fna";
							}
							move($fileFasta, $repositoryAssembly) or die "move failed: $!";
							
							# download genome report
							my $fileReport =  $gcfName."_assembly_report.txt";
							system("wget $assemblyReport");
							if (-e $fileReport) {
								system("gunzip $fileReport");
								my $fileInformations = $gcfName."_informations.xls";
							}
							move($fileReport, $repositoryAssembly) or die "move failed: $!";
						
							# download genbank files
							system("wget $genbankFile");
							my $fileGenbank = $gcfName."_genomic.gbff.gz";
							if (-e $fileGenbank) { system("gunzip $fileGenbank"); }
							$fileGenbank = $gcfName."_genomic.gbff";
							move($fileGenbank, $repositoryAssembly) or die "move failed: $!";
							
							# association report and fasta
							my $fileFasta_Genbank = $fileFasta . "," . $fileGenbank;
							$assemblyReportList{$fileReport} = $fileFasta_Genbank;
						}
					}
				}	
			}
			if (defined $quantity && $quantity == 0) {
				$quantity = undef;
				last;
			}
		}
		
		close(SUM) or die("close file error $!");
		
		if (!keys %assemblyReportList) {
			print "####################################################\n";
			print "Actual requirements are not found or are invalid for the database\n";
			print "It could be an error on the species name or the date\n";
			print "Check if the name is correctly write ex : Escherichia coli\n";
			print "The date could be too advances try again with an earlier date\n";
			print "Please try again with other requirements\n";
			print "####################################################\n\n";
			exit();
		}
		
		# write summary files 
		my @keysList = keys(%assemblyReportList);
		my $summary = "./summary.xls";
		my $specific_summary;
		my $plasmids_summary;
		my $chromosomes_summary;
		my $scaffolds_summary;
		my $contigs_summary;
		
		if ($components) { $specific_summary =  "./". $components . "_summary.xls"; }
		else {
			$plasmids_summary = "./plasmids_summary.xls";
			$chromosomes_summary = "./chromosomes_summary.xls";
			$scaffolds_summary = "./scaffolds_summary.xls";
			$contigs_summary = "./contigs_summary.xls";
		}
		my $file_report = "./" . $repositoryAssembly ."/" . $keysList[0]; 
		my @summary_list = ($plasmids_summary, $chromosomes_summary, $scaffolds_summary, $contigs_summary);
		my @header = ();
		
		open(FILE, $file_report) or die ("error $!");
		while(<FILE>) {
			chomp;
			$_ =~ s/^#*//;
			if($_ =~ /:/){
				my @ligne = split(':', $_);
				push(@header, $ligne[0]);
			}
		}
		close(FILE);
		open(HEAD, ">", $summary) or die (" error $!\n");
		foreach(@header) {
			print HEAD uc($_) . "\t";
		}
		
		print HEAD "PUBMED\tGC PERCENT\tENTROPY\tSPECIES\tGENUS\tFAMILY\tORDER\tCLASS\t". 
							"PHYLUM\tKINGDOM\tCOUNTRY\tA_PERCENT\tT_PERCENT\tG_PERCENT\tC_PERCENT\n";
		close(HEAD);
		
		if ($components) {
			open(SUM, ">>", $specific_summary) or die ("Could not open $!");
				print SUM "ID\tASSEMBLY\tDESCRIPTION\tLENGTH\tSTATUS\tLEVEL\t" . 
							"GC_PERCENT\tA_PERCENT\tT_PERCENT\tG_PERCENT\tC_PERCENT\n";	
			close(SUM);	
		} else {
			for my $sum(@summary_list) {
				open(SUM, ">>", $sum) or die ("Could not open $!");
				print SUM "ID\tASSEMBLY\tDESCRIPTION\tLENGTH\tSTATUS\tLEVEL\t" . 
							"GC_PERCENT\tA_PERCENT\tT_PERCENT\tG_PERCENT\tC_PERCENT\n";	
				close(SUM);
			}
		}	
		
		for my $file(@keysList) {
			my $file1 = "./" . $repositoryAssembly ."/". $file;
			my @fasta_genbank = split(",", $assemblyReportList{$file});
			my $ext_fasta = $fasta_genbank[0];
			my $ext_genbank = $fasta_genbank[1];
			my $file2 = "./" . $repositoryAssembly ."/". $ext_fasta;
			my $file3 = "./" . $repositoryAssembly ."/". $ext_genbank;
			
			write_assembly($file1, $file2, $file3, $summary, $repositoryAssembly,
				$chromosomes_summary, $plasmids_summary, $scaffolds_summary,
					$contigs_summary, $specific_summary, $components, @header);
		}
		
		move($summary, $repositoryAssembly) or die "move failed: $!";
		if ($components) {
			move($specific_summary, $specificRep) or die "move failed: $!";
			move($specificRep, $repositoryAssembly . "/" . $specificRep) or die "move failed: $!";
		}
		else {
			move($plasmids_summary, $plasmidsRep) or die "move failed: $!";
			move($chromosomes_summary, $chromosomesRep) or die "move failed: $!";
			move($scaffolds_summary, $scaffoldsRep) or die "move failed: $!";
			move($contigs_summary, $contigsRep) or die "move failed: $!";
			move($plasmidsRep, $repositoryAssembly . "/" . $plasmidsRep) or die "move failed: $!";
			move($chromosomesRep, $repositoryAssembly . "/" . $chromosomesRep) or die "move failed: $!";
			move($scaffoldsRep, $repositoryAssembly . "/" . $scaffoldsRep) or die "move failed: $!";
			move($contigsRep, $repositoryAssembly . "/" . $contigsRep) or die "move failed: $!";
		}
		move( $repositoryAssembly, $kingdomRep."/".$repositoryAssembly) or die "move failed: $!";
		clean_repository();
	} 
}


