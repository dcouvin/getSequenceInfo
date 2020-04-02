#!/usr/bin/perl

use strict;
use warnings;

use Archive::Tar;
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0) . '/lib';
use Bio::SeqIO;
use Bio::Species;
use Date::Calc qw(:all);
use File::Copy qw(cp move);
use File::Path qw(rmtree);
use Net::FTP;
use Treatment;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
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
print "## it under the terms of the GNU General Public License as published by\n";
print "## the Free Software Foundation, either version 3 of the License, or\n";
print "##  (at your option) any later version.\n";
print "##\n";
print "## This program is distributed in the hope that it will be useful,\n";
print "## but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
print "## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
print "## GNU General Public License for more details.\n";
print "##\n";
print "## You should have received a copy of the GNU General Public License\n";
print "## along with this program.  If not, see <https://www.gnu.org/licenses/>. 1\n";
print "##################################################################\n\n";

# parameters
my $version = "1.0";

my $directory = "refseq";
	
my $kingdom = ""; # kingdom of organism

my $releaseDate = "0000-00-00"; # sequence are download from this release date

my $components; # components of the assembly

my $species = ""; # species name

my $getSummary; # indicates if a new assembly report is required

my $representation = ""; # assembly level of the genome

my $quantity; # number of sequences to download

my $sequenceID;

my $ftpServor = "ftp.ncbi.nlm.nih.gov";

my $fldSep; # folder seperation change by OS 

my @availableKingdoms = (
	"archaea",
	"bacteria",
	"fungi",
	"invertebrate",
	"plant",
	"protozoa",	
	"vertebrate_mammalian",
	"vertebrate_other",
	"viral"
);  # list of available kingdoms

my $actualOS; # OS of the computer

my $mainFolder; # folder where the assembly are place

my $assemblyTaxid = ""; # taxid for assembly



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
    if ($ARGV[$i]=~/-kingdom/i or $ARGV[$i]=~/-k/i) {
		$kingdom = $ARGV[$i+1];
    }
    elsif ($ARGV[$i]=~/-directory/i or $ARGV[$i]=~/-dir/i) {
    	$directory = $ARGV[$i+1];
    }
    elsif ($ARGV[$i]=~/-date/i) {
    	$releaseDate = $ARGV[$i+1];
    }
	elsif ($ARGV[$i]=~/-getSummary/i or $ARGV[$i]=~/-get/i) {
		if (-e "assembly_summary.txt") { unlink "assembly_summary.txt" or die "$! fail"; }
    	$getSummary = 1;
    }
	elsif ($ARGV[$i]=~/-species/i or $ARGV[$i]=~/-s/i) {
		$species = $ARGV[$i+1];
	}
	elsif ($ARGV[$i]=~/-representation/i or $ARGV[$i]=~/-r/i) {
		$representation = $ARGV[$i+1];
	}
	elsif ($ARGV[$i]=~/-components/i or $ARGV[$i]=~/-c/i) {
		$components = $ARGV[$i+1];
	}
	elsif ($ARGV[$i]=~/-quantity/i or $ARGV[$i]=~/-q/i) {
		$quantity = int($ARGV[$i+1]);
	}
	elsif ($ARGV[$i]=~/-ena/i) {
		$sequenceID = $ARGV[$i+1];
	}
	elsif ($ARGV[$i]=~/-output/i or $ARGV[$i]=~/-o/i) {
		$mainFolder = $ARGV[$i+1];
	}
	elsif ($ARGV[$i]=~/-taxid/i) {
		$assemblyTaxid = $ARGV[$i+1];
	}
}

# define folcer separator and OS
if ($^O =~ /linux/) { 
	$fldSep = "/";
	$actualOS = "linux";
}
elsif ($^O =~ /MSWin32/) { 
	$fldSep = "\\";
	$actualOS	= "MSWin32";
}


print "Working ...\n"; 

if ($kingdom eq "viruses") { $kingdom = "viral"; }

if (grep(/^$kingdom$/i, @availableKingdoms)) {
	if ($species ne "") {
		my @speciesList = split(/,/, $species);
		
		foreach my $actualSpecies (@speciesList) {
			get_assembly_summary_species($releaseDate, $directory, $kingdom, $actualSpecies,
			$representation, $fldSep, $actualOS, $mainFolder, $assemblyTaxid);
		}
	}
	elsif ($assemblyTaxid ne "") {
		my @taxidList = split(/,/, $assemblyTaxid);
		
		foreach my $actualID (@taxidList) {
			get_assembly_summary_species($releaseDate, $directory, $kingdom, $species,
			$representation, $fldSep, $actualOS, $mainFolder, $actualID);
		}
	} 
	else {
		get_assembly_summary_species($releaseDate, $directory, $kingdom, $species,
		$representation, $fldSep, $actualOS, $mainFolder, $assemblyTaxid);
	}	
}


if ($sequenceID) { sequence_ena($sequenceID); }


my ($end_year,$end_month,$end_day, $end_hour,$end_min,$end_sec) = Today_and_Now();

my ($D_y,$D_m,$D_d, $Dh,$Dm,$Ds) =
      Delta_YMDHMS($start_year,$start_month,$start_day, $start_hour, $start_min, $start_sec,
                   $end_year, $end_month, $end_day, $end_hour,$end_min,$end_sec);

print "End of process Date: $start_year-$start_month-$start_day, $start_hour:$start_min:$start_sec\n";
print "Execution time: $Dh:$Dm:$Ds\n";

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
		
		-q quantity of assembly download the maximum possible if the number is to high
		perl $0 -k "XXX"  -q 10  -r "XXX" -date  yyyy-mm-dd -get
		
		-c specific component of the assembly avoid to download all type of
		component of the assembly
		perl $0 -k "XXX"  -q 00  -r "XXX" -c chromosome -date  yyyy-mm-dd -get
		
		it is possible to use the -c option when only kingdom is mention
		perl $0.pl -k XXX -q 00 -r "XXX" -c XXX  -date yyyy-mm-dd
		
		-ena allow to download report and fasta file by ena ID 
		perl $0.pl -ena XXX
		
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
	my ($releaseDate, $directory, $kingdom, $species, $representation, $fldSep, $actualOS, $mainFolder, $assemblyTaxid) = @_;
	
	# assembly_summary.txt file from NCBI FTP site
	my $assemblySummary = "/genomes/$directory/$kingdom/assembly_summary.txt"; 
	
	# lineage folder
	my $lineage_file = "/pub/taxonomy/new_taxdump/new_taxdump.tar.gz";

	# new assembly_summary information after modification
  	my $newAssemblySummary = "assembly_summary_$kingdom". $releaseDate.".tsv";

	# Assembly semicolon kingdom file
	my $newAssemblySemicolon = "assembly_semicolon_$kingdom".$releaseDate.".txt";
	
	# allow to check old summary download
	my $oldKingdom = ""; 
	
	# print "before download summary\n";
	
	# check assembly summary download
	if ($getSummary) { 
		download_file($ftpServor, $assemblySummary);
		open(KIN, ">", "kingdom.txt") or die "error open file $!:";
		print KIN $kingdom;
		close(KIN) or die "error close file $!:";
	}
	
	elsif ($kingdom ne "" && $species eq  "") {
		open(KIN, "<", "kingdom.txt") or die "error open file $!:";
		$oldKingdom = trim(<KIN>);
		close(KIN) or die "error close file $!:";
		
		if ($kingdom ne $oldKingdom) { 
			download_file($ftpServor, $assemblySummary);
			open(KIN, ">", "kingdom.txt") or die "error open file $!:";
			print KIN $kingdom;
			close(KIN) or die "error close file $!:";
		}
	}
	# print "after download summary\n";
	
	# create $kingdomRep
	# to check if Repository already exist
	my $oldRep = "";
	
	my $kingdomRep;

	if (defined $mainFolder) {
		$kingdomRep = $mainFolder; 
	}
	else { 
		$kingdomRep = $kingdom."_".$start_year."_".$start_month."_".$start_day; 
	}
	#my $kingdomRep = $kingdom."_".$start_year."_".$start_month."_".$start_day;
	mkdir $kingdomRep unless -d $kingdomRep;
	
	# Repository for Assembly
	my $repositoryAssembly = "assembly_repository_".$species."_".$assemblyTaxid."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;

	mkdir $repositoryAssembly unless -d $repositoryAssembly;
	
	$oldRep = "." . $fldSep . $kingdomRep . $fldSep . $repositoryAssembly;

	if (-d $oldRep) { rmtree($oldRep) }
	# print "after remove tree\n";
	
	my $specificRep; # Repository for required component
	
	my $plasmidsRep; # Repository for plasmids
	
	my $chromosomesRep; # Repository for chromosomes
	
	my $scaffoldsRep; # Repository for scaffolds
	
	my $contigsRep; # Repository for contigs
	
	if ($components) {
		$specificRep = $components."_".$species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $specificRep unless -d $specificRep;
	} else {
		$plasmidsRep = "plasmids_". $species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $plasmidsRep unless -d $plasmidsRep;
	
		$chromosomesRep = "chromosomes_" . $species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $chromosomesRep unless -d $chromosomesRep;
		
		$scaffoldsRep = "scaffolds_" . $species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $scaffoldsRep unless -d $scaffoldsRep;
		
		$contigsRep = "contigs_" . $species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $contigsRep unless -d $contigsRep;
	}	
	
	my %assemblyReportList;
	
	if (-e "assembly_summary.txt") {
		
		if ($actualOS eq "linux") {
		
			# initialiaze tar manipulation
			my $tar = Archive::Tar->new;
	
			# download taxdump folder
			download_file($ftpServor, $lineage_file);
			$tar->read("new_taxdump.tar.gz");
			$tar->extract_file("rankedlineage.dmp");
		}
		
		# Read file
		open (SUM, "assembly_summary.txt") or die "open assembly_summary.txt : $!";
		while(<SUM>) {
		
			chomp;
			my @tab = split('\t', $_);	
			
			if ($_ !~  m/^#/ && $tab[11] eq $representation && $tab[13] =~  m/Full/) {	
				
				my $indexInfo = 0;
				my $searchPattern = "";
				
				if ($species ne "") {
					$indexInfo = 7;
					$searchPattern = $species;
				}
				elsif ($assemblyTaxid ne "") { 
					$indexInfo = 6; 
					$searchPattern = $assemblyTaxid;
				}
				
				if (($tab[$indexInfo] =~ /^$searchPattern$/i) or ($kingdom ne "" && $searchPattern eq "")) { 
				
				
					my @gcfInfo = split(/\//, $tab[19]);  
					my $gcfName = pop(@gcfInfo);
					my $realDate = $tab[14];
					$realDate =~ s/\//-/g;
					
					my $genbankFile = $tab[19] . "/" . $gcfName . "_genomic.gbff.gz";
					my $dnaFile = $tab[19] . "/" . $gcfName . "_genomic.fna.gz";
					my $assemblyReport = $tab[19] . "/" . $gcfName . "_assembly_report.txt";
					
					if ($realDate gt $releaseDate) {
						
						$dnaFile = obtain_file($ftpServor, $dnaFile);
						$genbankFile = obtain_file($ftpServor, $genbankFile);
						$assemblyReport = obtain_file($ftpServor, $assemblyReport);
					
						download_file($ftpServor, $dnaFile);
						download_file($ftpServor, $genbankFile);
						download_file($ftpServor, $assemblyReport);
					
						# download sequences and check number of "N" characters
						my $fileFasta = $gcfName."_genomic.fna.gz";
						my $ucpFasta = $gcfName."_genomic.fna";
						if (-e $fileFasta) {
							gunzip $fileFasta => $ucpFasta or die "gunzip failed: $GunzipError\n";
							move($ucpFasta, $repositoryAssembly) or die "move failed: $!";
						}
				
						# download genome report
						my $fileReport =  $gcfName."_assembly_report.txt";
						if (-e $fileReport) {
							my $fileInformations = $gcfName."_informations.xls";
							move($fileReport, $repositoryAssembly) or die "move failed: $!";
						}
						
						# download genbank files
						my $fileGenbank = $gcfName."_genomic.gbff.gz";
						my $ucpGenbank = $gcfName."_genomic.gbff";
						if (-e $fileGenbank) {
							gunzip $fileGenbank => $ucpGenbank or die "gunzip failed: $GunzipError\n";
							move($ucpGenbank, $repositoryAssembly) or die "move failed: $!";
						}
				
						# association report and fasta
						my $fileFastaGenbank = $ucpFasta . "," . $ucpGenbank;
						$assemblyReportList{$fileReport} = $fileFastaGenbank;
					
						if ($quantity) { $quantity -= 1; }
						
					}
				}	
			}
			if (defined $quantity && $quantity == 0) {
				$quantity = undef;
				last;
			}
		}
		close(SUM) or die "close file error : $!";
		
		if (!keys %assemblyReportList) {
			print "##################################################################\n";
			print "Actual requirements are not found or are invalid for the database\n";
			print "It could be an error on the species name or the date\n";
			print "Check if the name is correctly write ex : Escherichia coli\n";
			print "The date could be too advances try again with an earlier date\n";
			print "Please try again with other requirements\n";
			print "##################################################################\n\n";
			
			if ($actualOS eq "linux") { unlink glob "*.dmp *.gz"  or die "for file *.dmp *.gz $!:"; }
			
			rmdir $kingdomRep or die "fail remove directory $!:";
			rmdir $repositoryAssembly or die "fail remove directory $!:"; 
			
			if ($components) {							
				rmdir $specificRep or die "fail remove directory $!:"; 
			} 
			else {
				rmdir $plasmidsRep or die "fail remove directory $!:";
				rmdir $chromosomesRep or die "fail remove directory $!:";
				rmdir $scaffoldsRep or die "fail remove directory $!:"; 
				rmdir $contigsRep or die "fail remove directory $!:";
			}
			exit(); 
		}
		
		# write summary files 
		my @keysList = keys %assemblyReportList;
		my $summary = "summary.xls";
		my $htmlSummary = "summary.html";
		my $specificSummary;
		my $plasmidsSummary;
		my $chromosomesSummary;
		my $scaffoldsSummary;
		my $contigsSummary;
		my @listComponents;
		
		if ($components) { 
			$specificSummary =  $components . "_summary.xls"; 
			push @listComponents, $components;
		} 
		else {
			$plasmidsSummary = "plasmids_summary.xls";
			$chromosomesSummary = "chromosomes_summary.xls";
			$scaffoldsSummary = "scaffolds_summary.xls";
			$contigsSummary = "contigs_summary.xls";
			
			push @listComponents, "plasmid";
			push @listComponents, "chromosome";
			push @listComponents, "scaffold";
			push @listComponents, "contig";
		}
		
		my $fileReport = ".".$fldSep. $repositoryAssembly . $fldSep . $keysList[0]; 
		
		my @summary_list = ($plasmidsSummary, $chromosomesSummary, $scaffoldsSummary, $contigsSummary);
		
		my @header = ();
		
		open(FILE, $fileReport) or die "error open file : $!";
		while(<FILE>) {
			chomp;
			
			$_ =~ s/^#*//;
			
			if($_ =~ /:/){
				my @ligne = split(':', $_);
				push(@header, $ligne[0]);
			}
		}
		close(FILE) or die "error close file : $!";
		
		open(HEAD, ">", $summary) or die " error open file : $!";
		foreach(@header) {
			print HEAD uc($_) . "\t";
		}
		
		print HEAD "PUBMED\tGC_PERCENT\tENTROPY\tSPECIES\tGENUS\tFAMILY\tORDER\tCLASS\t". 
								"PHYLUM\tKINGDOM\tCOUNTRY\tHOST\tISOLATION_SOURCE\tA_PERCENT\t".
									"T_PERCENT\tG_PERCENT\tC_PERCENT\n";
		close(HEAD) or die "error close file : $!";
		
		if ($components) {
			open(SUM, ">>", $specificSummary) or die "error open file : $!";
			print SUM "ID\tASSEMBLY\tDESCRIPTION\tLENGTH\tSTATUS\tLEVEL\t" . 
						"GC_PERCENT\tA_PERCENT\tT_PERCENT\tG_PERCENT\tC_PERCENT\n";	
			close(SUM) or die "error close file : $!";	
		} else {
			for my $sum(@summary_list) {
				open(SUM, ">>", $sum) or die "error open file : $!";
				print SUM "ID\tASSEMBLY\tDESCRIPTION\tLENGTH\tSTATUS\tLEVEL\t" . 
								"GC_PERCENT\tA_PERCENT\tT_PERCENT\tG_PERCENT\tC_PERCENT\n";	
				close(SUM) or die "error close file : $!"; 
			}
		}	
		
		# print "$actualOS\n";
		
		for my $file (@keysList) {
			my $file1 = $repositoryAssembly . $fldSep . $file;
			my @fastaGenbank = split(",", $assemblyReportList{$file});
			my $extFasta = $fastaGenbank[0];
			my $extGenbank = $fastaGenbank[1];
			my $file2 = $repositoryAssembly . $fldSep . $extFasta;
			my $file3 = $repositoryAssembly . $fldSep . $extGenbank;
			
			write_assembly($file1, $file2, $file3, $summary, $repositoryAssembly,
			$chromosomesSummary, $plasmidsSummary, $scaffoldsSummary,
			$contigsSummary, $specificSummary, $components, $kingdom,  $actualOS, @header);
		}
		
		write_html_summary($summary);
		
		my @listComponentFasta = create_component_sequence_file($fldSep, $repositoryAssembly, @listComponents);
		
		foreach my $componentFasta (@listComponentFasta) {
			if ($componentFasta =~ /plasmid/) { 
				move($componentFasta, $plasmidsRep) or die "move failed: $!"; 
			}
			elsif ($componentFasta =~ /chromosome/) {
				move($componentFasta, $chromosomesRep) or die "move failed: $!";
			}
			elsif ($componentFasta =~ /scaffold/) {
				move($componentFasta, $scaffoldsRep) or die "move failed: $!";
			}
			elsif ($componentFasta =~ /contig/) {
				move($componentFasta, $contigsRep) or die "move failed: $!";
			}
		}
		
		
		move($summary, $repositoryAssembly) or die "move failed: $!";
		move($htmlSummary, $repositoryAssembly) or die "move failed: $!";
		
		if ($components) {
			move($specificSummary, $specificRep) or die "move failed: $!";
			move($specificRep, $repositoryAssembly . $fldSep . $specificRep) or die "move failed: $!";
		}
		else {
			move($plasmidsSummary, $plasmidsRep) or die "move failed: $!";
			move($chromosomesSummary, $chromosomesRep) or die "move failed: $!";
			move($scaffoldsSummary, $scaffoldsRep) or die "move failed: $!";
			move($contigsSummary, $contigsRep) or die "move failed: $!";
			move($plasmidsRep, $repositoryAssembly . $fldSep . $plasmidsRep) or die "move failed: $!";
			move($chromosomesRep, $repositoryAssembly . $fldSep . $chromosomesRep) or die "move failed: $!";
			move($scaffoldsRep, $repositoryAssembly . $fldSep . $scaffoldsRep) or die "move failed: $!";
			move($contigsRep, $repositoryAssembly . $fldSep . $contigsRep) or die "move failed: $!";
		}
		move( $repositoryAssembly, $kingdomRep . $fldSep . $repositoryAssembly) or die "move failed: $!";
		
		if ($actualOS eq "linux") { unlink glob "*.dmp"  or die "for file *.dmp $!:"; }
		unlink glob "*.gz  *.dmp sequence.txt"  or die "$!: for file *.gz sequence.txt";

	} 
}


