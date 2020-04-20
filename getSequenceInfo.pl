#!/usr/bin/perl

use strict;
use warnings;

use Archive::Tar;
use Bio::SeqIO;
use Bio::Species;
use Date::Calc qw(:all);
use File::Copy qw(cp move);
use File::Path qw(rmtree);
use Net::FTP;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use LWP::Simple qw(get);
use POSIX qw(floor);

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

my $representation = "Complete Genome"; # assembly level of the genome

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
			get_assembly_summary_species($quantity, $releaseDate, $directory, $kingdom, $actualSpecies,
			$representation, $fldSep, $actualOS, $mainFolder, $assemblyTaxid);
		}
	}
	elsif ($assemblyTaxid ne "") {
		my @taxidList = split(/,/, $assemblyTaxid);
		
		foreach my $actualID (@taxidList) {
			get_assembly_summary_species($quantity, $releaseDate, $directory, $kingdom, $species,
			$representation, $fldSep, $actualOS, $mainFolder, $actualID);
		}
	} 
	else {
		get_assembly_summary_species($quantity, $releaseDate, $directory, $kingdom, $species,
		$representation, $fldSep, $actualOS, $mainFolder, $assemblyTaxid);
	}	
}


if ($sequenceID) { 
	my @sequenceIDList = split(/,/, $sequenceID);
	
	foreach my $enaID (@sequenceIDList) {
		sequence_ena($enaID);
	}
}


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
		
		-taxid specific taxid must be combin with -k option
		perl $0 -k "bacteria"  -taxid 9,24 -r "XXX" -get
		
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
	my ($quantity, $releaseDate, $directory, $kingdom, $species, $representation, $fldSep, $actualOS, $mainFolder, $assemblyTaxid) = @_;
	
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
	
	my $oldRep = "";
	
	my $kingdomRep;

	if (defined $mainFolder) {
		$kingdomRep = $mainFolder; 
	}
	else { 
		$kingdomRep = $kingdom."_".$start_year."_".$start_month."_".$start_day; 
	}
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
				my $regex = "";
				
				if ($species ne "") {
					$indexInfo = 7;
					$searchPattern = $species;
					$regex = qr/$searchPattern/i;
				}
				elsif ($assemblyTaxid ne "") { 
					$indexInfo = 6; 
					$searchPattern = $assemblyTaxid;
					$regex = qr/^$searchPattern$/i;
				}
				
				if (($tab[$indexInfo] =~ $regex) or ($kingdom ne "" && $searchPattern eq "")) { 
				
				
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
			print "$species $kingdom $assemblyTaxid\n";
			print "##################################################################\n\n";
			
			if ($actualOS eq "linux") { unlink glob "*.dmp *.gz"  or die "for file *.dmp *.gz $!:"; }
			
			my $emptyFolder = empty_folder($kingdomRep, $fldSep);
			if ($emptyFolder == 1) { rmdir $kingdomRep or die "fail remove directory $!:"; }
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
		}
		else {
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
			
			if ($components) {
				foreach my $componentFasta (@listComponentFasta) {
					move($componentFasta, $specificRep) or die "move failed: $!"; 
				}
			}
			else {
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
}
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
	
	my @header_report = keys %hashInformations;
	
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
	else {
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
		
		if ($desc[1]) { $level = $desc[1]; }
		
		if ($components) {
			if ($desc[0] =~ /$components/) {
				$status = $components;
				$info = $seqID . "\t" . $assembly_name ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t"
							. $gcpercent."\t". $a_percent ."\t". $t_percent ."\t". $g_percent ."\t". $c_percent . "\n";
				add_to_file($specific_summary, $info);
			}
		}	
		else {
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
	my($tax_id, $taxonomic_file) = @_;
	my $species = "na";
	my $genus = "na";
	my $family = "na";
	my $order = "na";
	my $class = "na";
	my $phylum = "na"; 
	
	# my ($species,$genus,$family,$order,$class,$phylum);
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
	my ($line, $htmlFile, $headers) = @_;
	
	my @assemblyHeader = split(/\t/, $headers);
	my @assemblyInfo = split(/\t/, $line);
	my %hashHeaderInfo;
	my $nbOfCell = 7;
	my $fullLine = floor(($#assemblyHeader + 1) / $nbOfCell);
	my $restCell = $#assemblyHeader + 1 - $fullLine * $nbOfCell;
	
	
	for (my $i = 0; $i < $#assemblyHeader + 1; $i++) {
		$hashHeaderInfo{trim($assemblyHeader[$i])} = $assemblyInfo[$i];
	}
	
	my @keysHeaderInfo = keys %hashHeaderInfo;
	my $cellIndex = 0;
	
	open(HTML, ">>", $htmlFile) or die "error open HTML summary $!";
	for (my $turn = 0; $turn < $fullLine; $turn++) {
	
		print HTML "   <tr>\n";
		for my $header (@assemblyHeader[$cellIndex..$cellIndex + $nbOfCell - 1]) {
			print  HTML "    <th>$header</th>\n";
		}
		print HTML "   </tr>\n";
		
		print HTML "   <tr>\n";
		for my $header (@assemblyHeader[$cellIndex..$cellIndex + $nbOfCell - 1]) {
			if ($header =~ /PUBMED/ && $hashHeaderInfo{$header} ne "na") {
				print HTML "   <td><a href=https://www.ncbi.nlm.nih.gov/pubmed/?term=".
				"$hashHeaderInfo{$header} target=\"_blank\">$hashHeaderInfo{trim($header)}</a></td>";
			}
			else {
				print  HTML "    <td>$hashHeaderInfo{trim($header)}</td>\n";
			}
		}
		print HTML "   </tr>\n";
		
		$cellIndex +=  $nbOfCell;
	}
	
	print HTML "   <tr>\n";
	for my $header(@assemblyHeader[$cellIndex..$#keysHeaderInfo]) {
		print  HTML "    <th>$header</th>\n";
	}
	print HTML "   <tr>\n";
	
	print HTML "   <tr>\n";
	for my $header(@assemblyHeader[$cellIndex..$#keysHeaderInfo]) {
		print  HTML "    <td>$hashHeaderInfo{trim($header)}</td>\n";		
	}
	print HTML "   <tr>\n";
	
	print HTML "  </table>\n";
	close(HTML) or die "error close HTML summary $!";
}
#------------------------------------------------------------------------------
#getTaxonomicRanks (function allowing to get taxonomic ranks from Genbank file)
sub get_taxonomic_rank_genbank {
	my ($genbank) = @_;
	my $species = "na";
	my $genus = "na";
	my $family = "na";
	my $order = "na";
	my $class = "na";
	my $phylum = "na"; 
	my $kingdomGB = "na";

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
#------------------------------------------------------------------------------
#add all sequences components to file
sub create_component_sequence_file {
	my ($fldSep, $repository, @listComponent) = @_;
	
	my @listFnaFile;
	
	opendir(my $dh, $repository) || die "Can't opendir $repository: $!";
	@listFnaFile = grep{/fna$/} readdir($dh);
	closedir $dh;
	
	my @listComponentFasta;

	foreach my $component (@listComponent) {
	
		my $componentFasta = $component.".fasta";
		
		foreach my $fnaFile (@listFnaFile) {
		
			# my $actualFile = $repository . $fldSep . $fnaFile;
			
			my $seq;
			my $seqIO = Bio::SeqIO->new(-format=>'Fasta', -file=>$repository . $fldSep . $fnaFile);
			
			while ($seq = $seqIO->next_seq()) {
			
				my $seqDesc = $seq->desc;
				
				if ($seqDesc =~ /$component/) {
					my $seqID = $seq->id;
					my $seqNuc = $seq->seq;
					my $shift = 60;
					my @seqArray = split //, $seqNuc;
					my $newSeqNuc = "";
					
					if (length $seqNuc <= $shift) {
						$newSeqNuc = $seqNuc;
					}
					else {
						for(my $i = 0; $i < $#seqArray + 1; $i ++) {
							$newSeqNuc .=  $seqArray[$i];
							if (($i + 1) % $shift == 0) { $newSeqNuc .= ","; }
						}
					}
					
					open(FASTA, ">>", $componentFasta) or die "error open file $!:";
					print FASTA ">$seqID $seqDesc\n";
					foreach my $subSeqNuc (split /,/, $newSeqNuc) {
						print FASTA "$subSeqNuc\n";
					}
					close(FASTA) or die "error close file $!:";
				}
			}
		}
		if (-e $componentFasta) { push @listComponentFasta, $componentFasta; }
	}
	return @listComponentFasta;
}
# remove back and front spaces
sub trim {
my ($string) = @_;
$string =~ s/^\s+//;
$string =~ s/\s+$//;
return $string;
}
#------------------------------------------------------------------------------
# use trim in array
sub trim_array {
	my (@array) = @_;
	foreach my $value (@array) {
		$value = trim($value);
	}
	return @array;
}
#------------------------------------------------------------------------------
# check if folder is empty
sub empty_folder {
	my ($folder, $fldSep) = @_;
	!<$folder.$fldSep.*>;
}
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
  return (log($n) / log(2));
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
  
	if (! $hashFlank{'G'}) { $hashFlank{'G'} = 0;}
	if (! $hashFlank{'C'}) { $hashFlank{'C'} = 0;}

	if(length($seq) == 0) {
		return 0;
	}
	else {
		return (($hashFlank{'G'} + $hashFlank{'C'}) / (length($seq))) * 100;
	}
	
}
#------------------------------------------------------------------------------
# download file from ftp protocol
sub download_file {
	my($servor, $file) = @_;
	
	my $ftp = Net::FTP->new($servor, Debug => 0)
	or die "Cannot connect to $servor";
	
	$ftp->login("anonymous", "-anonymous@")
		or die "Cannot login ", $ftp->message;
	$ftp->binary;	
	$ftp->get($file) or die "get failed ", $ftp->message;
	
	$ftp->quit;
}
#------------------------------------------------------------------------------
# obtain file directory
sub obtain_file {
	my($servor, $link) = @_;
	if ($link =~ /$servor(.*)/) { return ($1); }
}



