#!/usr/bin/perl

################################################################################
## "Copyright 2019 Vincent Moco and David Couvin"
## licence GPL-3.0-or-later
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>. 
################################################################################

use strict;
use warnings;

my $version = "1.0.1"; # version
# Date and time of the current day (Beginning)
my ($start_year,$start_month,$start_day, $start_hour,$start_min,$start_sec) = Today_and_Now();

print "##################################################################\n";
print "## ---> Welcome to $0 (version $version)!\n";
print "## Start Date (yyyy-mm-dd, hh:min:sec): $start_year-$start_month-$start_day, $start_hour:$start_min:$start_sec\n";
print "##################################################################\n\n";


my @modules = qw(
	BioPerl
	Archive::Tar
	Bio::SeqIO
	Bio::Species
	Date::Calc
	File::Copy
	File::Path
	Net::FTP
	IO::Uncompress::Gunzip
	LWP::Simple
	POSIX
	Tk
	Tk::ProgressBar
	utf8
	Shannon::Entropy
);

foreach my $module (@modules) {
	if (isModuleInstalled($module)) {
	  print "$module is.................installed!\n";
	} else {
	  print "$module was not installed.\nLet us install it\n";
	  system("cpan -i -f $module");
	}
}
print "\n";

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

####################################################################################################
##  A Perl script allowing to get sequence information from GenBank, RefSeq or ENA repositories.
##  
####################################################################################################

### main program
### parameters
my $directory = "genbank";
	
my $kingdom = ""; # kingdom of organism

my $releaseDate = "0000-00-00"; # sequence are download from this release date

my $components = "plasmid,chromosome,scaffold,contig";  # components of the assembly

my $species = ""; # species name

my $getSummary; # indicates if a new assembly report is required

my $assemblyLevel = "Complete Genome,Chromosome,Scaffold,Contig"; # assembly level of the genome

my $quantity; # number of sequences to download

my $sequenceID;

my $ftpServor = "ftp.ncbi.nlm.nih.gov";

my $enaFtpServor = "ftp.sra.ebi.ac.uk";

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

my $sraID; # SRA  sequence ID

my $assemblyPrjID; # assembly or prj ID

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

##requirements
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
	elsif ($ARGV[$i]=~/-level/i or $ARGV[$i]=~/-l/i) {
		$assemblyLevel = $ARGV[$i+1];
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
	elsif ($ARGV[$i]=~/-fastq/i) {
		$sraID = $ARGV[$i+1];
	}
	elsif ($ARGV[$i]=~/-assembly_or_project/i) {
		$assemblyPrjID = $ARGV[$i+1];
	}
}

#define folcer separator and OS
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
	
	my @levelList = split /,/, $assemblyLevel;
	
	if ($species ne "") {
		my @speciesList = split(/,/, $species);
		
		foreach my $actualSpecies (@speciesList) {
			get_assembly_summary_species($quantity, $releaseDate, $directory, $kingdom,
			$actualSpecies,\@levelList, $fldSep, $actualOS, $mainFolder, $assemblyTaxid);
		}
	}
	elsif ($assemblyTaxid ne "") {
		my @taxidList = split(/,/, $assemblyTaxid);
		
		foreach my $actualID (@taxidList) {
			get_assembly_summary_species($quantity, $releaseDate, $directory, $kingdom, 
			$species,\@levelList, $fldSep, $actualOS, $mainFolder, $actualID);
		}
	} 
	else {
		get_assembly_summary_species($quantity, $releaseDate, $directory, $kingdom, 
		$species,\@levelList, $fldSep, $actualOS, $mainFolder, $assemblyTaxid);
	}	
}


if ($sequenceID) { 
	my @sequenceIDList = split /,/, $sequenceID;
	
	foreach my $enaID (@sequenceIDList) {
		sequence_ena($enaID);
	}
}

if ($sraID) {
	my @sraIDList = split /,/, $sraID;
	
	foreach my $sra (@sraIDList) {
		download_ena_fastq($enaFtpServor, $sra);
	}	
}

if ($assemblyPrjID) {
	download_assembly_or_project($assemblyPrjID, $ftpServor, $fldSep, $directory);
}

my ($end_year,$end_month,$end_day, $end_hour,$end_min,$end_sec) = Today_and_Now();

my ($D_y,$D_m,$D_d, $Dh,$Dm,$Ds) =
      Delta_YMDHMS($start_year,$start_month,$start_day, $start_hour, $start_min, $start_sec,
                   $end_year, $end_month, $end_day, $end_hour,$end_min,$end_sec);

print "End Date: $start_year-$start_month-$start_day, $start_hour:$start_min:$start_sec\n";
print "Thank you for using getSequenceInfo!\n";
print "Execution time: $D_y years, $D_m months, $D_d days, $Dh:$Dm:$Ds (hours:minutes:seconds)\n";

### subroutine 
# display global help document
sub help_user_simple {
	my $programme = shift @_;
	print STDERR  "Usage : perl $programme -k XXX -s \"XXX\"  -r \"XXX\" -date yyyy-mm-dd -get  \n";
	print "type perl $programme -version or perl $programme -v to get actual version\n";
	print "type perl $programme -help or perl $programme -h to get full help\n";
}
#------------------------------------------------------------------------------
# display full help document
sub help_user_advance {
	print <<HEREDOC;
	
	Name: 
		$0
	
	Synopsis:
		A Perl script allowing to get sequence information from GenBank RefSeq or ENA repositories.
		
	Usage:
	  perl $0 [options]
	  examples: 
	     perl $0 -k bacteria -s "Helicobacter pylori" -l "Complete Genome" -date 2019-06-01 -get 
	     perl $0 -k viruses -q 5 -l "Complete Genome" -date 2019-06-01 -get
	     perl $0 -k "bacteria" -taxid 9,24 -get -q 10 -c plasmid -dir genbank -o Results
	     perl $0 -ena BN000065
	     perl $0 -fastq ERR818002
						 	
	Kingdoms:
		archaea
		bacteria
		fungi
		invertebrate
		plant
		protozoa
		vertebrate_mammalian
		vertebrate_other
		viral
	
	Assembly levels:
		"Complete Genome"
		Chromosome
		Scaffold
		Contig 
	
	General:
		-help or -h			displays this help 	
		-version or -v			displays the current version of the program
		
	Options:
		-get 				allows to obtain a new assembly summary	
		-k or -kingdom [XXX]		allows to indicate kingdom of the organism (see the examples above)
		-s or -species [XXX]		allows to indicate the species (must be combined with -k option)
		-taxid [XXX]			allows to indicate a specific taxid (must be combined with -k option)
		-assembly_or_project [XXX]	allows to indicate a specific assembly accession or bioproject (must be combined with -k option)
		-date [XXX]			indicates the release date (with format yyyy-mm-dd) from which sequence information are available
		-l or -level [XXX]		allows to select a specific assembly level (e.g. "Complete Genome")
		-o or output [XXX]		allows users to name the output result folder
		-q or -quantity [XXX]		allows to limit the total number of assemblies to be downloaded
		-c or -components [XXX]		allows to select specific components of the assembly (e.g. plasmid, chromosome, ...)
		-ena [XXX] 			allows to download report and fasta file given a ENA sequence ID 
		-fastq [XXX]			allows to download FASTQ sequences from ENA given a run accession (https://ena-docs.readthedocs.io/en/latest/faq/archive-generated-files.html)
		
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
	my ($quantity, $releaseDate, $directory, $kingdom, $species, $levelList, $fldSep, $actualOS, $mainFolder, $assemblyTaxid) = @_;
	
	# assembly_summary.txt file from NCBI FTP site
	my $assemblySummary = "/genomes/$directory/$kingdom/assembly_summary.txt"; 
	
	# lineage folder
	my $lineage_file = "/pub/taxonomy/new_taxdump/new_taxdump.tar.gz";
	
	# allow to check old summary download
	my $oldKingdom = ""; 
	
	# check assembly summary download
	if ($getSummary || ! -e $assemblySummary) { 
		download_file($ftpServor, $assemblySummary);
		open(KIN, ">", "kingdom.txt") or die "error open file $!:";
		print KIN $kingdom;
		close(KIN) or die "error close file $!:";
		
	} elsif ($kingdom ne "" && $species eq  "") {
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
	
	if ($actualOS eq "linux") {	
		# initialiaze tar manipulation
		my $tar = Archive::Tar->new;
	
		# download taxdump folder
		download_file($ftpServor, $lineage_file);
		$tar->read("new_taxdump.tar.gz");
		$tar->extract_file("rankedlineage.dmp");
	}
	
	my $oldAssemblyRep;
	my $kingdomRep;

	if ($mainFolder) {
		$kingdomRep = $mainFolder; 
	}
	else { 
		$kingdomRep = $kingdom."_".$start_year."_".$start_month."_".$start_day; 
	}
	
	
	mkdir $kingdomRep unless -d $kingdomRep;
	
	# Repository for request
	my $repositoryAssembly = "assembly_repository_".$species."_".$assemblyTaxid."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
	mkdir $repositoryAssembly unless -d $repositoryAssembly;
	
	$oldAssemblyRep = "." . $fldSep . $kingdomRep . $fldSep . $repositoryAssembly;
	if (-d $oldAssemblyRep) { rmtree($oldAssemblyRep) }
	
	#  Repository for fna file
	my $repositoryFNA = "Assembly";
	mkdir $repositoryFNA unless -e $repositoryFNA;
	
	# Repository for genbank file
	my $repositoryGenbank = "GenBank";
	mkdir $repositoryGenbank unless -e $repositoryGenbank;
	
	# Reposotiry for report file
	my $repositoryReport = "Report";
	mkdir $repositoryReport  unless -e $repositoryReport ;
	
	# Repositories for required components
	my %componentsRepHash;
	
	for my $component (split /,/, $components) {
		my $specificRep = $component."_".$species."_".$kingdom."_".$start_year."_".$start_month."_".$start_day;
		mkdir $specificRep unless -d $specificRep;
		$componentsRepHash{$component} = $specificRep;
	}
	
	my %assemblyReportList;
	my @assemblyRepresentationList = qw/Full Partial/;
	my @levelList = @{$levelList};
	my $checkCompleteGenome = grep(/complete genome/i, @levelList);
	
	if ($checkCompleteGenome > 0) {@assemblyRepresentationList = qw/Full/;}
	
	if (-e "assembly_summary.txt") {
		# Read file
		open (SUM, "assembly_summary.txt") or die "open assembly_summary.txt : $!";
		while(<SUM>) {
			chomp;
			if ($_ !~  m/^#/) {
				my @infoList = split /\t/, $_;
				my $foundAssemRep = grep (/$infoList[13]/i, @assemblyRepresentationList);
				my $checkLevel = grep(/$infoList[11]/i, @levelList);
				
				if ($foundAssemRep > 0 && $checkLevel > 0) {
					my $indexInfo = 0;
					my $searchPattern = "";
					my $regex = "";
					
					if ($species !~ //) {
						$indexInfo = 7;
						$searchPattern = $species;
						$regex = qr/$searchPattern/i;
					}
					elsif ($assemblyTaxid !~ //) { 
						$indexInfo = 6; 
						$searchPattern = $assemblyTaxid;
						$regex = qr/^$searchPattern$/i;
					}
					
					if (($infoList[$indexInfo] =~ $regex) or ($kingdom !~ // && $searchPattern =~ //)) { 
						my @gcfInfo = split(/\//, $infoList[19]);  
						my $gcfName = pop(@gcfInfo);
						my $realDate = $infoList[14];
						$realDate =~ s/\//-/g;
						
						my $genbankFile = $infoList[19] . "/" . $gcfName . "_genomic.gbff.gz";
						my $dnaFile = $infoList[19] . "/" . $gcfName . "_genomic.fna.gz";
						my $assemblyReport = $infoList[19] . "/" . $gcfName . "_assembly_report.txt";
						
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
								move($ucpFasta, $repositoryFNA) or die "move failed: $!";
							}
					
							# download genome report
							my $fileReport =  $gcfName."_assembly_report.txt";
							if (-e $fileReport) {
								my $fileInformations = $gcfName."_informations.xls";
								move($fileReport, $repositoryReport) or die "move failed: $!";
							}
							
							# download genbank files
							my $fileGenbank = $gcfName."_genomic.gbff.gz";
							my $ucpGenbank = $gcfName."_genomic.gbff";
							if (-e $fileGenbank) {
								gunzip $fileGenbank => $ucpGenbank or die "gunzip failed: $GunzipError\n";
								move($ucpGenbank, $repositoryGenbank) or die "move failed: $!";
							}
					
							# association report and fasta
							my $fileFastaGenbank = $ucpFasta . "," . $ucpGenbank;
							$assemblyReportList{$fileReport} = $fileFastaGenbank;
						
							if ($quantity) { $quantity -= 1; }
							
						}
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
			print "$species $kingdom $assemblyTaxid @levelList\n";
			print "##################################################################\n\n";
			
			if ($actualOS eq "linux") { unlink glob "*.dmp *.gz"  or die "for file *.dmp *.gz $!:"; }
			
			if (empty_folder($kingdomRep)) { rmdir $kingdomRep or die "fail remove directory $!:"; }
			rmdir $repositoryAssembly or die "fail remove directory $!:";
			rmdir $repositoryFNA or die "fail remove directory $!:";
			rmdir $repositoryGenbank or die "fail remove directory $!:"; 
			rmdir $repositoryReport or die "fail remove directory $!:"; 
			
			for my  $componentRep (values %componentsRepHash) {
				rmdir $componentRep or die "fail remove directory $!:"; 
			}
			
		} 
		else {
			# write summary files 
			my %componentsSumHash;
			my @keysList = keys %assemblyReportList;
			my $summary = "summary.xls";
			my $htmlSummary = "summary.html";
			
			for my $component (split /,/, $components) {
				my $specificSummary =  $component. "_summary.xls"; 
				$componentsSumHash{$component} = $specificSummary;
			}
		
			my $fileReport = "." . $fldSep. $repositoryReport . $fldSep . $keysList[0];
			my @header = ();
			
			open(FILE, $fileReport) or die "error open file : $!";
			while(<FILE>) {
				chomp;
				if($_ =~ /:/){
					$_ =~ s/^#*//;
					my @ligne = split /:/, $_;
					push(@header, $ligne[0]);
				}
			}
			close(FILE) or die "error close file : $!";
			
			open(HEAD, ">", $summary) or die " error open file : $!";
			foreach(@header) {
				print HEAD $_ . "\t";
			}
			
			print HEAD "Pubmed\tNucle score\tSpecies\tGenus\tFamily\tOrder\tClass\t"; 
			print HEAD "Phylum\tKingdom\tCountry\tHost\tIsolation source\tA percent\t";
			print HEAD "T percent\tG percent\tC percent\tN percent\tGC percent\t";
			print HEAD "ATGC ratio\tLength\tShape\n";
			close(HEAD) or die "error close file : $!";
			
	
			foreach my $componentSummary (values %componentsSumHash) {
				open(SUM, ">>", $componentSummary) or die "error open file : $!";
				print SUM "Id\tAssembly\tDescription\tLength\tStatus\tLevel\t";
				print SUM "GC percent\tA percent\tT percent\tG percent\tC percent\n";	
				close(SUM) or die "error close file : $!";
			}
		 
		 
			for my $file (@keysList) {
				my $reportFile =  $repositoryReport . $fldSep . $file;
				my @fastaGenbank = split /,/, $assemblyReportList{$file};
				my $extFasta = $fastaGenbank[0];
				my $extGenbank = $fastaGenbank[1];
				my $fnaFile = $repositoryFNA . $fldSep . $extFasta;
				my $genbankFile = $repositoryGenbank . $fldSep . $extGenbank;
				
				write_assembly($reportFile, $fnaFile, $genbankFile, $summary, $repositoryAssembly,\%componentsSumHash, $kingdom,  $actualOS, \@header);
			}
			
			write_html_summary($summary);
			
			my @componentList = keys %componentsSumHash;
			my %componentFastaHash = create_component_sequence_file($fldSep, $repositoryFNA, \@componentList);
			
			foreach my $component (keys %componentFastaHash) { 
				move($componentFastaHash{$component}, $componentsRepHash{$component}) or die "move failed: $!"; 
			}
			
			move($summary, $repositoryAssembly) or die "move failed: $!";
			move($htmlSummary, $repositoryAssembly) or die "move failed: $!";
			move($repositoryFNA, $repositoryAssembly . $fldSep . $repositoryFNA) or die "move failed: $!";
			move($repositoryGenbank, $repositoryAssembly . $fldSep . $repositoryGenbank) or die "move failed: $!";
			move($repositoryReport  , $repositoryAssembly . $fldSep . $repositoryReport) or die "move failed: $!";
			
			for my $component (keys %componentsSumHash) {
				move($componentsSumHash{$component}, $componentsRepHash{$component}) or die "move failed: $!";
				move($componentsRepHash{$component}, $repositoryAssembly . $fldSep . $componentsRepHash{$component}) or die "move failed: $!"
			}
			move($repositoryAssembly, $kingdomRep . $fldSep . $repositoryAssembly) or die "move failed: $!";
			
			if ($actualOS eq "linux") { unlink glob "*.dmp"  or die "for file *.dmp $!:"; }
			unlink glob "*.gz  *.dmp sequence.txt"  or die "$!: for file *.gz sequence.txt";
		}
	} 
}
#write general assembly file
sub write_assembly {
	my ($reportFile, $fnaFile, $genbankFile, $summary, $repositoryAssembly,
	$componentsSumHashRef, $kingdom,  $actualOS, $headerRef) = @_;
	
	my %componentsSumHash = %{$componentsSumHashRef};
	my @header = @{$headerRef};
	my %hashInformations = ();
	my $seq = "";
	my $genomeName = "";
	my $country = "na";
	my $GCpercent = -1;
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
	my $shape = "na";
	my $geneNumber = "na";
	
	open(REP, "<", $reportFile) or die "error open file $!:";
	while (<REP>) {
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
	close(REP) or die "error close file $!:";
	
	my @header_report = keys %hashInformations;
	
	open(SUM, ">>", $summary) or die "error open file $!:";
	foreach my $k(@header) {
		if (grep $_ eq $k, @header_report) {
			my $information = $hashInformations{$k};
			
			if ($k =~ /Assembly name/) { $genomeName = $information; }
			
			if (($information =~ /^\s*$/) || ($information eq "")) {
				print SUM "na\t";
			} else { 
				print SUM $information . "\t";
			}
		} else {
			print SUM "na\t";
		}
	}
	close(SUM) or die "error close file $!:";
	
	open(FNA, "<", $fnaFile) or die "Could not open $!:";
	while (<FNA>) {
		chomp;
		if ($_ !~ /^>/) { $seq .= $_; }
	}
	close(FNA)  or die "error close file :$!";
	
	if ($hashInformations{' Taxid'} !~ /\s+/) { $taxId = $hashInformations{' Taxid'} };

	if ($actualOS eq "linux") {
		($species, $genus, $family, $order, $class, $phylum) =  get_taxonomic_rank($taxId, "rankedlineage.dmp");
	}
	else {
		($species, $genus, $family, $order, $class, $phylum) = get_taxonomic_rank_genbank($genbankFile);
	}
	
	$GCpercent = gc_percent($seq);
	my ($ade, $thy, $gua, $cyt, $n, $length) = number_nuc_length_seq($fnaFile);
	my ($aPercent, $tPercent, $gPercent, $cPercent, $nPercent) = nucleotid_percent($ade, $thy, $gua, $cyt, $n, $length);
	my $atgcRatio =	 atgc_ratio($ade, $thy, $gua, $cyt);
	
	my @percentList = ($aPercent, $tPercent, $gPercent, $cPercent, $nPercent);
	
	my $variance = shift_data_variance(@percentList);
	my $nucleScore = nucle_score($variance, $GCpercent, $atgcRatio, $length);
	
	open(GBFF, "<", $genbankFile) or die "Error open file $!:";
	while(<GBFF>) {
		chomp;
		if ($_ =~ /\/country="(.*)"/) { $country = trim($1); }
		if ($_ =~ /PUBMED(.*)/) {  $pubmedId = trim($1); }
		if ($_ =~ /\/host="(.*)"/) {  $host = trim($1); }		
		if ($_ =~ /\/isolation_source="(.*)"/) {  $isoSource = trim($1); }
		if ($_ =~ /\(Genes \(total\)\s+::(.*)/) { $geneNumber = trim($1); }
		if ($_ =~ /LOCUS.*\s+([a-z]{1,})\s+[a-z]{1,}\s+[0-9]{2,}-[a-z]{1,}-[0-9]{4,}$/i) { $shape = trim($1); }
	}
	close(GBFF) or die "error close file $!:";

	
	open(SUM, ">>", $summary) or die "error open file $!:";
	print SUM $pubmedId . "\t" . $nucleScore . "\t" . $species . "\t" . $genus . "\t" . $family ."\t" ; 
	print SUM $order . "\t" . $class . "\t" . $phylum . "\t" . $kingdom . "\t" . $country . "\t" . $host . "\t"; 
	print SUM $isoSource  . "\t" . $aPercent . "\t" . $tPercent . "\t" . $gPercent . "\t" . $cPercent  ."\t" ;
	print SUM $nPercent . "\t" . $GCpercent ."\t". $atgcRatio ."\t". $length . "\t". $shape."\n"; 
	close(SUM) or die "error close file: $!";
	
	write_assembly_component($fnaFile, $genomeName, \%componentsSumHash);
}
#------------------------------------------------------------------------------
# get assembly component
sub write_assembly_component {
	my($fnaFile, $genomeName, $componentsSumHashRef) = @_;
	
	my %componentsSumHash = %{$componentsSumHashRef};
	my $status = "na";
	my $level = "na";
	my $gcpercent;
	my $tmp_fasta_file = "sequence.txt";
	my @desc = ();
	
	# check each sequence from (multi-)fasta file
	my ($seq, $inputfile);
	#my @tabNcounts = ();

	# extract sequence details
	my $seqIO = Bio::SeqIO->new(-format=>'Fasta', -file=>$fnaFile);
	
	while ($seq = $seqIO->next_seq()) {
		my $seqID = $seq->id; # ID of sequence
		my $seqDesc = $seq->desc; # Description of sequence
		my $globalSeq = $seq->seq;
		my $seqLength = $seq->length();
		
		open(TSEQ, ">", $tmp_fasta_file) or die "Error open file: $!";
		print TSEQ $globalSeq;
		close(TSEQ);
		
		my ($ade, $thy, $gua, $cyt, $n, $length) = number_nuc_length_seq($tmp_fasta_file);
		
		my ($aPercent, $tPercent, $gPercent, $cPercent, $nPercent) = nucleotid_percent($ade, $thy, $gua, $cyt, $n, $length);
		
		$gcpercent = gc_percent($globalSeq);
		
		@desc = split /,/, $seqDesc;
		
		if ($desc[1]) { $level = $desc[1]; }
		
		foreach my $component (keys %componentsSumHash) {
			if ($desc[0] =~ /$component/) {
				$status = $component;
				my $info = $seqID . "\t" . $genomeName ."\t" . $seqDesc . "\t" . $seqLength . "\t" . $status . "\t" . $level ."\t";
				$info.= $gcpercent."\t". $aPercent ."\t". $tPercent ."\t". $gPercent ."\t". $cPercent . "\n";
				add_to_file($componentsSumHash{$component}, $info);
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
	
	unlink "fichier.txt" or die "error delete file :$!";
	
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
	open(FILE, ">>", $report_file) or die "could not open: $!";
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
			if ($header =~ /PUBMED/i && $hashHeaderInfo{$header} ne "na") {
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
	my ($fldSep, $repository, $listComponentRef) = @_;
	
	my @listFnaFile;
	my @listComponent = @{$listComponentRef};
	
	opendir(my $dh, $repository) || die "Can't opendir $repository: $!";
	@listFnaFile = grep{/fna$/} readdir($dh);
	closedir $dh;
	
	my %componentFastaHash;

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
		if (-e $componentFasta) { $componentFastaHash{$component} = $componentFasta; }
	}
	return %componentFastaHash;
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
	my $dirname = shift;
    opendir(my $dholder, $dirname) or die "error not a directory";
	my $isEmpty = scalar(grep { $_ ne "." && $_ ne ".." } readdir($dholder));
	if ($isEmpty == 0) { return $isEmpty; }
}
#------------------------------------------------------------------------------
# number nucleotid and length
sub number_nuc_length_seq {
	my ($fastaFile) = @_;
	my $ade = 0;
	my $thy = 0;
	my $gua = 0;
	my $cyt = 0;
	my $n = 0;
	my $length = 0;
	
	open (FASTA, "<", $fastaFile) or die "Could not open $!";
	while (<FASTA>) {
		chomp;
		if ($_ !~ />/) {
			my @seq = split //, $_;
			
			for my $nuc (@seq) {
				$length +=1 ;
				if ($nuc =~ /a/i) {$ade+=1;}
				elsif ($nuc =~ /t/i) {$thy+=1;}
				elsif ($nuc =~ /g/i) {$gua+=1;}
				elsif ($nuc =~ /c/i) {$cyt+=1;}
				elsif ($nuc =~ /n/i) {$n+=1;}
			}
		}
	}
	close(FASTA) or die "Error close file :$!";
	return ($ade, $thy, $gua, $cyt, $n, $length);
	
}
#------------------------------------------------------------------------------
# compute percentage of nucleotid
sub nucleotid_percent {
	my($ade, $thy, $gua, $cyt, $n, $length) = @_;
	
	my $adePercent = $ade / $length * 100;
	my $thyPercent = $thy / $length * 100;
	my $guaPercent = $gua / $length * 100;
	my $cytPercent = $cyt / $length * 100;
	my $nPercent = $n / $length * 100;
	
	return ($adePercent, $thyPercent, $guaPercent, $cytPercent, $nPercent);
 
}
#------------------------------------------------------------------------------
# compute ATGC ratio 
sub atgc_ratio {
	my ($ade, $thy, $gua, $cyt) = @_;
	
	return (($ade + $thy) / ($gua + $cyt));
	
}
#------------------------------------------------------------------------------
# variance
sub shift_data_variance {
	my (@data) = @_;
	
	if ($#data + 1 < 2) { return 0.0; }

	my $K = $data[0];
	my ($n, $Ex, $Ex2) = 0.0;
	
	for my $x (@data) {
		$n = $n + 1;
		$Ex += $x - $K;
		$Ex2 += ($x - $K) * ($x - $K);
	}
	
	my $variance = ($Ex2 - ($Ex * $Ex) / $n) / ($n); ## ($n - 1)
	
	return $variance;

}
#------------------------------------------------------------------------------
# nucle score
sub nucle_score {
	my ($variance, $gcPercent, $atgcRatio, $length) = @_;
	
	return (($variance * $gcPercent * $atgcRatio) / $length);
}
#------------------------------------------------------------------------------
sub log2 {
  my $n = shift;
  return (log($n) / log(2));
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
#------------------------------------------------------------------------------
# download fastq file from ENA
sub download_ena_fastq {
	my ($enaFtpServor, $sraId) = @_;
	
	my $fastqDir = "/vol1/fastq/";
	my $dir1 = substr $sraId, 0, 6;
	my $dir2 = "000";
	my $digits = substr $sraId, 3;
	my $fastqRep =  $sraId . "_folder";
	

	if (length $digits == 6) {
		$dir2 = $sraId;
		$fastqDir .= $dir1 . "/" . $dir2 . "/";
	}
	elsif (length $digits  > 6) {
		my $digitsNumber = 0;
		my @digitsList = split //, (substr $digits, 6);
	
		foreach my $char (@digitsList) {
			if (length $dir2 >= 1) {
				$dir2 = substr $dir2, 0, (length $dir2) - 1;
				$digitsNumber += 1;
			}
		}
		$dir2 .= substr $digits, -$digitsNumber;
		$fastqDir .= $dir1 . "/" . $dir2 . "/" . $sraId . "/";
	}
	
	my $ftp = Net::FTP->new($enaFtpServor, Debug => 0)
	or die "Cannot connect to $enaFtpServor";
	
	$ftp->login("anonymous", "-anonymous@")
		or die "Cannot login ", $ftp->message;
	$ftp->binary;
	
	$ftp->cwd($fastqDir) 
		or die "maybe undefined sequence id, can't go to $fastqDir: ", $ftp->message;
	
	my @fastqFiles = $ftp->ls("$sraId*"); 
	
	if (!grep(/^$/, @fastqFiles)) {
		
		if (-d $fastqRep) { rmtree($fastqRep) }
		mkdir $fastqRep;
		
		foreach my $fastqFile (@fastqFiles) {
			$ftp->get($fastqFile) or die "get failed ", $ftp->message;
		
			my @baseAndExt = split /\./, $fastqFile;
			#my $unzipFastq = $baseAndExt[0] . ".fastq";
		
			#gunzip $fastqFile => $unzipFastq or die "gunzip failed: $GunzipError\n";
			move($fastqFile, $fastqRep) or die "move failed: $!"; # DC replaced $unzipFastq by $fastqFile
		} 
		#unlink glob "*fastq.gz"  or die "$!: for file *fastq.gz";
	}
	$ftp->quit;
}
#------------------------------------------------------------------------------
# download fastq file from ENA
sub get_assembly_or_project {
	my ($file, $sequence, $ftpServor, $fldSep) = @_;
	
	my $pattern;
	my $indexInfo;
	my %folderHash;
	
	#  Repository for fna file
	my $repositoryFNA = "Assembly";
	
	# Repository for genbank file
	my $repositoryGenbank = "GenBank";
	
	# Reposotiry for report file
	my $repositoryReport = "Report";
	
	# global repository
	my $repositorySequence = $sequence;

	
	if ($sequence =~ /^GC[AF]_(.*)/) {
		$indexInfo = 0;
		$pattern = $1;
	}
	elsif ($sequence =~ /^PRJ/) {
		$indexInfo = 1;
		$pattern = $sequence;
	}
	
	open(SUM, $file) or die "error open file $!:";
	while(<SUM>) {
		chomp;
		if ($_ !~ /^#/) {
			my @infoList = split /\t/, $_;
			if ($infoList[$indexInfo] =~ /$pattern/) {
				my @gcfInfo = split(/\//, $infoList[19]);  
				my $gcfName = pop(@gcfInfo);
				
				
				my $genbankFile = $infoList[19] . "/" . $gcfName . "_genomic.gbff.gz";
				my $dnaFile = $infoList[19] . "/" . $gcfName . "_genomic.fna.gz";
				my $assemblyReport = $infoList[19] . "/" . $gcfName . "_assembly_report.txt";
				
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
					$folderHash{$ucpFasta} = $repositoryFNA;
				}
				
				# download genome report
				my $fileReport =  $gcfName."_assembly_report.txt";
				if (-e $fileReport) {
					$folderHash{$fileReport} = $repositoryReport;
				}
				
				# download genbank files
				my $fileGenbank = $gcfName."_genomic.gbff.gz";
				my $ucpGenbank = $gcfName."_genomic.gbff";
				if (-e $fileGenbank) {
					gunzip $fileGenbank => $ucpGenbank or die "gunzip failed: $GunzipError\n";
					$folderHash{$ucpGenbank} = $repositoryGenbank;
				}
				
				
			}
		}
	}
	close(SUM) or die "error close file $!";	
	
	if (keys %folderHash) {
		if (-e $repositorySequence) {rmtree($repositorySequence);}
		
		mkdir $repositorySequence;
		mkdir $repositoryFNA;
		mkdir $repositoryGenbank;
		mkdir $repositoryReport;
		
		for my $ucpFile (keys %folderHash) {
			move($ucpFile, $folderHash{$ucpFile}) or die "error move file $!:";
		}
		move($repositoryFNA, $repositorySequence . $fldSep. $repositoryFNA) or die "error move file $!:";
		move($repositoryGenbank, $repositorySequence . $fldSep. $repositoryGenbank) or die "error move file $!:";
		move($repositoryReport, $repositorySequence . $fldSep. $repositoryReport) or die "error move file $!:";
		unlink glob "*.gz"  or die "for file *.gz $!:";
	} 
	
}
sub download_assembly_or_project {
	my ($sequenceId, $ftpServor, $fldSep, $directory) = @_;
	
	my $assemblySummary;
	my @sequenceIdList = split /,/, $sequenceId;
	
	if ($directory =~ /refseq/) {
		$assemblySummary = "assembly_summary_refseq.txt";
	} elsif ($directory =~ /genbank/) {
		$assemblySummary = "assembly_summary_genbank.txt";
	}
	
	my $assemblySummaryPath = "/genomes/ASSEMBLY_REPORTS/".$assemblySummary; 
	download_file($ftpServor, $assemblySummaryPath);
	
	foreach my $sequence (@sequenceIdList) {
		get_assembly_or_project($assemblySummary, $sequence, $ftpServor, $fldSep);
	}	
}

sub isModuleInstalled {
  my $mod = shift;

  #eval("use $mod");
  my $commandModule = `perldoc -l $mod`;
  
  if ($commandModule) {
    return(1);
  } else {
    return(0);
  }
}
