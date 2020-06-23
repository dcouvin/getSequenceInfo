#!/usr/bin/perl 


use warnings; 
use strict; 
use Tk;
use Tk::ProgressBar; 
use utf8;
 

# variable declaration
my $kingdom;
my $species;
my $date;
my $getSummary;
my $levelAssembly;
my $component;
my $quantity;
my $enaID;
my $output;
my $taxID;
my $fastqID;
my $directoryNcbi;
my $checkChromosome;
my $checkPlasmid;
my $checkContig;
my $checkScaffold;
my $specificComponent = "keyword";
my $components;
my $assemblyPrjID;
my $outputFile;

# window creation 
my $mw = new MainWindow( 
	-title      => 'getSequenceInfo', 
	-background => 'white', 
); 

# Taille de ma fenêtre 
# $mw->minsize(600,620); 
  

# a Frame you'll never see
my $frame1 = $mw->Frame(-background => 'white');
my $frame2 = $mw->Frame(-background => 'white');
my $frame3 = $mw->Frame(-background => 'white'); 
my $frame4 = $mw->Frame(-background => 'white');
my $frame5 = $frame3->Frame(-background => 'white', -borderwidth => 5, -relief => 'groove'); 

my $homeMessage = "Thanks for using getSequenceInfo tool\n\n"; 

#Affichage d'un texte 
$frame1 = $mw->Label( 
  -text       => $homeMessage, 
  -background => 'white', 
)->pack();




## left frame side
#database where tto download sequences 
$frame2->Label(
	-text => 'NCBI sequence repository: ',
	-background => 'white', 
)->pack();


my @directoryList = qw/refseq genbank/;

foreach my $directory (@directoryList) {
	$frame2->Radiobutton(
		-text => $directory,
		-value => $directory,
		-variable => \$directoryNcbi,
		-indicatoron => 0,
		-width => 20,
	)->pack();
}


#  label and listbox for kingdom
$frame2->Label(
	-text => 'Select the Kingdom: ',
	-background => 'white',
)->pack();


my @kingdomList = qw/archaea bacteria fungi invertebrate plant protozoa 
									vertebrate_mammalian vertebrate_other viral/;


foreach my $subKingdom(@kingdomList) {
	 $frame2->Radiobutton(
		-text => $subKingdom,
		-value => $subKingdom,
		-variable => \$kingdom,
		-indicatoron => 0,
		-width => 20,
	)->pack();
}

#label and entry for species
$frame2->Label( 
	-text       => 'Species: ', 
	-background => 'white', 
)->pack(); 

$frame2 ->Entry(
	-textvariable => \$species,
)->pack();

#label and entry for date
$frame2->Label( 
	-text       => 'Assemblies from this date (yyyy-mm-dd): ', 
	-background => 'white', 
)->pack(); 

$frame2->Entry(
	-textvariable => \$date,
)->pack(); 

#label and listbox for summary
$frame2->Label( 
	-text       => 'Get the latest version of the summary: ', 
	-background => 'white',
)->pack(); 

$frame2->Checkbutton(
	-text => 'yes',
	-onvalue => 1,
	-offvalue => 0,
	-relief => 'raised',
	-width => 17,
	-variable => \$getSummary,
)->pack();

#label and entry for taxid
$frame2->Label( 
	-text       => 'NCBI Taxonomy ID (taxID): ', 
	-background => 'white',
)->pack();

$frame2->Entry(
	-textvariable => \$taxID,
)->pack();  


# label and entry to rename the output
$frame2->Label( 
	-text       => 'Output folder name: ', 
	-background => 'white', 
)->pack(); 

$frame2->Entry(
	-textvariable => \$output,
)->pack(); 


# #right frame side
#label and entry for representation
$frame3->Label( 
	-text       => 'Assembly level   : ', 
	-background => 'white',
)->pack(); 


my @levelList = ("Complete Genome", "All"); 


foreach my $level (@levelList) {
	$frame3->Radiobutton(
		-text => $level,
		-value => $level,
		-variable => \$levelAssembly,
		-width => 20,
		-indicatoron => 0,
	)->pack();
}

# label and entry for componenents
$frame3->Label( 
	-text       => 'Component: ',
	-background => 'white', 
)->pack(); 


my %checkComponentHash = (
	chromosome => $checkChromosome,
	plasmid => $checkPlasmid,
	contig => $checkContig,
	scaffold => $checkScaffold,
);


foreach my $component (keys %checkComponentHash) {
	$frame3->Checkbutton(
	-text => $component,
	-onvalue => 1,
	-offvalue => 0,
	-relief => 'raised',
	-width => 17,
	-variable => \$checkComponentHash{$component},
	)->pack(); 
}

$frame3->Entry( 
	-textvariable => \$specificComponent,
)->pack();  

# label and entry for quantity
$frame3->Label( 
	-text       => 'Number of assemblies (limit): ', 
	-background => 'white', 
)->pack(); 

$frame3->Entry( 
	-textvariable => \$quantity,
)->pack();  

# label and entry for quantity
$frame3->Label( 
	-text       => 'Assembly or Project ID: ', 
	-background => 'white', 
)->pack(); 

$frame3->Entry( 
	-textvariable => \$assemblyPrjID,
)->pack();

$frame3->Label( 
	-text       => 'Get a log file: ', 
	-background => 'white',
)->pack(); 

$frame3->Checkbutton(
	-text => 'yes',
	-onvalue => 1,
	-offvalue => 0,
	-relief => 'raised',
	-width => 17,
	-variable => \$outputFile,
)->pack();

# label and entry for ena
$frame5->Label( 
	-text       => 'ENA', 
	-background => 'white', 
)->pack(); 

# label and entry for ena
$frame5->Label( 
	-text       => 'ENA sequence ID (enaID): ', 
	-background => 'white', 
)->pack(); 

$frame5->Entry( 
	-textvariable => \$enaID,
)->pack(); 


# label and entry to rename the output
$frame5->Label( 
	-text       => 'FASTQ run accession: ', 
	-background => 'white', 
)->pack(); 

$frame5->Entry(
	-textvariable => \$fastqID,
)->pack(); 




## bottom frame
# start search button
$frame4->Button( 
	-text    => 'Search', 
	-command => \&search, 
)->pack();

#progress bar
my $progress = $frame4->ProgressBar( 
  -from   => 0, 
  -to     => 100, 
  -width  => 25,
  -length => 150,
  -colors => [ 0, 'blue',], 
)->pack(); 

$frame1->pack(-side => 'top', -fill => 'both');
$frame2->pack(-side => 'left', -fill => 'both', -expand => 1);
$frame3->pack(-side => 'right', -fill => 'both', -expand => 1);
$frame5->pack(-expand => 1);
$frame4->pack(-side => 'bottom', -fill => 'both');


MainLoop;


sub search {
	my $i = 0;
	
	my $command = "perl getSequenceInfo.pl";
	
	
	foreach my $component (keys %checkComponentHash) {
		if ($checkComponentHash{$component}) {
			$components.= $component . ",";
		}
	}
	
	if ($levelAssembly && $levelAssembly =~ /all/i) {
		$levelAssembly = "Complete Genome,Chromosome,Scaffold,Contig";
	}
	
	if ($specificComponent && $specificComponent !~ /keyword/) {$components.= $specificComponent;}
	
	if (defined $enaID) {
		$command .= " -ena $enaID";
		if ($outputFile) { $command .= " -log"; }
		print "$command\n";	
		system("$command");
		$i += 100;
		$progress->value($i); 
		$mw->update();
		sleep 1; 
		
		undef $enaID;
	}
	elsif (defined $fastqID) {
		$command .= " -fastq $fastqID";
		if ($outputFile) { $command .= " -log"; }
		print "$command\n";	
		system("$command");
		$i += 100;
		$progress->value($i); 
		$mw->update();
		sleep 1;
		
		undef $fastqID;
	}
	elsif (defined $assemblyPrjID) {
		$command .= " -assembly_or_project $assemblyPrjID";
		if ($outputFile) { $command .= " -log"; }
		print "$command\n";	
		system("$command");
		$i += 100;
		$progress->value($i); 
		$mw->update();
		sleep 1;
		
		undef  $assemblyPrjID;
	}
	else {
		my @optionCharList = ( '-k', '-s', '-taxid', '-q', '-l', '-c', '-o', '-dir', '-date');
		
		my %optionHash = (
			'-k' => $kingdom,
			'-s' => $species,
			'-taxid' => $taxID,
			'-q' => $quantity,
			'-l' => $levelAssembly,
			'-c' => $components,
			'-o' => $output,
			'-dir' => $directoryNcbi,
			'-date' => $date
		);
		
		foreach my $option (@optionCharList) {
			if (defined  $optionHash{$option}) { $command .= " $option \"$optionHash{$option}\""; }
			$i += 10;
			$progress->value($i); 
			$mw->update();
			sleep 1;
		}
	
		if ($outputFile) { $command.= " -log"; }
		if ($getSummary) { $command .= " -get"; }
		
		print "$command\n";
		system("$command");
		
		undef $getSummary;
		undef $components;
		undef $outputFile;
		
		$i += 10;
		$progress->value($i); 
		$mw->update();
		sleep 1;
	}
}
