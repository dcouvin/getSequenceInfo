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
my $representation;
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
	-text => 'directory where to download sequences : ',
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
	-text => 'select the name of the kingdom : ',
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
	-text       => 'name of the species to search : ', 
	-background => 'white', 
)->pack(); 

$frame2 ->Entry(
	-textvariable => \$species,
)->pack();

#label and entry for date
$frame2->Label( 
	-text       => 'assembly from this date (yyyy-mm-dd) : ', 
	-background => 'white', 
)->pack(); 

$frame2->Entry(
	-textvariable => \$date,
)->pack(); 

#label and listbox for summary
$frame2->Label( 
	-text       => 'obtain last updates from the servor : ', 
	-background => 'white',
)->pack(); 

$frame2->Checkbutton(
	-text => 'yes',
	-onvalue => 'yes',
	-offvalue => 'no',
	-relief => 'raised',
	-variable => \$getSummary,
)->pack();

#label and entry for taxid
$frame2->Label( 
	-text       => 'taxid of the species to search : ', 
	-background => 'white',
)->pack();

$frame2->Entry(
	-textvariable => \$taxID,
)->pack();  


# label and entry to rename the output
$frame2->Label( 
	-text       => 'name of the folder : ', 
	-background => 'white', 
)->pack(); 

$frame2->Entry(
	-textvariable => \$output,
)->pack(); 


# #right frame side
#label and entry for representation
$frame3->Label( 
	-text       => 'assembly level   : ', 
	-background => 'white',
)->pack(); 


my @representationList = ("Complete Genome", "All"); 


foreach my $subRepresentation (@representationList) {
	$frame3->Radiobutton(
		-text => $subRepresentation,
		-value => $subRepresentation,
		-variable => \$representation,
		-width => 20,
		-indicatoron => 0,
	)->pack();
}

# label and entry for componenents
$frame3->Label( 
	-text       => 'component : ',
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
	-text       => 'quantity of assembly : ', 
	-background => 'white', 
)->pack(); 

$frame3->Entry( 
	-textvariable => \$quantity,
)->pack();  

# label and entry for ena
$frame5->Label( 
	-text       => 'ENA', 
	-background => 'white', 
)->pack(); 

# label and entry for ena
$frame5->Label( 
	-text       => 'download sequences (ena ID)  : ', 
	-background => 'white', 
)->pack(); 

$frame5->Entry( 
	-textvariable => \$enaID,
)->pack(); 


# label and entry to rename the output
$frame5->Label( 
	-text       => 'run accession id for fastq : ', 
	-background => 'white', 
)->pack(); 

$frame5->Entry(
	-textvariable => \$fastqID,
)->pack(); 




## bottom frame
# start search button
$frame4->Button( 
	-text    => 'start search', 
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
	
	if (defined $getSummary) { $getSummary = $getSummary eq "yes" ? "-get" : ""; }
	
	foreach my $component (keys %checkComponentHash) {
		if ($checkComponentHash{$component}) {
			$components.= $component . ",";
		}
	}
	
	if ($representation && $representation =~ /all/i) {
		$representation = "Complete Genome,Chromosome,Scaffold,Contig";
	}
	
	if ($specificComponent && $specificComponent !~ /keyword/) {$components.= $specificComponent;}
	
	if (defined $enaID) {
		print "$command -ena $enaID\n";
		system("$command -ena $enaID");
		$i += 100;
		$progress->value($i); 
		$mw->update();
		sleep 1; 
	}
	elsif (defined $fastqID) {
		print "$command -fastq $fastqID\n";
		system("$command -fastq $fastqID");
		$i += 100;
		$progress->value($i); 
		$mw->update();
		sleep 1;
	}
	else {
		my @optionCharList = ( '-k', '-s', '-taxid', '-q', '-r', '-c', '-o', '-dir', '-date');
		
		my %optionHash = (
			'-k' => $kingdom,
			'-s' => $species,
			'-taxid' => $taxID,
			'-q' => $quantity,
			'-r' => $representation,
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
	
		if (defined  $getSummary) { $command .= " $getSummary"; }
		
		print "$command\n";
		system("$command");
		
		undef $getSummary;
		undef $components;
		undef $enaID;
		
		$i += 10;
		$progress->value($i); 
		$mw->update();
		sleep 1;
	}
}
