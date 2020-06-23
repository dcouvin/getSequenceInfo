#!/usr/bin/perl 


use warnings; 
use strict; 
use Tk;
use Tk::ProgressBar;
use Tk::Labelframe;
use utf8;

# variable declaration
my $kingdom;
my $species;
my $date;
my $getSummaries;
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
 
 $mw->minsize(550, 550); 

# a Frame you'll never see
my $frame1 = $mw->Labelframe(-background => 'white');
my $frame2 = $mw->Labelframe(-background => 'white',  -text => 'NCBI', -relief => 'groove', -borderwidth => '2');
my $firstSubFrame2 = $frame2->Frame(-background => 'white');
my $subSubFrame1 = $firstSubFrame2->Frame(-background => 'white');
my $subSubFrame2 = $firstSubFrame2->Frame(-background => 'white');
my $directoryFrame = $subSubFrame1->Frame(-background => 'white', -borderwidth => '5');
my $kingdomFrame = $subSubFrame2->Frame(-background => 'white', -borderwidth => '5');
my $levelFrame = $subSubFrame2->Frame(-background => 'white',  -borderwidth => '5');
my $componentFrame = $subSubFrame1->Frame(-background => 'white', -borderwidth => '5');
my $secondSubFrame2 = $frame2->Frame(-background => 'white');
my $frame3 = $mw->Labelframe(-background => 'white', -text => 'ENA');
my $frame4 = $mw->Labelframe(-background => 'white');

my $homeMessage = "Thanks for using getSequenceInfo tool\n\n"; 

#Affichage d'un texte 
$frame1 = $mw->Label( 
  -text       => $homeMessage, 
  -background => 'white', 
)->pack();


##subframe1
$directoryFrame->Label(
	-text => 'directory : ',
	-background => 'white', 
)->pack(-expand => 1, -fill => 'x');

my @directoryList = qw/refseq genbank/;

foreach my $directory (@directoryList) {
	$directoryFrame->Radiobutton(
		-text => $directory,
		-value => $directory,
		-variable => \$directoryNcbi,
	)->pack(-expand => 1, -fill => 'x');
}

#label and entry for representation
$levelFrame->Label( 
	-text       => 'assembly level   : ', 
	-background => 'white',
)->pack(-expand => 1, -fill => 'x'); 

my @levelList = ("Complete Genome", "All"); 

foreach my $level (@levelList) {
	$levelFrame->Radiobutton(
		-text => $level,
		-value => $level,
		-variable => \$levelAssembly,
	)->pack(-expand => 1, -fill => 'x');
}

#  label and listbox for kingdom
$kingdomFrame->Label(
	-text => 'Kingdom : ',
	-background => 'white',
)->pack(-expand => 1, -fill => 'x');


my $lst = $kingdomFrame->Scrolled("Listbox", -scrollbars => 'oe')->pack(-expand => 1, -fill => 'x'); 
$lst->configure(-height => 5);

my @kingdomList = qw/archaea bacteria fungi invertebrate plant protozoa 
	vertebrate_mammalian vertebrate_other viral/;
	
$lst->insert('end', @kingdomList);

# label and entry for componenents
$componentFrame->Label( 
	-text       => 'component : ',
	-background => 'white', 
)->pack(-expand => 1, -fill => 'x'); 

my %checkComponentHash = (
	chromosome => $checkChromosome,
	plasmid => $checkPlasmid,
	contig => $checkContig,
	scaffold => $checkScaffold,
);

foreach my $component (keys %checkComponentHash) {
	$componentFrame->Checkbutton(
	-text => $component,
	-onvalue => 1,
	-offvalue => 0,
	-relief => 'raised',
	-width => 17,
	-variable => \$checkComponentHash{$component},
	)->pack(-expand => 1, -fill => 'x'); 
}

$componentFrame->Entry( 
	-textvariable => \$specificComponent,
)->pack(-expand => 1, -fill => 'x');

## subframe2
# label and entry for species
$secondSubFrame2->Label( 
	-text       => 'Species : ', 
	-background => 'white', 
)->grid($secondSubFrame2->Entry(-textvariable => \$species),
-sticky => 'nsew');

#label and entry for date
$secondSubFrame2->Label( 
	-text       => 'Assembly date (yyyy-mm-dd) : ', 
	-background => 'white', 
)->grid($secondSubFrame2->Entry(-textvariable => \$date),
-sticky=> 'nsew'); 

#label and entry for taxid
$secondSubFrame2->Label( 
	-text       => 'Taxid : ', 
	-background => 'white',
)->grid($secondSubFrame2->Entry(-textvariable => \$taxID),
-sticky =>'nsew');  

# label and entry to rename the output
$secondSubFrame2->Label( 
	-text       => 'Name of the folder : ', 
	-background => 'white', 
)->grid($secondSubFrame2->Entry(-textvariable => \$output),
-sticky => 'nsew'); 

# label and entry for quantity
$secondSubFrame2->Label( 
	-text       => 'quantity of assembly : ', 
	-background => 'white', 
)->grid($secondSubFrame2->Entry(-textvariable => \$quantity),
-sticky => 'nsew');  

# label and entry for quantity
$secondSubFrame2->Label( 
	-text       => 'assembly or project ID : ', 
	-background => 'white', 
)->grid($secondSubFrame2->Entry(-textvariable => \$assemblyPrjID),
-sticky => 'nsew');

# label and entry for database
$secondSubFrame2->Label( 
	-text       => 'Database to save : ', 
	-background => 'white', 
)->grid($secondSubFrame2->Entry(-textvariable => \$getSummaries),
-sticky => 'nsew');


## ENA frame
# label and entry for ena
$frame3->Label( 
	-text       => 'download sequences (ena ID)  : ', 
	-background => 'white', 
)->grid($frame3->Entry(-textvariable => \$enaID),
-sticky => 'nsew'); 


$frame3->Label( 
	-text       => 'run accession id for fastq : ', 
	-background => 'white', 
 )->grid( $frame3->Entry(-textvariable => \$fastqID),
 -sticky => 'nsew'); 
 
## bottom frame
# start search button
$frame4->Button( 
	-text    => 'Search',
	-command => \&search, 
)->grid(my $progress = $frame4->ProgressBar(
	-from   => 0, 
	-to => 100,   
	-length => 300, 
	-colors => [ 0, 'blue',]),
    -sticky => "nsew"
);

my @frameList = (
	$frame3,
	$frame4,
	$secondSubFrame2
);

foreach my $frame (@frameList) {
	my ($columns, $rows) = $frame->gridSize();
	
	for (my $i = 0; $i < $columns; $i++) {
		$frame->gridColumnconfigure($i, -weight => 1); 
	}
	for (my $i = 0; $i < $rows; $i++) { 
		$frame->gridRowconfigure($i, -weight => 1); 
	}
}
 

$frame1->pack(-expand => 1, -fill => 'x');
$frame2->pack(-expand => 1, -fill => 'x');
$firstSubFrame2->pack(-expand => 1, -fill => 'x');
$subSubFrame1->pack(-side => 'left', -expand => 1, -fill => 'both');
$subSubFrame2->pack(-side => 'right', -expand => 1, -fill => 'both');
$directoryFrame->pack(-expand => 1, -fill => 'x');
$kingdomFrame->pack(-expand => 1, -fill => 'x');
$levelFrame->pack(-expand => 1, -fill => 'x');
$componentFrame->pack(-expand => 1, -fill => 'x');
$secondSubFrame2->pack(-expand => 1, -fill => 'x');
$frame3->pack(-expand => 1, -fill => 'x');
$frame4->pack(-expand => 1, -fill => 'x');


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
		my @listIndex = $lst->curselection;
		
		foreach (@listIndex) {
			$kingdom = $kingdomList[$_];
		}
		
		my @optionCharList = ( '-k', '-s', '-taxid', '-q', '-l', '-c', '-o', '-dir', '-date', '-get');
		
		my %optionHash = (
			'-k' => $kingdom,
			'-s' => $species,
			'-taxid' => $taxID,
			'-q' => $quantity,
			'-l' => $levelAssembly,
			'-c' => $components,
			'-o' => $output,
			'-dir' => $directoryNcbi,
			'-date' => $date,
			'-get' => $getSummaries
		);
		
		foreach my $option (@optionCharList) {
			if (defined  $optionHash{$option}) { $command .= " $option \"$optionHash{$option}\""; }
			$i += 10;
			$progress->value($i); 
			$mw->update();
			sleep 1;
		}
	
		if ($outputFile) { $command.= " -log"; }
		
		print "$command\n";
		system("$command");
		
		undef $components;
		undef $outputFile;
		
		$i += 10;
		$progress->value($i); 
		$mw->update();
		sleep 1;
	}
}
