#!/usr/bin/perl 
use warnings; 
use strict; 
use Tk;   
use utf8;
 
## Main program

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

# window creation 
my $window = new MainWindow( 
	-title      => 'getSequenceInfo', 
	-background => 'white', 
); 
  
# Taille minimale de ma fenêtre 
$window->minsize( 650, 600 ); 
  
my $homeMessage = "Thanks for using getSequenceInfo tool\n\n"; 
  
# Affichage d'un texte 
my $homeLabel = $window->Label( 
  -text       => $homeMessage, 
  -background => 'white', 
)->pack(); 

#  label and listbox for kingdom
my $kingdomLabel = $window->Label(
	-text => 'select the name of the kingdom : ',
	-background => 'white',
)->place(-x => 80, -y => 30);


my @kingdomList = qw/archaea bacteria fungi invertebrate plant protozoa 
									vertebrate_mammalian vertebrate_other viral/;

my $kingdomPos = 50;

foreach my $subKingdom(@kingdomList) {
	my $kingdomRadio = $window->Radiobutton(
		-text => $subKingdom,
		-value => $subKingdom,
		-variable => \$kingdom,
		-width => 20,
		-indicatoron => 0,
	)->place(-x => 80, -y => $kingdomPos);
	$kingdomPos += 22;
}

	
#label and entry for species
my $speciesLabel = $window->Label( 
	-text       => 'name of the species to search : ', 
	-background => 'white', 
)->place(-x => 80, -y => 270); 

my $speciesEntry = $window->Entry(
	-textvariable => \$species,
)->place(-x => 80, -y => 290); 

#label and entry for date
my $dateLabel = $window->Label( 
	-text       => 'assembly from this date (yyyy-mm-dd) : ', 
	-background => 'white', 
)->place(-x => 80, -y => 315); 

my $dateEntry = $window->Entry(
	-textvariable => \$date,
)->place(-x => 80, -y => 335); 

#label and listbox for summary
my $summaryLabel = $window->Label( 
	-text       => 'obtain last updates from the servor : ', 
	-background => 'white',
)->place(-x => 80, -y => 360); 

my $summaryCheckButton = $window->Checkbutton(
	-text => 'yes',
	-onvalue => 'yes',
	-offvalue => 'no',
	-variable => \$getSummary,
	-background => 'white'
)->place(-x => 80, -y => 380);

#label and entry for taxid
my $taxidLabel = $window->Label( 
	-text       => 'taxid of the species to search : ', 
	-background => 'white',
)->place(-x => 80, -y => 400);

my $taxidEntry = $window->Entry(
	-textvariable => \$taxID,
)->place(-x => 80, -y => 420);  

# label and entry for representation
my $representationLabel = $window->Label( 
	-text       => 'assembly level   : ', 
	-background => 'white',
)->place(-x => 420, -y => 30); 


my @representationList = ("Complete Genome", "Chromsome", "Scaffold", "Contig"); 

my $representationPos = 50;

foreach my $subRepresentation (@representationList) {
	my $representationRadio = $window->Radiobutton(
		-text => $subRepresentation,
		-value => $subRepresentation,
		-variable => \$representation,
		-width => 20,
		-indicatoron => 0,
	)->place(-x => 420, -y => $representationPos);
	$representationPos += 22;
}

# label and entry for componenents
my $componenentsLabel = $window->Label( 
	-text       => 'component : ',
	-background => 'white', 
)->place(-x => 420, -y => 145); 

my @componentList = qw/chromosome plasmid contig scaffold/; 

my $componentPos = 165;

foreach my $subComponent (@componentList) {
	my $componentRadio = $window->Radiobutton(
		-text => $subComponent,
		-value => $subComponent,
		-variable => \$component,
		-width => 20,
		-indicatoron => 0,
	)->place(-x => 420, -y => $componentPos);
	$componentPos += 22;
}

# label and entry for quantity
my $quantityLabel = $window->Label( 
	-text       => 'quantity of assembly : ', 
	-background => 'white', 
)->place(-x => 420, -y => 255); 

my $quantityEntry = $window->Entry( 
	-textvariable => \$quantity,
)->place(-x => 420, -y => 275);  

# label and entry for ena
my $enaLabel = $window->Label( 
	-text       => 'download sequences (ena ID)  : ', 
	-background => 'white', 
)->place(-x => 420, -y => 295); 

my $enaEntry = $window->Entry( 
	-textvariable => \$enaID,
)->place(-x => 420, -y => 315); 

# label and entry to rename the output
my $outputLabel = $window->Label( 
	-text       => 'name of the folder : ', 
	-background => 'white', 
)->place(-x => 420, -y => 340); 

my $outputEntry = $window->Entry(
	-textvariable => \$output,
)->place(-x => 420, -y => 360);  

# Affichage d'un bouton pour fermer la fenêtre 
my $button = $window->Button( 
	-text    => 'start search', 
	-command => \&search, 
)->place(-x => 300, -y => 440); 


MainLoop;


sub search {
	my $command = "perl getSequenceInfo.pl";
	
	if (defined $getSummary) {$getSummary = $getSummary eq "yes" ? "-get" : ""; }
	
	if (defined $enaID) {
		system("$command -ena $enaID");
	}
	else {
		if (defined $kingdom) { $command .= " -k $kingdom"; }
		if (defined $species) { $command .= " -s \"$species\""; }
		if (defined $taxID) { $command .= " -taxid $taxID"; }
		if (defined $quantity) { $command .= " -q $quantity"; }
		if (defined $representation) { $command .= " -r \"$representation\""; }
		if (defined $component) { $command .= " -c $component";}
		if (defined $date) { $command .= " -date $date"; }
		print "$command $getSummary\n";
		system("$command $getSummary");
	}
}



