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
); 

#  label and listbox for kingdom
my $kingdomLabel = $window->Label(
	-text => 'select the name of the kingdom : ',
	-background => 'white',
);


my @kingdomList = qw/archaea bacteria fungi invertebrate plant protozoa
										vertebrate_mammalian vertebrate_other 
										viral/;


my $kingdomListbox = $window->Scrolled(
	"Listbox",
	-scrollbars => 'oe',
	-selectmode => 'single'
); 

$kingdomListbox->insert("end", @kingdomList); 

$kingdomListbox->configure(-height => 5);

my $kingdomButton = $window->Button(
	-text => 'add',
	-command => \&addKingdom,
);

#label and entry for species
my $speciesLabel = $window->Label( 
	-text       => 'name of the species to search : ', 
	-background => 'white', 
); 

my $speciesEntry = $window->Entry(
	-textvariable => \$species,
); 

#label and entry for date
my $dateLabel = $window->Label( 
	-text       => 'assembly from this date (yyyy-mm-dd) : ', 
	-background => 'white', 
); 

my $dateEntry = $window->Entry(
	-textvariable => \$date,
); 

#label and listbox for summary
my $summaryLabel = $window->Label( 
	-text       => 'obtain last updates from the servor : ', 
	-background => 'white',
); 

my $summaryCheckButton = $window->Checkbutton(
	-text => 'yes',
	-onvalue => 'yes',
	-offvalue => 'no',
	-variable => \$getSummary,
	-background => 'white'
);

#label and entry for taxid
my $taxidLabel = $window->Label( 
	-text       => 'taxid of the species to search : ', 
	-background => 'white',
);

my $taxidEntry = $window->Entry(
	-textvariable => \$taxID,
);  

# label and entry for representation
my $representationLabel = $window->Label( 
	-text       => 'assembly level   : ', 
	-background => 'white',
); 


my $representationListbox = $window->Scrolled( 
	"Listbox",
	-scrollbars => 'oe',
	-selectmode => 'single'
);

my @representationList = ("Complete Genome", "Chromsome", "Scaffold", "Contig"); 

$representationListbox->insert('end', @representationList);

$representationListbox->configure(-height => 3);

my $representationButton = $window->Button(
	-text    => 'add', 
	-command => \&addRepresentation, 
); 

# label and entry for componenents
my $componenentsLabel = $window->Label( 
	-text       => 'component : ',
	-background => 'white', 
); 

my @componentList = qw/chromosome plasmid contig scaffold/; 

my $componentsListbox = $window->Scrolled( 
	"Listbox",
	-scrollbars => 'oe',
	-selectmode => 'single'
);

$componentsListbox->insert('end', @componentList);

$componentsListbox->configure(-height => 3);

my $componentButton = $window->Button(
	-text => 'add',
	-command => \&getComponent,
);

# label and entry for quantity
my $quantityLabel = $window->Label( 
	-text       => 'quantity of assembly : ', 
	-background => 'white', 
); 

my $quantityEntry = $window->Entry( 
	-textvariable => \$quantity,
);  

# label and entry for ena
my $enaLabel = $window->Label( 
	-text       => 'download sequences (ena ID)  : ', 
	-background => 'white', 
); 

my $enaEntry = $window->Entry( 
	-textvariable => \$enaID,
); 

# label and entry to rename the output
my $outputLabel = $window->Label( 
	-text       => 'name of the folder : ', 
	-background => 'white', 
); 

my $outputEntry = $window->Entry(
	-textvariable => \$output,
);  

# Affichage d'un bouton pour fermer la fenêtre 
my $button = $window->Button( 
	-text    => 'search', 
	-command => \&search, 
); 

# placement
$homeLabel->pack();
$kingdomLabel->place(-x => 80, -y => 30);
$kingdomListbox->place( -x => 80, -y => 50);
$kingdomButton->place(-x => 80, -y => 145);
$speciesLabel->place(-x => 80, -y => 175);
$speciesEntry->place(-x => 80, -y => 195);
$dateLabel->place(-x => 80, -y => 215);
$dateEntry->place(-x => 80, -y => 235);
$summaryLabel->place(-x => 80, -y => 255);
$summaryCheckButton->place(-x => 80, -y => 275);
$taxidLabel->place(-x => 80, -y => 297);
$taxidEntry->place(-x => 80, -y => 315);
$representationLabel->place(-x => 420, -y => 30);
$representationListbox->place(-x => 420, -y => 50);
$representationButton->place(-x => 420, -y => 110);
$componenentsLabel->place(-x => 420, -y => 145);
$componentsListbox->place(-x => 420, -y => 165);
$componentButton->place(-x => 420, -y => 225);
$quantityLabel->place(-x => 420, -y => 255);
$quantityEntry->place(-x => 420, -y => 275);
$enaLabel->place(-x => 420, -y => 295);
$enaEntry->place(-x => 420, -y => 315);
$outputLabel->place(-x => 420, -y => 340);
$outputEntry->place(-x => 420, -y => 360);
$button->place(-x => 300, -y => 400);

MainLoop;

sub addKingdom {
	my @selection = $kingdomListbox-> curselection;
	my $index = pop @selection;
	$kingdom = $kingdomList[$index];
}

sub addRepresentation {
	my @selection = $representationListbox-> curselection;
	my $index = pop @selection;
	$representation = $representationList[$index];
}

sub getComponent {
	my @selection = $componentsListbox-> curselection;
	my $index = pop @selection;
	$component = $componentList[$index];
}

sub search {
	my $command = "perl getSequenceInfo.pl";
	$getSummary = $getSummary eq "yes" ? "-get" : "";
	
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



