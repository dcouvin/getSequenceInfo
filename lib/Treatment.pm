package Treatment;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT = qw(trim trim_array empty_folder);


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
sub empty_folder {!<$_[0]/*>}

1;
