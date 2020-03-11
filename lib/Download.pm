package Download;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT = qw(download_file obtain_file);

# connect and download file from ftp servor
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


1;