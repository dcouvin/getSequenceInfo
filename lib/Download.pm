package Download;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT = qw(download_file obtain_file clean_repository);

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
	if ($link =~ /$servor(.*)/) {return ($1);}
}
#------------------------------------------------------------------------------
#   delete all compress folder download
sub clean_repository {
	unlink glob "*.gz *.dmp"  or die "$!: for file *.gz *.dmp";
}

1;