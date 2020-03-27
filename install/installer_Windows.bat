@ECHO OFF
ECHO Installation of needed Perl Modules:
ECHO Please note that ActivePerl (https://www.activestate.com/products/activeperl/downloads/) must be installed on your Microsoft Windows computer before running the following commands.
ppm install Date::Calc
ppm install Bio::SeqIO
ppm install LWP::Simple
ppm install Data::Dumper
ppm install IO::Uncompress::Gunzip
ppm install IO::File
ppm install Getopt::Long
ppm install Net::FTP
ppm install Tk
ECHO Installation has been executed successfully. Please press "Enter" to close this window.
PAUSE