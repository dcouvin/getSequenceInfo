ECHO OFF
ECHO Please note that Perl (http://strawberryperl.com/) must be installed on your computer before running the following commands.
ECHO Installation of needed Perl Modules:
set PATH=C:\strawberry\c\bin;%PATH%
cpan -f -i Tk
cpan -f -i BioPerl
cpan -f -i Date::Calc
cpan -f -i Bio::SeqIO
cpan -f -i LWP::Simple
cpan -f -i Data::Dumper
cpan -f -i IO::Uncompress::Gunzip
cpan -f -i IO::File
cpan -f -i Getopt::Long
cpan -f -i Net::FTP
ECHO Finished... Please press "Enter" to close this window.
PAUSE
