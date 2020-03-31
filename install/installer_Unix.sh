#!/usr/bin/bash

echo 'Welcome to getSequenceInfo Unix installer !';
echo 'cpan must first be installed on your computer';
echo '----------------------------------------------------------';
sudo cpan -i Date::Calc;
sudo cpan -i Bio::SeqIO
sudo cpan -i LWP::Simple;
sudo cpan -i Data::Dumper;
sudo cpan -i IO::Uncompress::Gunzip;
sudo cpan -i IO::File;
sudo cpan -i Getopt::Long;
sudo cpan -i Net::FTP;
sudo cpan -i Tk;
echo 'end of install';

