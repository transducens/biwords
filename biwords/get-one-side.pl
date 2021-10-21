use strict;
use warnings;

# Getting command line arguments:
use Getopt::Long;
# Documentation:
use Pod::Usage;
# I/O Handler
use IO::Handle;

#use locale;
#use POSIX qw(locale_h);

#print STDERR "LOCALE: ", setlocale(LC_ALL,""), "\n";

#Vars for the command line options
my($alignfile, $side, $gzip);

# Command line arguments
GetOptions( 'align|a=s' => \$alignfile, 'side|s=n' => \$side, "gzip|z" => \$gzip);

exit print STDERR "Sintax: $0 -align|-a alignments -side|-s 1|2|3 [-gzip|-z]\n" unless ($alignfile && $side);

print STDERR "Alignment file: '$alignfile'\n";

my $cmd="cat";
$cmd = "zcat" if ($gzip);

exit print STDERR "Error: Cannot open input file '$alignfile': $!\n" unless (open(ALIGN, "$cmd $alignfile|"));

while (<ALIGN>) {
  chomp;

  my @fields = split /\|/;
  $_ = $fields[$side];
  
  s/^[ ]+//g;
  s/[ ]+$//g;

  print $_, "\n";
}

close(ALIGN);

