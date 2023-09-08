#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use Pod::Usage qw( pod2usage );
use Getopt::Long qw( :config gnu_getopt );
#use version; my $VERSION = qv('0.0.1');
my $VERSION = '0.0.1';
use English qw( -no_match_vars );
use List::Util qw( max min sum );
use POSIX qw( floor );

my %config = (tail => 10);
GetOptions(
   \%config, 'usage', 'help', 'man', 'version',
   qw(
     mind|min-data-points|m=i noise|N=s percentual|p!
     numeric|n! step|s|d=s include_zero|include-zero|zero|z!
     tail|t=i
     )
);
pod2usage(message => "$0 $VERSION", -verbose => 99, -sections => '')
  if $config{version};
pod2usage(-verbose => 99, -sections => 'USAGE') if $config{usage};
pod2usage(-verbose => 99, -sections => 'USAGE|EXAMPLES|OPTIONS')
  if $config{help};
pod2usage(-verbose => 2) if $config{man};

# Script implementation here
if ($config{noise}) {
   my ($noise_level, $perc) = $config{noise} =~ m{\A ([\d.]+) (%?) \z}mxs
     or pod2usage(-verbose => 99, -sections => 'USAGE');
   $config{noise}     = $noise_level;
   $config{noiseperc} = $perc;
} ## end if ($config{noise})

$config{step} = 1 if $config{numeric} && !$config{step};

my %hits;
while (<>) {
   chomp;
   $_ = floor($_ / $config{step})
     if $config{step} && $config{step} != 1;
   ++$hits{$_};
} ## end while (<>)

my @labels;
my @frequencies;
if ($config{numeric}) {
   my $min = min keys %hits;
   $min = min 0, $min if $config{include_zero};
   my $max = max keys %hits;
   $max = max 0, $max if $config{include_zero};
   while ($min <= $max) {
      push @labels, $min * $config{step};
      push @frequencies, exists($hits{$min}) ? $hits{$min} : 0;
      ++$min;
   }
} ## end if ($config{numeric})
else {

   # Order by frequency, descreasing
   @labels = reverse sort { $hits{$a} <=> $hits{$b} } keys %hits;
   @frequencies = @hits{@labels};
} ## end else [ if ($config{numeric})

$_ ||= 0 for @frequencies;    # eliminate undef's
my $sum;
$sum = sum(@frequencies) if $config{percentual};

if ($config{noise}) {
   my $noise_level =
     $config{noiseperc}
     ? (max(@frequencies) * $config{noise} / 100)
     : $config{noise};

   my $noise_index;
   for my $index (0 .. $#frequencies) {
      if ($frequencies[$index] >= $noise_level) {
         $noise_index = undef;
      }
      elsif (!defined $noise_index) {
         $noise_index = $index;
      }
   } ## end for my $index (0 .. $#frequencies)

   if (defined $noise_index) {
      $noise_index =
        max($noise_index + $config{tail}, $config{mind} || 0);
      splice @frequencies, $noise_index 
         if $noise_index < scalar @frequencies;
   } ## end if (defined $noise_index)
} ## end if ($config{noise})

if ($config{percentual}) {
   my $ratio = 100 / $sum;
   $_ *= $ratio for @frequencies;
}

print {*STDOUT} $labels[$_], "\t", $frequencies[$_], "\n"
  for 0 .. $#frequencies;

__END__

=head1 NAME

histogram - calculate frequency histogram out of lists of stuff

=head1 VERSION

Ask the version number to the script itself, calling:

   shell$ histogram --version

=head1 USAGE

   histogram [--usage] [--help] [--man] [--version]
   histogram [--include-zero|--zero|-z] [--min-data-points|-m <num>]
             [--noise|-N <threshold>] [--numeric|-n]
             [--percentual|-p] [--step|-s <step>] [--tail|-t <length>]

=head1 EXAMPLES

   shell$ histogram

   # Generating histogram's data
   shell$ grep 'interesting' file.txt | gawk '{print $3}' | histogram

   # If you have numbers you can keep ordering a divide into "classes"
   shell$ histogram --step 10 --numeric data-column.txt

=head1 DESCRIPTION

This utility produces histograms out of input data. Every line in input
(except the optional newline) is regarded as an item, and a count
is kept for each different item. After the counting phase is terminated,
the label-count pairs are printed in output, ordered by count, descending.
This is the basic work mode.

If you happen to know that your inputs are numbers, and you care about
keeping them in order, you can specify C<--numeric|-n>. This will make
sure that you'll have something resembling a distribution, and also that
all the gaps between will be filled (at integer intervals). If also
want 0 to be included, however far it may be, just pass option
C<--include-zero|--zero|-z>.

Moreover, if your data are numeric, and you'd rather group them by
steps (e.g. 0-9, 10-19, ecc.) you can pass option C<--step|-s>. Steps
start all from 0, and need not be integer ones.

=head1 OPTIONS


=over

=item --help

print a somewhat more verbose help, showing usage, this description of
the options and some examples from the synopsis.

=item B<< --include-zero | --zero | -z >>

in numeric mode (see C<--numeric|-n>) ensure that 0 is included in the
x axis of the output distribution.

=item --man

print out the full documentation for the script.

=item B<< --min-data-points | -m <num> >>

in numeric mode, set the minimum number of data points in output.

=item B<< --noise | -N <threshold> >>

set the noise level which can be used to cut the final output. If you
set a noise level, tails below the noise level are cut out and not put
in the output, yielding some kind of zoom on the area "where the meat
is". See also C<--tail|-t>.

=item B<< --numeric | -n >>

set numeric mode. Each input label is actually a number, so the output
is ordered by label (numerically) and not by frequency. Moreover,
see options C<--step|-s>, C<--include-zero|--zero|-z>.

=item B<< --percentual | -p >>

output percentual frequencies instead of absolute values. Percentuals
are calculated over the whole data set, i.e. including elements that
could be cut out by C<--noise|-n>, so they could not sum up to 100.

=item B<< --step | -s <step> >>

in numeric mode, divide the x axis by steps this long, and work on
"classes" instead of values. Useful to group stuff and have tighter
histograms.

=item B<< --tail | -t <length> >>

when a noise threshold is set, output at least C<length> elements
after the cut point, so that a bit of the tail is shown.

=item --usage

print a concise usage line and exit.

=item --version

print the version of the script.

=back

=head1 CONFIGURATION AND ENVIRONMENT

histogram requires no configuration files or environment variables.

=head1 DEPENDENCIES

None, apart L<version> that's core starting from 5.10. You can safely
modify the relevant line if you don't have it and you don't want to
install it.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests through http://rt.cpan.org/

Currently, the tail is defined only for rightmost elements, not for
leftmost. This can be a problem in numeric mode, where there could
be two tails (left and right).

=head1 AUTHOR

Flavio Poletti C<flavio@polettix.it>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2008, Flavio Poletti C<flavio@polettix.it>. All rights reserved.

This script is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>
and L<perlgpl>.

Questo script &#65533; software libero: potete ridistribuirlo e/o
modificarlo negli stessi termini di Perl stesso. Vedete anche
L<perlartistic> e L<perlgpl>.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=head1 NEGAZIONE DELLA GARANZIA

Poiché questo software viene dato con una licenza gratuita, non
c'è alcuna garanzia associata ad esso, ai fini e per quanto permesso
dalle leggi applicabili. A meno di quanto possa essere specificato
altrove, il proprietario e detentore del copyright fornisce questo
software "così com'è" senza garanzia di alcun tipo, sia essa espressa
o implicita, includendo fra l'altro (senza però limitarsi a questo)
eventuali garanzie implicite di commerciabilità e adeguatezza per
uno scopo particolare. L'intero rischio riguardo alla qualità ed
alle prestazioni di questo software rimane a voi. Se il software
dovesse dimostrarsi difettoso, vi assumete tutte le responsabilità
ed i costi per tutti i necessari servizi, riparazioni o correzioni.

In nessun caso, a meno che cià non sia richiesto dalle leggi vigenti
o sia regolato da un accordo scritto, alcuno dei detentori del diritto
di copyright, o qualunque altra parte che possa modificare, o redistribuire
questo software così come consentito dalla licenza di cui sopra, potrò
essere considerato responsabile nei vostri confronti per danni, ivi
inclusi danni generali, speciali, incidentali o conseguenziali, derivanti
dall'utilizzo o dall'incapacità di utilizzo di questo software. Ciò
include, a puro titolo di esempio e senza limitarsi ad essi, la perdita
di dati, l'alterazione involontaria o indesiderata di dati, le perdite
sostenute da voi o da terze parti o un fallimento del software ad
operare con un qualsivoglia altro software. Tale negazione di garanzia
rimane in essere anche se i dententori del copyright, o qualsiasi altra
parte, è stata avvisata della possibilità di tali danneggiamenti.

Se decidete di utilizzare questo software, lo fate a vostro rischio
e pericolo. Se pensate che i termini di questa negazione di garanzia
non si confacciano alle vostre esigenze, o al vostro modo di
considerare un software, o ancora al modo in cui avete sempre trattato
software di terze parti, non usatelo. Se lo usate, accettate espressamente
questa negazione di garanzia e la piena responsabilità per qualsiasi
tipo di danno, di qualsiasi natura, possa derivarne.
=cut
