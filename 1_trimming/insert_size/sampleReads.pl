use strict;
use warnings;


my $num = 1000000;
open READS, "gzip -cd $ARGV[0] |" or die "$!.\n";
while (my $id = <READS>, my $sq = <READS>, <READS>, my $qu = <READS>)
{
    last if $num == 0;

    printf "%s%s+\n%s", $id, $sq, $qu;

    $num --;
}
close READS;
