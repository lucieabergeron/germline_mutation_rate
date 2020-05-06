use strict;
use warnings;
use File::Basename 'basename';
use Statistics::Descriptive;


my $samtools = '/com/extra/samtools/0.1.19/bin/samtools';

my @suffix1 = ('_1.clean.smp.bam', '_1.smp.bam');
my @suffix2 = ('_2.clean.smp.bam', '_2.smp.bam');


die "USAGE: $0 <read1_SE.bam> <read2_SE.bam>\n" if not @ARGV == 2;
die "The file with read2 is not the mate of the file with read1\nUSAGE: $0 <read1_SE.bam> <read2_SE.bam>\n"
 if (basename $ARGV[0], @suffix1) ne (basename $ARGV[1], @suffix2);


my $lane = basename $ARGV[0], @suffix1;

my %box;
open IN, '-|', "$samtools view -X -F 3844 $ARGV[0]" or die "$!\n";
while (<IN>)
{
    chomp;

    my ($read, $flag, $chr, $pos, $seq) = (split /\t/)[0 .. 3, 9]; my $len = length $seq;
    $box{$read}{'r1'} = [$chr, $pos, $flag, $len];
}
close IN;

open IN, '-|', "$samtools view -X -F 3844 $ARGV[1]" or die "$!\n";
while (<IN>)
{
    chomp;

    my ($read, $flag, $chr, $pos, $seq) = (split /\t/)[0 .. 3, 9]; my $len = length $seq;
    $box{$read}{'r2'} = [$chr, $pos, $flag, $len];
}
close IN;


my @insertSize;
foreach (keys %box)
{
    if (defined $box{$_}{'r1'} and defined $box{$_}{'r2'})
    {
        next if $box{$_}{'r1'} -> [0] ne $box{$_}{'r2'} -> [0];

        next if $box{$_}{'r1'} -> [2] =~ /r/ and $box{$_}{'r2'} -> [2] =~ /r/;
        next if $box{$_}{'r1'} -> [2] !~ /r/ and $box{$_}{'r2'} -> [2] !~ /r/;

        next if $box{$_}{'r1'} -> [2] =~ /r/ and $box{$_}{'r1'} -> [1] < $box{$_}{'r2'} -> [1];
        next if $box{$_}{'r2'} -> [2] =~ /r/ and $box{$_}{'r2'} -> [1] < $box{$_}{'r1'} -> [1];

        push @insertSize, ($box{$_}{'r1'} -> [1] + $box{$_}{'r1'} -> [3] - 1) - $box{$_}{'r2'} -> [1] + 1 if $box{$_}{'r1'} -> [2] =~ /r/;
        push @insertSize, ($box{$_}{'r2'} -> [1] + $box{$_}{'r2'} -> [3] - 1) - $box{$_}{'r1'} -> [1] + 1 if $box{$_}{'r2'} -> [2] =~ /r/;
    }
}


(my $stat = Statistics::Descriptive::Full -> new) -> add_data ((sort {$a <=> $b} @insertSize)[(int (@insertSize - 1) * 0.05) .. (int (@insertSize - 1) * 0.95)]);
printf "#LANE\tNUM_of_PAIRS\tMEAN\tMEDIAN\tSD\n%s\t%d\t%d\t%d\t%.2f\n",
       $lane, $stat -> count, $stat -> mean, $stat -> median, $stat -> standard_deviation;
