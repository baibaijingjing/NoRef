#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use POSIX qw (floor ceil);

my $usage = <<__EOUSAGE__;


########################################################################

 usage: $0 [opts] cluster_1 cluster_2 ...

########################################################################
#
#  --width_per_plot <int>      default: 800
#  --height_per_plot <int>     default: 600
#  --title <str>               default: Genes
#
#  --plot_format <str>         default: png     png | bmp | jpeg | tiff
#
########################################################################

__EOUSAGE__

    ;


my $width_per_plot = 800;
my $height_per_plot = 600;

my $plot_format = 'png';
my $title = 'Genes';

my $help_flag;

&GetOptions( 'h' => \$help_flag,
             'width_per_plot=i' => \$width_per_plot,
             'height_per_plot=i' => \$height_per_plot,
             'plot_format=s' => \$plot_format,
             'title=s' => \$title,
             );


my @cluster_files = @ARGV;

if ($help_flag) {
    die $usage;
}

if ($plot_format ne 'png' and $plot_format ne 'bmp' and $plot_format ne 'jpeg' and $plot_format ne 'tiff') {
    die $usage;
}

unless (@cluster_files) {
    die $usage;
}

main: {

    ## ensure each cluster file can be found as a file
    foreach my $file (@cluster_files) {
        unless (-s $file) {
            die "Error, cannot find file \"$file\" ";
        }
                
    }


    my $R_script = "__tmp_plot_clusters.R";

    open (my $ofh, ">$R_script") or die "Error, cannot write to $R_script";
    print $ofh "files = c(\"" . join("\",\"", @cluster_files) . "\")\n";

    print $ofh "for (i in 1:length(files)) {\n";
#   print $ofh "    png(file= paste(files[i], '.png ', sep=''), width = $width_per_plot, height = $height_per_plot, )\n";
    print $ofh "    $plot_format(file= paste(files[i], '.$plot_format', sep=''), width = $width_per_plot, height = $height_per_plot, res= $width_per_plot * 100 / 800, units='px', )\n";
    print $ofh "    data = read.table(files[i], header=T, row.names=1, check.names=F)\n";
    print $ofh "    ymin = min(data); ymax = max(data);\n";
    print $ofh "    plot_label = paste( length(data[,1]), \" $title\", sep='')\n";
    print $ofh "    op <- par(mar = c(max(nchar(colnames(data)))/2+3, 4, 4, 2)); on.exit(par(op))\n";
    print $ofh "    plot(as.numeric(colMeans(data)), type='o', ylim=c(ymin,ymax), main=\"\", col='black', xaxt='n', xlab='', ylab='')\n";
    print $ofh "    title(main = list(plot_label, cex=2.25) )\n";
    print $ofh "    title(ylab = list('Expression Level', cex=1.35))\n";
    print $ofh "    axis(side=1, at=1:length(data[1,]), labels=colnames(data), las=2, cex.axis = 1.25)\n";
#    print $ofh "    for(r in 2:length(data[,1])) {\n";
#    print $ofh "        points(as.numeric(data[r,]), type='l', col='lightgray')\n";
#    print $ofh "    }\n";
#    print $ofh "    points(as.numeric(colMeans(data)), type='o', col='blue')\n";
    print $ofh "    dev.off()\n";
    print $ofh "}\n";
    
    close $ofh;
    

    &process_cmd("R --vanilla -q < $R_script");
    system "rm $R_script";

    exit(0);
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}
