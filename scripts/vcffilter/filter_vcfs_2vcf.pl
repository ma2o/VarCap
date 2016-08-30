#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Data::Dumper;

# script filters vcf due to absolute read count and percentage of variant reads
if (scalar(@ARGV) == "0"){
	die "Usage: perl filter_vcfs_2vcf.pl <file.vcf> <lower_limit_reads:int> <lower_limit_percentage:int>\n";
}

my $vcf_name = $ARGV[0];
my $base_file = basename($ARGV[0]);
$base_file =~ s/.vcf$//;
my $cutoff_abs = $ARGV[1];
my $cutoff_pct = $ARGV[2];


# main
read_vcf_files($vcf_name);


# read in input files
sub read_vcf_files {
	my ($vcf_filename) = @_;
	my @vcf_file = open_vcf($vcf_filename);
	my $counter_line = "0";
	foreach my $line (@vcf_file){
    		chomp ($line);
    		
    		if ($line !~ /^#/){
    			$counter_line++;
    			my @split_line = split('\t',$line);
    			my $set_filter = "0";# several filters will update this variable in order to determine if the line should be printed
    			# --- filter set ---
    			# filter out all calls within the first 15 bases
    			# if ($split_line[1] < "15"){
    			#	$set_filter = "1";
    			# }
    			# filter breakdancer ITX and delly DUP shorter than 300bp
    			my $info_field = $split_line[7];
    			if ( $info_field =~ m/delly/ && $info_field =~ m/SVLEN=[-]{0,1}[1-2]{1}[0-9]{2}\D/){
    				$set_filter = "1";
    				# print  $line."\n";
    			}
    			# if ( $info_field =~ m/breakdancer/ && $info_field =~ m/SVLEN=[-]{0,1}[1-2]{1}[0-9]{2}/){
    			#	$set_filter = "1";
    			#	# print  $line."\n";
    			# }
    			# filter below MRA threshold
    			my $cov_field = $split_line[9];
    			my @cov_vals = split(":", $cov_field);
    			if ($cov_vals[3] < $cutoff_pct){
    				$set_filter = "1";
    				# print  $line."\n";
    			}
    			# filter below MAA threshold
    			if ($cov_vals[2] < $cutoff_abs){
    				$set_filter = "1";
    				# print  $line."\n";
    			}
    			# filter homopolymers
    			$info_field =~ m/HP_LEN=([0-9]+)/;
    			my $hp_length = $1 || "0";
    			if ($hp_length > 8 && $cov_vals[3] < 5) {
    				$set_filter = "1";
    			}
    			# filter SNPs or PH_S from cortex because of missing strand bias filter
    			if ( $info_field =~ m/SVTYPE=PH_S/){
    				$set_filter = "1";
    			}
    			
    			# tag set
    			# tag repetitive regions
    			
    			
    			# check for filter var, and print if not set
    			if ($set_filter == "0"){
    				print $line."\n";
    			}
    			
    		} else {
    			# print header
    			print $line."\n";
    		}
    }
}

# open and process file
sub open_vcf {
	my ($ref_name) = @_;
	my @ref_file;
	open (FILE, "<", $ref_name) || die "file not found";
	@ref_file = <FILE>;
	close FILE;
	return (@ref_file);
}

