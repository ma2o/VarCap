#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


my @vcf_input;

if (@ARGV ge "2" || die "usage: perl get_coverage.pl <vcf> <variant.coverage> <average.coverage/chrom>  ") {
  chomp (@vcf_input = @ARGV);
}

# input vars
my $vcf_name = $vcf_input[0];
my $cov_name = $vcf_input[1];
my $cov_av = $vcf_input[2];

# main
my %chrom_cov;
if ( defined($cov_av) ){
	%chrom_cov = average_cov($cov_av);
}

vcf_coverage();

# build hashtable for chrom average coverage
sub average_cov {
	my $cov_sub = @_;
	my @cov_average = open_vcf($cov_av);
	my %chrom_cov_sub;
	foreach my $line (@cov_average){
		my @cov_details = split("\t",$line);
		my $chrom = $cov_details[0];
		my $cov = $cov_details[1];
        $chrom_cov_sub{$chrom} = $cov;
	}
	return (%chrom_cov_sub);
}

sub cov_dev {
	my ($cov,$chrom,$dev) = @_;
  	my $av_cov = $chrom_cov{$chrom};
  	$dev = ( $dev/100);
  	my $cov_deviation = ( $cov/$av_cov );
  	my $cov_ret = "NA";
  	if ( $cov_deviation < ( 1 - $dev ) ){
  		$cov_ret = "COL";
  	} elsif ( $cov_deviation > ( 1 + $dev ) ){
  		$cov_ret = "COH";
  	}
  	return ($cov_ret);
}

# run through vcf and replace missing coverage
sub vcf_coverage {
	
	my (@vcf_raw, @cov_pos);
	@vcf_raw = open_vcf($vcf_name);
	@cov_pos = open_vcf($cov_name);
	
	
	foreach my $line (@vcf_raw){
		if ($line !~ m/^#/){
  			# if there is a nd in the end, get coverage and percentage for that position
  			if ($line =~ m/nd$/){
  				my @line_fields = split("\t",$line);
  				$line_fields[9] =~ m/:([0-9]+):nd$/;
  				my $reads = $1;
  				my $chrom = $line_fields[0];
  				my $pos = $line_fields[1];
  				my $regex = $chrom."\t".$pos;
  				$regex =~ tr/\\|/./;
  				# $regex =~ s/([ab])/$map{$1}/g;
  				# print "covtest,".$chrom.",".$pos."\n";
  				# print "regex:".$regex."\n";
  				# get coverage of position
  				my @cov_pos_line = grep(/$regex/, @cov_pos);
  				# print Dumper (@cov_pos_line);
  				# create entry for no coverage, array is zero
  				my $cov_pos_line_size = @cov_pos_line;
  				# print "size:".$cov_pos_line_size."\n";
  				if ($cov_pos_line_size < 1 ){
  					# samtools does not report zero coverage, so mimic zero coverage
  					@cov_pos_line = ("$chrom.\t.$pos.\t0\n");
  				}
  				chomp(@cov_pos_line);
  				my @cov_pos_line_sp = split("\t", $cov_pos_line[0]);
  				my $coverage = $cov_pos_line_sp[2];
  				# for percentage calculation check if coverage is 0, bec division by 0, return 100 percent
  				my $percentage;
  				if ($coverage == 0){
  					$percentage = "100";
  				} else {
  					$percentage = (int( $reads * 1000 / $coverage )) / 10;
  				}
  				# add coverage tag to filter column if coverage deviates more than 20% from the average
  				my $allowed_deviation = "20";
  				my $coverage_deviation = cov_dev($coverage,$chrom,$allowed_deviation);
  				if ( $coverage_deviation eq "COL" || $coverage_deviation eq "COH" ){
  					if ( $line_fields[6] eq "." ){
  						 $line_fields[6] = $coverage_deviation;
  					} else {
  						$line_fields[6] = $line_fields[6].",".$coverage_deviation;
  					}
  				}
  				
  				# update variant base
  				if ($line_fields[4] =~ m/[0-9,nd]+/){
  					$line_fields[4] = ".";
  				}
  				# update genotype info
  				my $genotype_info = $coverage.":"."nd".":".$reads.":".$percentage."\n";
  				pop (@line_fields);
  				push (@line_fields, $genotype_info);
  				# print updated line
  				print join("\t",@line_fields);
  				# print $chrom."\t".$pos."\t".$reads."\n";
  				# print $coverage.":"."nd".":".$reads.":".$percentage."\n";
  			} else {
  				# filter and add cov tag if no mod to total cov is needed
  				# add coverage tag to filter column if coverage deviates more than 20% from the average
  				my @line_fields = split("\t",$line);
  				$line_fields[9] =~ m/^([0-9]+):/;
  				my $coverage = $1;
  				my $chrom = $line_fields[0];
  				my $allowed_deviation = "20";
  				my $coverage_deviation = cov_dev($coverage,$chrom,$allowed_deviation);
  				if ( $coverage_deviation eq "COL" || $coverage_deviation eq "COH" ){
  					if ( $line_fields[6] eq "." ){
  						 $line_fields[6] = $coverage_deviation;
  					} else {
  						$line_fields[6] = $line_fields[6].",".$coverage_deviation;
  					}
  				}
  				print join("\t",@line_fields);
  				# print $line;
  			}
		} else {
			# print header lines
			print $line;
		}
	}
	
}

# open and process reference
sub open_vcf {
	my ($ref_name) = @_;
	my @ref_file;
	open (FILE, "<", $ref_name) || die "File not found: ".$ref_name;
	@ref_file = <FILE>;
	close FILE;
	return (@ref_file);
}
