#!/usr/bin/perl
#$-q all.q@cube[ab]*

use strict;
use warnings;
use File::Basename;
use Data::Dumper;

# script filters vcf due to absolute read count and percentage of variant reads
if (scalar(@ARGV) == "0"){
	die "Usage: perl convert_vcf2mysql.pl <file.vcf>\n";
}

my $vcf_name = $ARGV[0];
my $base_file = basename($ARGV[0]);
$base_file =~ s/.vcf$//;

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
    			# update vcf entries, to contain one value per column (splitting info and genotype column)
    			# info field update
    			my $info_field = $split_line[7];
    			my @info_format = reformat_info($info_field);
    			# update snpeff annotation from info field if present
    			my @snpeff_format;
    			if ($info_field =~ "EFF=" || $info_field =~ "ANN="){
    				@snpeff_format = reformat_snpeff($info_field);
    			} else {
    				@snpeff_format = ("NOT_ANNOTATED",".");
    			}
				# gentype field update
    			my $genotype_field = $split_line[9];
    			my @genotype_format = reformat_genotype($genotype_field);
    			# set new format
    			# format: CHROM POS ID REF ALT QUAL FILTER REPEAT_POS SVTYPE SVLEN CALLER TOT_COV REF_COV ALT_COV ALT_PC
    			splice(@split_line,-3);# remove last 3 columns
    			my @reformat_line = (@split_line,@info_format,@genotype_format,@snpeff_format);# add info and genotype array
    			# add experiment/project info column
    			# unshift(@reformat_line,"project_id");
    			# print updated/reformated line
    			print join("\t",@reformat_line)."\n";
    		}
	}
}

# reformat info entries to column view
sub reformat_info {
	my ($info_field) = @_;
	my @info_vals = split(";", $info_field);
    my %info_format;
    foreach my $info_line (@info_vals){
    	if ($info_line =~ m/=/){
    		my @info_2 = split("=", $info_line);
    		$info_format{$info_2[0]} = $info_2[1];
    	} else {
    		$info_format{$info_line} = "nd";
    	}
    }
    # further decompose SVCALL
    my @caller_info = split(",",$info_format{"SVCALL"});
    my $caller = $caller_info[0];
    # add repeat column
    my $repeat = ".";
    if ($info_format{"REPEAT"}){
    	$repeat = $info_format{"REPEAT"};
    }
    # add homopolymer length and sequence
    my $hp_len = ".";
    my $hp_seq = ".";
    if ($info_format{"HP_LEN"}){
    	$hp_len = $info_format{"HP_LEN"};
    }
    if ($info_format{"HP_SEQ15"}){
    	$hp_seq = $info_format{"HP_SEQ15"};
    }
    # bring hashtable to new column format: REPEAT_POS SVTYPE SVLEN CALLER
    my @info_column = ($repeat,$info_format{"SVTYPE"},$info_format{"SVLEN"},$caller,$hp_len,$hp_seq);
    # while (my ($k,$v)=each %info_format){print "$k $v\n"}
    return (@info_column);
}

# reformat genotyping/read coverage
sub reformat_genotype {
	my ($geno_field) = @_;
	my @geno_vals = split(":", $geno_field);
	# format: total_coverage:ref_coverage:var_coverage:percentage_of_variant_to_total
    my @geno_column = @geno_vals;
    return (@geno_column);
}

# reformat snpeff annotation tags within info column
sub reformat_snpeff {
	my ($info_field) = @_;
	my @snpeff_column;
	my @snpeff_locations_temp;
	my @snpeff_gene_temp;
        my @info_vals;
        if ( $info_field =~ "EFF=" ){
	      @info_vals = split("EFF=", $info_field);
	    } else {
          @info_vals = split("ANN=", $info_field);
        }
	my @snpeff_entries = split(",", $info_vals[1]);
	foreach my $snpeff_value (@snpeff_entries){
		# loop through snpeff annotation (most positions will only have one entry)
		my @snpeff_details = split("\\|",$snpeff_value);
		my $snpeff_location = $snpeff_details[1];
		my $snpeff_gene = $snpeff_details[3];
		
		if ( scalar @snpeff_locations_temp < "1"){
			@snpeff_locations_temp = ($snpeff_location);
		} else {
			push(@snpeff_locations_temp,$snpeff_location);
		}
		# add affected gene name(s)
		if ( scalar @snpeff_gene_temp < "1"){
			@snpeff_gene_temp = ($snpeff_gene);
		} else {
			push(@snpeff_gene_temp,$snpeff_gene);
		}
	}
	if ($snpeff_gene_temp[0] eq ""){
		@snpeff_gene_temp = (".");
	}
	@snpeff_column = ( join(",",@snpeff_locations_temp),join(",",@snpeff_gene_temp) );
	return (@snpeff_column);
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
