#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Data::Dumper;

# script filters vcf due to absolute read count and percentage of variant reads
if (scalar(@ARGV) == "0"){
	die "Usage: perl filter_homopolymers_2vcf.pl <file.vcf> <reference.fasta>\n";
}

my $vcf_name = $ARGV[0];
my $base_file = basename($ARGV[0]);
$base_file =~ s/.vcf$//;
my $ref_fasta = $ARGV[1];



# main
my @vcf_file = open_file($vcf_name);
my @ref_file = open_file($ref_fasta);
my %ref_lookup = ref_table();
@ref_file = [];
find_homopolymers(5);

# iterate through small indel ans snp positions and find out if they lie within homopolymer stretches
sub find_homopolymers {
	my ($hp_lenght) = @_;
	foreach my $line (@vcf_file){
    	chomp ($line);
    	if ($line !~ m/^#/){
    		my @line_fields = split('\t',$line);
    		if ($line_fields[7] =~ m/(DEL|INS|IND|SNP)/ && $line_fields[7] =~ m/SVLEN=[-]{0,1}[0-5]{1}\D/ ){
    			# get sequence at Indel pos
    			my $sequence = substr($ref_lookup{$line_fields[0]},$line_fields[1],15);
    			# calculate homopolymer length
    			my $starting_base = substr($sequence,0,1);
    			$sequence =~ m/($starting_base+)/;
    			my $seq_length = length($1);
    			# if homopolmer length is at or above 8 bases, write tag to filter column
    			if ( $seq_length >= "8" ){
    			  if ( $line_fields[6] eq "."){
    				  $line_fields[6] = "HOP";
    			  } else {
    				  $line_fields[6] = $line_fields[6].",HOP";
    			  }
    			}
    			# update info entry
    			my $info_update = $line_fields[7].";"."HP_SEQ15=".$sequence.";"."HP_LEN=".$seq_length;
    			$line_fields[7] = $info_update;
    			# print updated line
    			print join("\t",@line_fields)."\n";
    			# print $line_fields[1].":".$line_fields[7].":".$sequence.":".$starting_base.":".$seq_length."\n";
    			
    		} else {
    			print $line."\n";
    		}
    		
    	} else {
    		print $line."\n";
    	}
	}
}

# convert reference input to hash lookup
sub ref_table {
	my %ref_table;
	my $chrom_name;
	my $seq = "";
	my $line_counter = "1";
	# init table with first entry
	foreach my $line (@ref_file){
    	chomp ($line);
    	if ($line =~ m/^>/){# chromosome name
    		# $line =~ m/^>([\w,\d,\.]+)[ ]*/;
                $line =~ m/^>([\w,\d,\.,\|]+)[ ]*/;
    		$chrom_name = $1;
    		$seq = "";
    		# print $chrom_name."\n";
    	} else {
    		if ($seq ne ""){# check if seq exists
    			$seq = $seq.$line;
    		} else {# else initialize it
    			$seq = $line;
    		}
    		$ref_table{$chrom_name} = $seq;
    	}
    	$line_counter++;
	}
	# print Dumper %ref_table;
	return (%ref_table);
}

# open and process file
sub open_file {
	my ($ref_name) = @_;
	my @ref_file;
	open (FILE, "<", $ref_name) || die "file not found";
	@ref_file = <FILE>;
	close FILE;
	return (@ref_file);
}
