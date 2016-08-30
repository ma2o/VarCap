#!/usr/local/bin/perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

my @vcf_input;

if (@ARGV eq "2" || die "usage: perl vcf_contains_pos2.pl <file.vcf> <ref.vcf>") {
  chomp (@vcf_input = @ARGV);
}

# input vars
my $vcf1_name = $vcf_input[0];
my $vcf2_name = $vcf_input[1];

# define vars
my (@vcf_1,%vcf_2);
my $header;

### main
@vcf_1 = open_vcf($vcf1_name);
%vcf_2 = read_vcf_files($vcf2_name,2);
# print Dumper %vcf_2;
# print "compare and update vcfs\n";
compare_vcfs();


### subroutines start

# check if pos(and length) of first vcf id contained in pos(and length) of second vcf file
sub compare_vcfs {
	# iterate through vcf entries and tag(update them) as repetitive when found within the second vcf lookup table
	foreach my $line (@vcf_1){

			if ($line !~ m/^#/){
				chomp($line);
				my @splitline = split("\t",$line);
				my $chrom = $splitline[0];
				my $pos = $splitline[1];
				my $vcf_filter = $splitline[6];
				my $vcf_info = $splitline[7];
				# get length
				$vcf_info =~ s/"non"/0/g;
				$vcf_info =~ m/SVLEN=[-]{0,1}([0-9]*)\W/;
				my $info_length = $1 || "0";
				# restrict length of first file to avoid DUP/ITX lengths
				if( $info_length > "3000"){
					$info_length = "3000";
				}
				# print "chrom:".$chrom."pos:".$pos."length=".$info_length."\n";# print
				# print Dumper \%vcf_2;
				
				# use subroutine to check if pos+length lies within a repetitive region:
				# iterate through keys and lengths of lookup table vcf_2, to identify repeat overlaps
				my ($check_rep,$check_info) = check_pos_rep($chrom,$pos,$info_length);
				if ($check_rep eq "REP"){
					# update filter and info column of splitline array
					if ( $splitline[6] eq "."){
						$splitline[6] = $check_rep;
					} else {
						$splitline[6] = $splitline[6].",".$check_rep;
					}
					
					$splitline[7] = $vcf_info.";".$check_info;
				}
				print join("\t",@splitline)."\n";
				
			} else {
				# print header
				# $header = $header.$line;
				# print $line;
			}
	}
}

sub check_pos_rep {
	my ($chrom,$pos,$info_length) = @_;
	my $check_rep = "NO";
	my $check_info;
	
	if (exists ($vcf_2{$chrom})) {
	# if chromosome is present, loop through keys to find overlap
		foreach my $key_2 ( keys %{$vcf_2{$chrom}} ){
			my $val_2 = $vcf_2{$chrom}{$key_2};
			my @val_2a = @{$val_2};
			my $length2 = "0";
			my $rep_entry;
			foreach my $val_2a_entry (@val_2a){
				# get length
				$val_2a_entry =~ m/LENGTH=\W*([0-9]*)\W/;
				my $val_2a_length = $1;
				$val_2a_entry =~ m/POS=([0-9]*)\W/;
				my $val_2a_pos = $1;
				if ( $val_2a_length > $length2 ){
					$length2 = $val_2a_length;
				}
				# check for each rep pos if actual pos lies within and update REPEAT tag
				if ( ($pos <= $key_2 && $key_2 <= ($pos + $info_length)) || ($key_2 <= $pos && $pos <= ($key_2 + $length2)) ){
					# overlap: write chrom,pos
					$check_rep = "REP";
					# update info rep entry
					if ($check_info){
						$check_info = $check_info.",".$val_2a_pos.":".$val_2a_length;
					} else {
						$check_info = "REPEAT=".$val_2a_pos.":".$val_2a_length;
					}
				}
			}
		}
	}
	return ($check_rep,$check_info);
}

# open and process reference
sub read_vcf_files {
	my ($vcf_name,$nr) = @_;
	my %vcf_x;
	open(FILE, "<", $vcf_name) || die "cannot open $vcf_name: $!\n";
	my @vcf = <FILE>;
	close FILE;
	# create hashtable for chrom
	foreach my $line (@vcf){
		if ($line !~ m/^#/){
			chomp($line);
			my @splitline = split("\t",$line);
			my $chrom = $splitline[0];
			my $pos = $splitline[1];
			# create hashtable for chrom and pos including array entries of the line content (%chrom{%pos{@line_infos}}) 
			my @vcf_xm;
			if( exists $vcf_x{$chrom}){
				if(exists $vcf_x{$chrom}{$pos}){
					push @{$vcf_x{$chrom}{$pos}},$line;
				} else {
					push (@vcf_xm,$line);
					$vcf_x{$chrom}{$pos} = [ @vcf_xm ];
				}
			} else {
				push (@vcf_xm,$line);
				$vcf_x{$chrom} = {
					$pos => [ @vcf_xm ],
				}
			}
			
		}
	}
	return (%vcf_x);
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
