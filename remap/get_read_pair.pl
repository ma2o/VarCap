#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
# use IO::Zlib;
# use PerlIO::gzip;

my @fastq_input;

if ( scalar @ARGV == "2" || die "usage: perl get_read_pair.pl <fastq1 query> <fastq2 target lib> ") {
  chomp (@fastq_input = @ARGV);
}

# input vars
my $fq = $fastq_input[0];
my $qu = $fastq_input[1];

my %ta_lookup;

# main
my $fq_lookup = create_ta_lookup($fq);
# print Dumper($fq_lookup);
# my $qu_lookup = qu_lookup_list($qu);
qu_lookup_list($qu);


### subroutines
# get list of read names and serch for entry
sub qu_lookup_list {
	my ($qut) = @_;
	my @temp_arr;
	my @temp_reads;
	print "Lookup read names: \n";
	# open file incl gz
	if ($qut =~ /.gz$/) {
		open(FILE, "gunzip -c $qut |") || die "can't open pipe to $qut";
	} else {
		open(FILE, $qut) || die "can't open pipe to $qut";
	}
	
	# iterate through query file
	while (my $qline = <FILE>) {
		chomp $qline;
		# print "$qline\n";
		# dereferencing hash keys
		if ( exists ${$fq_lookup}{$qline} ){
			print "paired read found: $qline\n";
			# extract fastq entry for target read and write to fastq file2/rearange name and write to fastq file2
			my @read_full = @{${$fq_lookup}{$qline}};
			splice @read_full, 2, 0, '+';
			my $read_full_line = join ("\n",@read_full);
			$read_full_line = $read_full_line."\n";
			# print "$read_full_line\n";
			# store reads in temp array and write to file
			push (@temp_reads,$read_full_line);
			if ( scalar @temp_reads > "1000" ){
				write_reads(@temp_reads);
				splice (@temp_reads);
			}
		} else {
			print "paired read not found: $qline\n";
		}
	}
	# print remaining reads
	write_reads(@temp_reads);
}

# get fastq entry and search for hash key
sub qu_lookup_fastq {
	my ($qut) = @_;
	my @temp_arr;
	print "Lookup read names: \n";
	# open file incl gz
	if ($qut =~ /.gz$/) {
		open(FILE, "gunzip -c $qut |") || die "can't open pipe to $qut";
	} else {
		open(FILE, $qut) || die "can't open pipe to $qut";
	}
	# iterate through query file
	while (my $qline = <FILE>) {
		chomp $qline;
		if ( scalar @temp_arr == "4" ){
			# if array is full/one fastq block read, find readname in target fastq
			my @tkey = grep {/$temp_arr[0]/} keys %{$fq_lookup};
			print "name:$temp_arr[0]\tkey=$tkey[0]\n";
			if ( scalar @tkey == "1"){
				print "paired read found: $tkey[0]\n";
				# extract fastq entry for target read and write to fastq file2/rearange name and write to fastq file2
				
			} else {
				print "paired read not found: $temp_arr[0]\n";
			}
			# empty/initialize array 
			splice (@temp_arr,0);
			push (@temp_arr,$qline);
		} else {
			# fill array
			push (@temp_arr,$qline);
			
		}
	}
}

# create lookup table from target fastq
sub create_ta_lookup {
	my ($fqt) = @_;
	my %ta_lookup_temp;
	my $counter = "1";
	my @temp_arr;
	print "Create hashtable ... ";
	# open file incl gz
	if ($fqt =~ /.gz$/) {
		open(FILE, "gunzip -c $fqt |") || die "can't open pipe to $fqt";
	} else {
		open(FILE, $fqt) || die "can't open pipe to $fqt";
	}
	# iterate through fastq in 4 line blocks (temp stored in array)
	while (my $fline = <FILE>) {
		chomp $fline;
		# if ($counter >= "99"){
		#	last;
		# }
		# print Dumper ($fline);
		# add 4 fastq values to array and add it to hashtable
		if ( scalar @temp_arr == "4" ){
			# print "array size:".scalar @temp_arr."\n";
			# create substring of name
			my $name_sub = ( split / /, $temp_arr[0], 2 )[0];
			# print "name_sub:$name_sub\n";
			$ta_lookup_temp{$name_sub} = [ ($temp_arr[0],$temp_arr[1],$temp_arr[3]) ];# [ splice (@temp_arr,1) ]
			splice (@temp_arr,0);
			push (@temp_arr,$fline);
		} else {
			push (@temp_arr,$fline);
			
		}
		$counter++;
    }
	# print Dumper(%ta_lookup_temp);
	print "done.\n\n";
	close FILE;
	return (\%ta_lookup_temp);
}

# write to file
sub write_reads() {
	my (@reads) = @_;
	my $outfile = "extracted_reads.fastq";
	open (FILE1, ">> $outfile") || die "problem opening $outfile\n";
	print FILE1 @reads;
	close(FILE1);
}

