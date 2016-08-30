#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

my @vcf_input;

if ( @ARGV == "2" || die "Usage: perl compare_vcfs_stats.pl <ref vcf1> <alt vcf2> <threshhold>\n"){
	chomp (@vcf_input = @ARGV);
}

my $base_file1 = basename($vcf_input[0]);
$base_file1 =~ s/.vcf$//;
my $base_file2 = basename($vcf_input[1]);
$base_file2 =~ s/.vcf$//;
my $out_both = ${base_file1}."_".${base_file2}.".true_pos.vcf";
my $false_neg = ${base_file1}."_".${base_file2}.".false_neg.vcf";
my $false_pos = ${base_file1}."_".${base_file2}.".false_pos.vcf";
my $thresh = $vcf_input[2];

my $readlength = "100";
my $insertsize = "300";

my (@both,@fp,@fn);
my @output_files = ($out_both,$false_neg,$false_pos);

my %vcf_1;# reference vcf
my %vcf_2;# variant vcf

### main
read_vcf_files($vcf_input[0],1);
read_vcf_files($vcf_input[1],2);
# print "compare vcfs\n";
compare_vcfs();
print_vcf_files($out_both,join("\n",@both));
print_vcf_files($false_neg,join("\n",@fn));
print_vcf_files($false_pos,join("\n",@fp));

### subsroutines

# compare vcf files and write to arrays: @both,@fn,@fp
sub compare_vcfs {
	# extract position and construct hash table
	my $counter_ref_SNP = "0";
	my $counter_ref_Indels = "0";
	keys %vcf_1;
	while(my($key_1, $val_1) = each %vcf_1){
		### prepare reference vcf
		my $pos_ref = $key_1;
		# get type and length for key1(position1) from reference vcf
		# in order to compare: position, type, length
		my @val_1_array = @{$val_1};
		my @val_1_ref = split("\t",$val_1_array[0]);
		my @val_1_info = split(";",$val_1_ref[7]);
		# info column converted to array, convert to hash
		my %info;
		foreach my $info_entry (@val_1_info){
			my @key_val = split("=",$info_entry);
			$info{$key_val[0]} = $key_val[1];
		}
		my $type_ref = $info{"TYPE"};
		my $length_ref = $info{"LENGTH"};
		if ($type_ref eq "SNP"){
			$length_ref = "0";
		}
		my $type_sv = "none";
		if (exists $info{"SVINFO"}) {
			my @svinfo = split(",",$info{"SVINFO"});
			$type_sv = $svinfo[0];
		}
		# print $pos_ref." ".$type_ref." ".$length_ref." ".$type_sv."\n";
		
		# set different offsets according to variant type, the default (also used for SNPs is 0,0), for Indels and SVs see settings below
		my $offset_low = "0";
		my $offset_high = "0";
		if ($type_ref eq "DEL" || $type_ref eq "INS" ){
			# Indels
			if ($type_sv eq "none"){
				# small Indels <= 10
				if ($length_ref <= "10"){
					$offset_low = "5";
					$offset_high = "5";
				} else {
					# large Indels	
					$offset_low = ($readlength * 2);
					$offset_high = ($length_ref + $readlength);
				}
			}
			# ITX
			if ($type_sv eq "ITX"){
				$offset_low = ($readlength * 2);
				$offset_high = ($readlength * 2);
			}
			# DUP
			if ($type_sv eq "DUP"){
				$offset_low = ($readlength * 2);
				$offset_high = ($readlength * 2);
			}
		}
		# INV
		if ($type_ref eq "INV"){
			$offset_low = ($insertsize / 2);
			$offset_high = ($offset_low + $length_ref);
		}
		
		# compare the position,type,length with the variant vcf file positions,type,length with possible offset
		my ($pos_min,$pos_max);
		$pos_min = ($pos_ref - $offset_low);
		$pos_max = ($pos_ref + $offset_high);
		# print $pos_ref." ".$type_ref." ".$length_ref." ".$pos_min." ".$pos_max."\n";
		# iterate through variant positions with different offsets
		for (my $i = $pos_min;$i <= $pos_max;$i++){
			if($vcf_2{$i}){
				# search for position and check for type and length
				# iterate through array
				my @pos_calls = @{$vcf_2{$i}};
				# print Dumper @pos_calls;
				my $j = "0";
				foreach my $pos_call (@pos_calls){
					# print Dumper $pos_call;
					# for each call at position get type and length
					my @pos_call_alt = split("\t",$pos_call);
					# print Dumper @pos_call_alt;
					my @pos_call_info = split("\;",$pos_call_alt[7]);
					# info column converted to array, convert to hash
					my %info_var;
					foreach my $info_entry (@pos_call_info){
						my @key_val_alt = split("=",$info_entry);
						$info_var{$key_val_alt[0]} = $key_val_alt[1];
					}
					my $type_alt = $info_var{"SVTYPE"};
					my $length_alt = $info_var{"SVLEN"};
					if ($length_alt ne "non"){
						$length_alt = abs($length_alt);
					}
					my $svcall_alt = $info_var{"SVCALL"};
					my @svcall_alt_a = split(",",$svcall_alt);
					my $svcaller = $svcall_alt_a[0];
					# get length difference of ref and alt bases
					my $ref_base_len = length $pos_call_alt[3];
					my $alt_base_len = length $pos_call_alt[4];
					my $ref_alt_diff = "0";
					$ref_alt_diff = $alt_base_len - $ref_base_len;
					if ( $type_alt eq "SNP" &&  $ref_alt_diff != 0 ){
						$length_alt = $ref_alt_diff;
						if ( $ref_alt_diff > 0 ){
							$type_alt = "INS";
						} else {
							$type_alt = "DEL";
						}
					}
					
					# compare type and length for ref and alt
					# and categorize variants as tp, fn, fp
					# normal variant type comparison
					if ($type_ref eq $type_alt || $type_sv eq $type_alt){
						# if types match, add them to true positives (both)
						my @tp_a = ($pos_ref,$type_ref,$type_sv,$length_ref,$i,$type_alt,$length_alt,$svcaller);
						my $tp_s = join("\t",@tp_a);
						@both = (@both,$tp_s);
						# delete key,value from vcf1(reference) since variant was found
						delete $vcf_1{$pos_ref};
					# add IND to INS/DEL
					} elsif ((($type_ref eq "INS" &&  $length_ref <= "50") || ($type_ref eq "DEL" &&  $length_ref <= "50")) && ($type_alt eq "IND" || $type_alt eq "RPL")){
						# if types match, add them to true positives (both)
						my @tp_a = ($pos_ref,$type_ref,$type_sv,$length_ref,$i,$type_alt,$length_alt,$svcaller);
						my $tp_s = join("\t",@tp_a);
						@both = (@both,$tp_s);
						# delete key,value from vcf1(reference) since variant was found
						delete $vcf_1{$pos_ref};
					# add gatk INS/DEL
					} elsif ((($type_ref eq "INS" &&  $length_ref <= "3000") || ($type_ref eq "DEL" &&  $length_ref <= "3000")) && ($type_alt eq "DEL" && $svcaller eq "gatk" || $type_alt eq "INS" && $svcaller eq "gatk")){
						# if types match, add them to true positives (both)
						my @tp_a = ($pos_ref,$type_ref,$type_sv,$length_ref,$i,$type_alt,$length_alt,$svcaller);
						my $tp_s = join("\t",@tp_a);
						@both = (@both,$tp_s);
						# delete key,value from vcf1(reference) since variant was found
						delete $vcf_1{$pos_ref};
					# add ITX/DUP
					} elsif (($type_sv eq "DUP" || $type_sv eq "ITX") && ($type_alt eq "DUP" || $type_alt eq "ITX")){
						# if types match, add them to true positives (both)
						my @tp_a = ($pos_ref,$type_ref,$type_sv,$length_ref,$i,$type_alt,$length_alt,$svcaller);
						my $tp_s = join("\t",@tp_a);
						@both = (@both,$tp_s);
						# delete key,value from vcf1(reference) since variant was found
						delete $vcf_1{$pos_ref};
					# add LI (large insertion from pindel)
					} elsif (($type_ref eq "INS") && ($type_alt eq "LI")){
						# if types match, add them to true positives (both)
						my @tp_a = ($pos_ref,$type_ref,$type_sv,$length_ref,$i,$type_alt,$length_alt,$svcaller);
						my $tp_s = join("\t",@tp_a);
						@both = (@both,$tp_s);
						# delete key,value from vcf1(reference) since variant was found
						delete $vcf_1{$pos_ref};
					# add BP (break position from pindel)
					} elsif (($type_ref eq "INS" || $type_ref eq "DEL" || $type_ref eq "INV") && ($type_alt eq "BP")){
						# if types match, add them to true positives (both)
						my @tp_a = ($pos_ref,$type_ref,$type_sv,$length_ref,$i,$type_alt,$length_alt,$svcaller);
						my $tp_s = join("\t",@tp_a);
						@both = (@both,$tp_s);
						# delete key,value from vcf1(reference) since variant was found
						delete $vcf_1{$pos_ref};
					# add COMPLEX call from cortex ()
					} elsif (($type_ref eq "INS" || $type_ref eq "DEL" || $type_ref eq "INV") && ($type_alt eq "COMPLEX")){
						# if types match, add them to true positives (both)
						my @tp_a = ($pos_ref,$type_ref,$type_sv,$length_ref,$i,$type_alt,$length_alt,$svcaller);
						my $tp_s = join("\t",@tp_a);
						@both = (@both,$tp_s);
						# delete key,value from vcf1(reference) since variant was found
						delete $vcf_1{$pos_ref};
					} else {
						
						# TODO add special cases for IND,LI,BP
						
						# if found at same position but different type, call and add as fp
						@fp = (@fp,$pos_call);
					}
					
				}
				# remove entry from variant hash vcf2
				delete $vcf_2{$i};
			}
		}# variant iteration end

	}
	# post processing of tables
	# true positives(@both) is ready
	
	# false negatives are left in %vcf_1, so convert values to array
	foreach my $val (values %vcf_1){
		my @val_a = @{$val};
		@fn = (@fn,@val_a);
		# push(@fn,$val);
	}
	# convert remaining entries of vcf2 to array (assumed as fp) and add to previous recorded fp
	my @vcf2_fp;
	foreach my $val (values %vcf_2){
		my @val_a = @{$val};
		@vcf2_fp = (@vcf2_fp,@val_a);
		# push(@both,$val);
	}
	@fp = (@vcf2_fp,@fp);
	
	my @fn_sort = sort { (split '\t', $a)[1] <=> (split '\t', $b)[1] } @fn;
	my @fp_sort = sort { (split '\t', $a)[1] <=> (split '\t', $b)[1] } @fp;
	my @both_sort = sort { (split '\t', $a)[0] <=> (split '\t', $b)[0] } @both;
	@fn = @fn_sort;
	@fp = @fp_sort;
	@both = @both_sort;

}

# read in input files to hashtable
sub read_vcf_files {
	my ($vcf_name,$nr) = @_;
	my %vcf_x;
	open(FILE, "<", $vcf_name) || die "cannot open $vcf_name: $!\n";
	my @vcf = <FILE>;
	close FILE;
	# create hashtable of pos => vcf_line entries
	foreach my $line (@vcf){
		if ($line !~ m/^#/){
			chomp($line);
			my @splitline = split("\t",$line);
			my $pos = $splitline[1];
			if($vcf_x{$pos}){
				my @vcf_xm = @{$vcf_x{$pos}};
				push (@vcf_xm,$line);
				$vcf_x{$pos} = [@vcf_xm];
			} else {
				my @vcf_xm;
				push (@vcf_xm,$line);
				$vcf_x{$pos} = [@vcf_xm];
			}
		}
	}
	if ($nr == "1"){
		%vcf_1 =  %vcf_x;
	} else {
		%vcf_2 =  %vcf_x;
	}
}

# print to vcf files
sub print_vcf_files {
	my ($file, $body) = @_;
	my $header = "## VCF version 4.1\n##compare_vcfs_stats\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n";
	open FILE, ">", $file or die $!;
	print FILE $header;
	print FILE $body;
	close FILE;
}
