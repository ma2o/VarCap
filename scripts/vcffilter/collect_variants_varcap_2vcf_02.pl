#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Env;
require "sub_variants_2vcf.pl";
require "sub_variants_ref2vcf.pl";

# updated file that reads variant files from centralized folder
# and incorporates basic tagging [BP,SAR,COV]

### definitions
my $ref_parsed;
my %sum_calls;
my @vcf_calls;
my @vcf_calls_sorted;
my ($path_vcfs,$filename_base,$repeats,$filter_min_calls_per_caller,$filter_min_caller_per_var,$ref_file);
# set minimum read count = minimum absolute abundance (MAA)
my $min_reads = "3";# everything at and below $min_reads will not taken into consideration

# check, if path to variant.config file is defined, else assume it is in the same directory
my $config_file;
if ( @ARGV ){
	$config_file = $ARGV[0]."/variant.config";
} else {
	$config_file = "variant.config";
}

# define hash of arrays for caller data
# caller => [path to caller,[suffix(es) of result files]]
my %caller_data = (
	samtools => ["samtools",[".snp.flt.vcf"]],
	gatk => ["gatk",[".snps.filtered.vcf"]],
	varscan => ["varscan",[".filter.snp",".filter.indel"]],
	lofreq => ["lofreq",["_lofreq_filter"]],
	lofreq2 => ["lofreq2",["_lofreq2_filter.vcf"]],
	breakdancer => ["breakd",[".ctx"]],
	delly => ["delly",[".del.txt",".dup.txt",".inv.txt",".jmp.txt"]],
	delly_072 => ["delly_072",[".del.vcf",".dup.vcf",".ins.vcf",".tra.vcf",".inv.vcf"]],
	pindel => ["pindel",["_D.vcf","_SI.vcf","_INV.vcf","_TD.vcf","_LI","_BP"]],
	pindel025 => ["pindel025",["_ALL.vcf"]],
	cortex => ["cortex",["_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf"]],
);


##### start block
define_ext_variables();
print "Load reference\n";
$ref_parsed = parse_reference();
print "collecting variant calls from input files\n";
caller_iteration();
# print Dumper(\%sum_calls);
print "creating vcf structure\n";
convert_to_vcf_raw();
# print Dumper(\@vcf_calls);
print "writing vcf entries\n";
generate_vcf_files(@vcf_calls_sorted);



### subroutines

sub define_ext_variables {
	# define perl variables from bash export(ed) bash variables
	$path_vcfs = $ENV{PATH_CALLER_COLLECT};
	$ref_file = $ENV{REF_FA};
	
	# $path_data = $ENV{PATH_DATA};
	$filename_base = $ENV{BAM_NAME_BASE};
	$repeats = $ENV{REPEATS};
	$filter_min_calls_per_caller = $ENV{MIN_CPC};
	$filter_min_caller_per_var = $ENV{MIN_CPV};
}

# iterates over all caller available in %caller_data
sub caller_iteration {
	foreach my $caller (keys %caller_data){
		filename_iteration($caller);
	}
	
}

# open and prepare reference file for parsing
sub parse_reference {
	# first open and parse reference file (remove newline)
	open(my $ref_handle, "<", $ref_file) || die "cannot open $ref_file: $!\n";
	my @ref_a = <$ref_handle>;
	chomp(@ref_a);
	my $ref_load = $ref_a[0]."\n";# add my
	for(my $i = 1; $i < @ref_a; $i++) {
		$ref_load = $ref_load.$ref_a[$i];
	}
	return($ref_load);
}

sub filename_iteration {
	my ($caller) = @_;
	
	# first get path and construct filename
	my @caller_suffixes = @{$caller_data{$caller}[1]};
	my $caller_name = $caller_data{$caller}[0];
	for my $caller_suffix (@caller_suffixes){		
		# print "CSUFF".$caller_suffix."\n";
	    
	    for (my $i = 1; $i <= $repeats; $i++){
		  # my $iter_name = $caller_name."_"."v".$i.$filename_base_mod."v".$i.$caller_suffix;
		  my $iter_name = $caller_name."_"."v".$i;
		  #print $iter_name."\n";
		  my $infile = join "/", $path_vcfs,$iter_name;
		  # check with pattern matching
		  my @input_files = glob("$infile*$caller_suffix");
		  # print "$infile $caller_suffix\n";
		  if ( scalar @input_files == 0){
		  	print "File not found: $infile $caller_suffix - ignoring!\n";
		  	next;
		  }
		  if ( scalar @input_files > 1){
		  	print "Fileprefix or suffix not unique: $infile $caller_suffix - ignoring!\n";
		  	next;
		  }
		  print "Input file:".$input_files[0]."\n";
		  $infile=$input_files[0];
		  # check if file exists
		  if (-e $infile ){
		  	
		    if ($caller eq "samtools"){
		  	  #print "ex_samtools: \n";
		  	  extract_samtools($infile, $i, $caller, $min_reads);
		    }
		    if ($caller eq "gatk"){
		  	  #print "ex_gatk: \n";
		  	  extract_gatk($infile, $i, $caller, $min_reads);
		    }
		    if ($caller eq "varscan"){
		  	  #print "ex_varscan: \n";
		  	  extract_varscan($infile, $i, $caller, $min_reads);
		    }
		    if ($caller eq "lofreq"){
		  	  #print "ex_lofreq: \n";
		  	  extract_lofreq($infile, $i, $caller, $min_reads);
		    }
		    if ($caller eq "lofreq2"){
		  	  #print "ex_lofreq2: \n";
		  	  extract_lofreq2($infile, $i, $caller, $min_reads);
		    }
		    if ($caller eq "breakdancer"){
		  	  #print "ex_breakdancer: \n";
		  	  extract_breakdancer($infile, $i, $caller, $ref_parsed, $min_reads);
		    }
		    if ($caller eq "delly"){
		  	  #print "ex_delly: \n";
		  	  extract_delly($infile, $i, $caller, $ref_parsed, $min_reads);
		    }
		    if ($caller eq "delly_072"){
		  	  #print "ex_delly: \n";
		  	  extract_delly_072($infile, $i, $caller, $ref_parsed, $min_reads);
		    }
		    if ( $caller eq "pindel" || $caller eq "pindel025"){
		  	  #print "ex_pindel: \n";
		  	  extract_pindel($infile, $i, $caller, $ref_parsed, $min_reads);
		    }
		    if ($caller eq "cortex"){
		  	  #print "ex_cortex: \n";
		  	  extract_cortex($infile, $i, $caller, $min_reads);
		    }
	      }
	    }
		
	}
}
# subs extract_caller are in sub_variants.pl file


# feed data to hashtable
# check if key exists, if not create and add value, else just add value
# structure is: hashtable, hashtable, hashtable, hashtable, array, array
# {$rep_i {$chrom {$position {$caller [[$repeat, $type, $length, $ref, $var,$reads_ref, $reads_alt],[...]]}}}}
sub add_data {
  my($position, $caller, $repeat, $type, $length, $ref, $var, $reads_ref, $reads_alt, $chrom) = @_;
  my @variant_details = ( $repeat, $type, $length, $ref, $var, $reads_ref, $reads_alt );
  # print Dumper($position, $caller, @variant_details);
  # check for file_number (repeat)
  if (exists $sum_calls{$repeat}) {
  	
  	# check for chromosome
  	if (exists $sum_calls{$repeat}{$chrom}) {
    	# check if position key exists, else add/create new
    	if (exists $sum_calls{$repeat}{$chrom}{$position}) {
      	#check if caller for position exists, else create new
      	if (exists $sum_calls{$repeat}{$chrom}{$position}{$caller}) {
    		push @{ $sum_calls{$repeat}{$chrom}{$position}{$caller} },[ @variant_details ];
      	} else {
    		#print "add caller \n";
    		$sum_calls{$repeat}{$chrom}{$position }{$caller} = [[ @variant_details ]];
      	}# caller end
      
    	} else {
      	#print "create new chrom \n";
      	$sum_calls{$repeat}{$chrom}{$position } = {
        	$caller => [[ @variant_details ]],
      	}
    	}# pos end
  	} else {
  		$sum_calls{$repeat}{$chrom}{$position }{$caller} = [[ @variant_details ]];
  	}# chrom end
  
  } else {
  	$sum_calls{$repeat}{$chrom}{$position }{$caller} = [[ @variant_details ]];
  	
  }# repeat end
}

sub convert_to_vcf_raw {
	# for all file iterations
	for my $rep_i (sort {$a <=> $b} keys %sum_calls){
		my @rep_calls;
		
		# for all chromosomes
		for my $chrom (sort {$a cmp $b} keys %{$sum_calls{$rep_i}}){
			# for all positions on each chromosome
			for my $pos (sort {$a<=>$b} keys %{$sum_calls{$rep_i}{$chrom}}){
				
				# summarize calls for same position
				my $id = $rep_i;
				# my $id = ".";
				my $qual = ".";
				my $filter =".";
				my @vcf_pos;
				my $vcf_pos_line;
				my ($sv_len,$sv_type,$sv_info,$ref_i,$alt_i);
				for my $caller (sort {$a cmp $b} keys %{$sum_calls{$rep_i}{$chrom}{$pos}}){
					# for all caller for this position
					foreach my $caller_i (@{$sum_calls{$rep_i}{$chrom}{$pos}{$caller}}){
						my @caller_entry = @{$caller_i};
						# $repeat, $type, $length, $ref, $var, $reads_ref, $reads_alt
						my $type = $caller_entry[1];
						my $length = $caller_entry[2];
						my $ref_base = $caller_entry[3];
						my $alt_base = $caller_entry[4];
						my $ref_reads = $caller_entry[5];
						my $alt_reads = $caller_entry[6];
						# print "$caller, $type, $length, $ref_base, $alt_base, $ref_reads, $alt_reads\n";
						my $sv_caller = $caller.",".$type.",".$length.",".$ref_reads.",".$alt_reads;
						# extract/add type
						$type =~ s/SNP_FROM_COMPLEX/SNPc/;
						$type =~ s/INDEL_FROM_COMPLEX/INDc/;
						$type =~ s/PH_SNPS/PH_S/;
					
						# create info entry
						my $info = "SVTYPE=".$type.";SVLEN=".$length.";SVCALL=".$sv_caller;
						# create format entry
						my $tcov = "-1";
						my $pct = "-1";
						if ($ref_reads eq "nd"){
							$tcov = "nd";
							$pct = "nd";
						} else {
							$tcov = ($ref_reads + $alt_reads);
							$pct = (int(($alt_reads * 100 / $tcov) * 10)) / 10;
						}
						my $format = "TAA:RAA:VAA:VRA\t".$tcov.":".$ref_reads.":".$alt_reads.":".$pct;
						# generate vcf entry for each input file
						$vcf_pos_line = $chrom."\t".$pos."\t".$id."\t".$ref_base."\t".$alt_base."\t".$qual."\t".$filter."\t".$info."\t".$format;
						# print $vcf_pos_line."\n";
						push @rep_calls,$vcf_pos_line;
						
					}# end caller entries
				}# end caller
			}# end pos
		}# end chrom
		push @vcf_calls,[@rep_calls];
	}# end rep_i
	# sort the array by 1, 2, 3 column
	foreach my $input_files (@{vcf_calls}){
		my @file_entry = @{$input_files};
		my @file_entry_sorted = sort { (split '\t', $a)[0] cmp (split '\t', $b)[0] || (split '\t', $a)[1].(split '\t', $a)[2] <=> (split '\t', $b)[1].(split '\t', $b)[2] } @file_entry;
		push @vcf_calls_sorted, [@file_entry_sorted];
	}
}


sub generate_vcf_files {
	my (@vcf_entries) = @_;
	# iterate over array
	my $counter_rep = "1";
	foreach my $call_iteration (@vcf_entries) {
		my @call_iteration_a = @{$call_iteration}; 
		# print header
		my @header_entries = header();
		my @call_entries;
		foreach my $call_line (@call_iteration_a) {
			# print $call_line."\n";
			push @call_entries,$call_line;
		}
		my @header_and_calls = (@header_entries,@call_entries);
		# create filename
		my $filename_total = $filename_base."_".$counter_rep.".vcf";
		# write output(s) per iteration to file
		write_2_file($filename_total,@header_and_calls);
		$counter_rep++;
	}
}

sub write_2_file {
	my ($filename, @vcf_file) = @_;
	my $vcf_content = join("\n",@vcf_file);
	# print $vcf_content;
	open FILE, ">", $filename or die $!;
	print FILE $vcf_content;
	close FILE;
}

sub header {
	#create header for vcf file
	my $vcf_format = "##fileformat=VCFv4.2";
	my $vcf_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	my $info_svcall = "##INFO=<ID=VCALL,Number=1,Type=String,Description=\"Listing of caller,type,length,ref reads, var reads responsible for the variant call\">";
	my $info_svtype = "##INFO=<ID=VTYPE,Number=1,Type=String,Description=\"Type of structural variant\">";
	my $info_svlen = "##INFO=<ID=VLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">";
	my $filter_col = "##FILTER=<ID=COL/COH,Number=1,Type=String,Description=\"Low or high total coverage at the position\">";
	my $filter_sar = "##FILTER=<ID=SAR,Number=1,Type=String,Description=\"SNP accumulating region\">";
	my $filter_bpa = "##FILTER=<ID=BPA,Number=1,Type=String,Description=\"Break position associated\">";
	my $format_taa = "##FORMAT=<ID=TAA,Number=1,Type=String,Description=\"Total absolute abundance (reads)\">";
	my $format_raa = "##FORMAT=<ID=RAA,Number=1,Type=String,Description=\"Reference absolute abundance (reads)\">";
	my $format_vaa = "##FORMAT=<ID=VAA,Number=1,Type=String,Description=\"Variant absolute abundance (reads)\">";
	my $format_vra = "##FORMAT=<ID=VRA,Number=1,Type=String,Description=\"Variant relative abundance (percentage)\">";
	
	my @header_entries = ($vcf_format, $info_svcall, $info_svtype, $info_svlen, $filter_col, $filter_sar, $filter_bpa, $format_taa, $format_raa, $format_vaa, $format_vra, $vcf_header);
	return (@header_entries);
}

# print Dumper(\%sum_calls);
# print Dumper(\%vcf_calls);
