#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Data::Dumper;

# script tags SNPs of vcf file based on accumulation of number of allowed SNPs per region (normally readlength)
if (scalar(@ARGV) == "0"){
	die "Usage: perl filter_multi_snps_2vcf.pl <file.vcf> <number of allowed SNPs> <region>\n";
}

my $vcf_name = $ARGV[0];
my $base_file = basename($ARGV[0]);
$base_file =~ s/.vcf$//;
my $number_snps = $ARGV[1];
my $length_region = $ARGV[2];


# main

read_vcf_files($vcf_name);


# read in input files
sub read_vcf_files {
	my ($vcf_filename) = @_;
	my @vcf_file = open_vcf($vcf_filename);
	my $counter_line = "0";
	my @region_vars;# stores all variants that span the region
	my %region_vars_pos_unique;# stores uniques for snp positions as key and the times the SNP is encountered at the position as value
	my $region_snps_counter = "0";
	
	# loop through file
	foreach my $line (@vcf_file){
    		chomp ($line);
    		
    		if ($line !~ /^#/){
    			$counter_line++;
    			my @splitline = split('\t',$line);
    			my $chrom = $splitline[0];
				my $pos = $splitline[1];
				my $vcf_filter = $splitline[6];
				my $vcf_info = $splitline[7];
				# get type
				$vcf_info =~ m/SVTYPE=([a-z,A-Z,_]*)\W/;
				my $sv_type = $1;
    			
    			# push entries until readlength distance, array of arrays (pos, line)
    			if (scalar(@region_vars) > "0"){# check if array is initialized, else initialize
    				# print Dumper \@region_vars;
    				my $first_array_pos = $region_vars[0][0];
    				my $first_array_content = $region_vars[0][1];
    				
    				# if new position is outside region, loop as long until new position is within the region
    				# and print entries
    				if ( int($pos - $first_array_pos) > $length_region ){
    					# 1. if variants span area larger than defined region, then check if SNP count for SAR region is fullfilled
    					if ( $region_snps_counter > $number_snps){
    						# if fullfilled, then iterate through region vars array and update all filter tags for SNPS
    						my @region_vars_upd;
    						foreach my $reg_line (@region_vars){
    							my @reg_line_vals = @{$reg_line};		
    							my @reg_line_split = split("\t",$reg_line_vals[1]);
    							my $reg_line_filter = $reg_line_split[6];
    							my $reg_info = $reg_line_split[7];
								# get type
								$reg_info =~ m/SVTYPE=([a-z,A-Z,_]*)\W/;
								my $reg_type = $1;
								if ($reg_type eq "SNP"){
									if ($reg_line_filter eq "."){
										$reg_line_split[6] = "SAR";
    								} else {
    									if ($reg_line_split[6] !~ m/SAR/){
    										$reg_line_split[6] = $reg_line_filter.","."SAR";
    									}
    								}
								}
								my $reg_line_upd = join("\t",@reg_line_split);
								push (@region_vars_upd,[$reg_line_vals[0],$reg_line_upd]);
    						}
    						@region_vars = @region_vars_upd;
    					}# else do nothing, as SNP threshhold not reached
    					
    					# 2. then add new entry
    					my @pos_line = ([$pos,$line]);
    					push (@region_vars,@pos_line);
    					if ($sv_type eq "SNP"){
    						if (exists ($region_vars_pos_unique{$pos})){
    							$region_vars_pos_unique{$pos}++;
    							$region_snps_counter++;
    							$region_vars_pos_unique{$pos} = "1";
    						} else {
    							$region_snps_counter++;
    							$region_vars_pos_unique{$pos} = "1";
    						}
    				    }
    					
    					# 3. write out all pos from region_var array until all variants lie within one region, or one is left
    					
    					until ($pos - $first_array_pos <= $length_region || scalar(@region_vars) == "1"){
    						# if pos range is bigger than region, print out region_var array, but allways keep one (last entry)
    						# decrement SNP counter
    						my $array_pos = $region_vars[0][0];
    						# print "WRITE: array_pos:".$array_pos."region_snps_counter:".$region_snps_counter."\n";
    						my @array_cont = split("\t",$region_vars[0][1]);
    						my $var_info = $array_cont[7];
    						$var_info =~ m/SVTYPE=([a-z,A-Z,_]*)\W/;
							my $reg_type = $1;
							# print "WRITE";
							# print Dumper \%region_vars_pos_unique;
							if ($reg_type eq "SNP"){
								# print "region_snps_counter:".$region_snps_counter."\n";
								if ($region_vars_pos_unique{$array_pos} > "1"){
									$region_vars_pos_unique{$array_pos}--;
								} else {
									if ($region_vars_pos_unique{$array_pos} == "1"){
										$region_snps_counter--;
										delete $region_vars_pos_unique{$array_pos};
									}
								}
								
							}
    						
    						# print top entry line, then delete array entry
    						print join("\t",$region_vars[0][1])."\n";
    						shift(@region_vars);
    						$first_array_pos = $region_vars[0][0];
    					}# until end
    					
    					
    				} else {
    					# if inside region, add enry and increase snps counter (if SNP)
    					my @pos_line = ([$pos,$line]);
    					push (@region_vars,@pos_line);
    					if ($sv_type eq "SNP"){
    						if (exists ($region_vars_pos_unique{$pos})){
    							$region_vars_pos_unique{$pos}++;
    							# print "SVTYPE:".$sv_type.$pos." inc1_array:".$region_vars_pos_unique{$pos}."\n";
    							# print Dumper \%region_vars_pos_unique;
    						} else {
    							$region_snps_counter++;
    							$region_vars_pos_unique{$pos} = "1";
    							# print "SVTYPE:".$sv_type.$pos." set1_array:".$region_vars_pos_unique{$pos}."\n";
    							# print Dumper \%region_vars_pos_unique;
    						}
    				    }
    				}
    				
    			} else {
    				# initialize array
    				@region_vars = ([$pos,$line]);
    				if ($sv_type eq "SNP"){
    					$region_snps_counter++;
    					$region_vars_pos_unique{$pos} = "1";
    					
    				}
    			}
    			
    		} else {
    			# print header
    			print $line."\n";
    		}
    }
    # print out remaining entries in region_vars array
    foreach my $rem_reg_var (@region_vars){
    	my @rem_reg_var_vals = @{$rem_reg_var};
    	print $rem_reg_var_vals[1]."\n";
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

