#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Data::Dumper;

# add tags for numberof callers per position + add margin for SV
if (scalar(@ARGV) == "0"){
	die "Usage: perl filter_caller2pos.pl <file.vcf> <number of callers/variant>\n";
}

my $vcf_name = $ARGV[0];
my $base_file = basename($ARGV[0]);
$base_file =~ s/.vcf$//;
my $numberCaller = $ARGV[1];
my $differenceSV = "150";
my $snpBp_dist = "100";


# main
my @header;
my @input_lines = read_vcf_files($vcf_name);
# print "input".scalar @input_lines;
# print Dumper(@input_lines);
my @pos_exact = comp_pos_exact($numberCaller,$differenceSV);
# print "exact".scalar @pos_exact;
# print Dumper(@pos_exact);
my @pos_approx = comp_pos_approx($numberCaller,$differenceSV);
# print "approx".scalar @pos_approx;
# print Dumper(@pos_approx);
my @pos_bp_sar = snp_bp_match($snpBp_dist);
## print "pos_bp".scalar @pos_bp_sar."\n\n";
## print Dumper(@pos_bp_sar);
print_header(@header);
print_array(@pos_approx);

### subs ###

# gets chrom, pos, type, len and vartype of vcf arrayline
sub var_line {
  	my @line = @_;
  	# set new position values
	my $n_chrom = $line[0];
	my $n_pos = $line[1];
	$line[7] =~ /SVTYPE=(\w+)/;
	my $n_type = $1;
	$line[7] =~ /SVLEN=-{0,1}(\d+)/;
	my $n_len = $1;
	$line[7] =~ /SVCALL=([A-Z0-9a-z_]+),/;
	my $n_call = $1;
	# print "n".$n_chrom.",".$n_pos.",".$n_call."\n";
	# assign vartype: either exact for SNPs and small Indels or approx for SV and large Indels
	my $vartype = "exact";
	if ($n_type =~ /ITX|DUP|LI|BP|INV|COMPLEX/ || ($n_type =~ /DEL|INS/ && $n_len > "10" ) ){
		$vartype = "approx";
	}
	return($n_chrom,$n_pos,$n_type,$n_len,$vartype,$n_call);
}

# gets chrom, pos, type, len and vartype of temp array
sub var_temp {
  	my ($i,@temp_pos) = @_;
  	# set current position values
	my $o_chrom = $temp_pos[$i][0];
	my $o_pos = $temp_pos[$i][1];
	$temp_pos[$i][7] =~ /SVTYPE=(\w+)/;
	my $o_type = $1;
	$temp_pos[$i][7] =~ /SVLEN=-{0,1}(\d+)/;
	my $o_len = $1;
	$temp_pos[$i][7] =~ /SVCALL=([A-Z0-9a-z_]+),/;
	my $o_call = $1;
	# print "o".$o_chrom.",".$o_pos."\n";
	# assign vartype: either exact for SNPs and small Indels or approx for SV and large Indels
	my $vartype = "exact";
	if ($o_type =~ /ITX|DUP|LI|BP|INV|COMPLEX/ || ($o_type =~ /DEL|INS/ && $o_len > "10" ) ){
		$vartype = "approx";
	}
	return($o_chrom,$o_pos,$o_type,$o_len,$vartype,$o_call);
}

# compare callers at the exact same position
sub comp_pos_exact {
	my ($nrCall,$diffSV) = @_;
	my @temp_pos;
	my @output_exact;
	# my @temp_sv;
	# my (@exact,@approx);
	
	foreach my $line (@input_lines){
		# check if data exists
		my @line = @{$line};
		# add to header
		if ($line =~ /^#.*/){
			push (@header,[@line]);
			next;
		}
		if ( exists $temp_pos[0]){
			# get values from line (new)
			my ($n_chrom,$n_pos,$n_type,$n_len,$n_vartype,$n_call) = var_line(@line);
			my $first_val = "0";
			
			# get values from array (temp store)
			my ($o_chrom,$o_pos,$o_type,$o_len,$o_vartype,$o_call) = var_temp($first_val,@temp_pos);
			# print "o".$o_chrom.",".$o_pos."\n";
			
			
			# compare if positon is the same, if then add value
			if ($o_chrom eq $n_chrom && $o_pos == $n_pos && $n_type ne "SNP" && $o_type ne "SNP" ){
				push(@temp_pos,$line);
			} elsif ($o_chrom eq $n_chrom && $o_pos == $n_pos && $n_type eq "SNP" && $o_type eq "SNP" ){
				push(@temp_pos,$line);
			} elsif ( ($o_chrom eq $n_chrom && $o_pos < $n_pos) || ($o_chrom ne $n_chrom) ){
				# if new pos is greater, write temp array and add new pos
				# print Dumper(@temp_pos);
				my $callCounts = @temp_pos;
				foreach my $line2 (@temp_pos){
					# create filter tag and modify filter column
					my $filterCallPerPos = "CV".$callCounts;
					my @line2 = @{$line2};
					if ( $line2[6] eq "."){
						$line2[6] = $filterCallPerPos;
					} else {
						$line2[6] = $line2[6].",".$filterCallPerPos;
					}
					$line2 = join("\t",@line2);
					push (@output_exact,[@line2]);
					# print $callCounts.",".$o_pos.",".$o_vartype."\n";
					# shift @temp_pos;
				}
				splice (@temp_pos, 0, $callCounts);
				push(@temp_pos,$line);
			}
			
		} else {
			# initialize/read values from array
			my ($n_chrom,$n_pos,$n_type,$n_len,$vartype,$n_call) = var_line(@line);
			push(@temp_pos,$line);
			
		}

	}
	# print remaining lines
	my $callCountsRem = @temp_pos;
	foreach my $line3 (@temp_pos){
		my $filterCallPerPos = "CV".$callCountsRem;
		my @line3 = @{$line3};
		if ( $line3[6] eq "."){
			$line3[6] = $filterCallPerPos;
		} else {
			$line3[6] = $line3[6].",".$filterCallPerPos;
		}
		$line3 = join("\t",@line3);
		# print "add_test".$line3."\n";
		push (@output_exact,[@line3]);
		# print $callCountsRem.",".@{$line3}[1]."\n";
	}
	return(@output_exact);
}

# compare callers within a window e.g. structural variations
sub comp_pos_approx {
	my ($nrCall,$diffSV) = @_;
	my @temp_pos;
	my @output_approx;
	my $svCount = 1;
	
	foreach my $line (@pos_exact){
		# check if data exists
		my @line = @{$line};
		my ($n_chrom,$n_pos,$n_type,$n_len,$n_vartype,$n_call) = var_line(@line);
		if ( $n_vartype eq "approx"){
		# check if temp array is initialized/empty
		if ( exists $temp_pos[0]){
			my $first_val = "0";
			# get index of last array entry
			my $ar_size = $#temp_pos;
			# get values from array (temp store)
			my ($o_chrom,$o_pos,$o_type,$o_len,$o_vartype,$o_call) = var_temp($ar_size,@temp_pos);			
			
			# compare if positon is within range, if then add value to temp_array
			if ($o_chrom eq $n_chrom && $o_pos >= $n_pos - $diffSV){
				push(@temp_pos,$line);
			# if position is greater than the region or an another chromosome
			} elsif ( ($o_chrom eq $n_chrom && $o_pos < $n_pos - $diffSV) || ($o_chrom ne $n_chrom) ){
				# if new pos is greater, write temp array and add new pos
				# print Dumper(@temp_pos);
				
				# get number of hits / variant window
				my $callCounts = @temp_pos;
				# get different callers / window
				my %callerTemp;
				my $callerCount = "0";
				foreach my $line2 (@temp_pos){
					my @line2 = @{$line2};
					# compare if 2 different callers call this structural variant
					my ($n_chrom,$n_pos,$n_type,$n_len,$n_vartype,$n_call) = var_line(@line2);
					$callerTemp{$n_call} = $n_type;
				}
				$callerCount = keys %callerTemp;
				
				# modify vcf entry
				foreach my $line2 (@temp_pos){
					my @line2 = @{$line2};
					# create filter tag counts
					my $filterCallPerPos = "CSV".$callCounts;
					# create filter tag caller number
					my $filterCaller = "CNSV".$callerCount;
					# modify filter column of line
					if ( $line2[6] eq "."){
						$line2[6] = $filterCallPerPos.",".$filterCaller.",SVID".$svCount;
					} else {
						$line2[6] = $line2[6].",".$filterCallPerPos.",".$filterCaller.",SVID".$svCount;
					}
					$line2[6] =~ s/,{0,1}CV[\d+]//;
					# $line2 = join("\t",@line2);
					push (@output_approx,[@line2]);
					# print $callCounts.",".$o_pos.",".$o_vartype."\n";
					# shift @temp_pos;
				}
				splice (@temp_pos, 0, $callCounts);
				push(@temp_pos,$line);
				$svCount+=1;
			}
			
		} else {
			# initialize/read values from array
			my ($n_chrom,$n_pos,$n_type,$n_len,$vartype,$n_call) = var_line(@line);
			push(@temp_pos,$line);
			
		}

	} else {
		# print SNPS
		push (@output_approx,[@line]);
	}
	
	}
	# print remaining lines
	my $callCountsRem = @temp_pos;
	my %callerTempRem;
	my $callerCountRem = "0";
	
	foreach my $line3 (@temp_pos){
			my @line3 = @{$line3};
			# compare if 2 different callers call this structural variant
			my ($n_chrom,$n_pos,$n_type,$n_len,$n_vartype,$n_call) = var_line(@line3);
			$callerTempRem{$n_call} = $n_type;
	}
	$callerCountRem = keys %callerTempRem;
				
	foreach my $line3 (@temp_pos){
		my $filterCallPerPos = "CSV".$callCountsRem;
		my $filterCaller = "CNSV".$callerCountRem;
		my @line3 = @{$line3};
		if ( $line3[6] eq "." || $line3[6] =~ /CV[0-9]?/){
			$line3[6] = $filterCallPerPos.",".$filterCaller.",SVID".$svCount;
		} else {
			$line3[6] = $line3[6].",".$filterCallPerPos.",".$filterCaller.",SVID".$svCount;
		}
		$line3[6] =~ s/,{0,1}CV[\d+]//;
		$line3 = join("\t",@line3);
		# print "add_test".$line3."\n";
		push (@output_approx,[@line3]);
		# print $callCountsRem.",".@{$line3}[1]."\n";
	}
	# sort array to get correct pos index
	my @output_approx_sorted = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @output_approx;
	
	# foreach my $lp (@output_approx_sorted){
	# 	my @red = splice(@{$lp},0,2);
	# 	print join("\t",@red)."\n";
	# }
	return(@output_approx_sorted);
}

sub snp_bp_match {
	my ($snp_dist) = @_;
	my @temp_pos_up;
	my @temp_pos_down;
	my @output_bps;
	my $tag = "BPS";
	
	foreach my $line (@pos_approx){
		# check if data exists
		my @line = @{$line};
		my $dist = $snpBp_dist;
		my ($n_chrom,$n_pos,$n_type,$n_len,$n_vartype,$n_call) = var_line(@line);
		# print "$n_chrom,$n_pos,$n_type,$n_len,$n_vartype\n";
		# if variant is BP check if SNPs lie within dist window up or/and downstream and tag both with BPS
		# current BP pos is stored in @temp_pos_down, if new pos is BP but still within dist($temp_pos_down[0] exists), then
		# ignore upstream part and proceed downstream with new BP as $temp_pos_down[0]
		if ($n_type eq "BP" && !exists $temp_pos_down[0]){
			# upstream
			my $outdist = "0";
			# while the new variants are within the dist, outdist is 0 and pos lines are written to the temp array
			# when one variant is out of dist, outdist is 1 and the temp array will be tagged and added to the outarray
			# also check if you are at the start and no more lines are available
			while ( $outdist == "0" && exists $output_bps[0]){
				# print "BP outdist 0 and output[0]\n";
				my @line_up = @{pop(@output_bps)};
				# print Dumper(@line_up);
				# my @line_up = @{$line_up};
				my ($up_chrom,$up_pos,$up_type,$up_len,$up_vartype,$up_call) = var_line(@line_up);
				if ($up_chrom eq $n_chrom && ($up_pos + $dist) >= $n_pos){
					# print "within dist $dist\n";
					# if within distance tag and add to top of temp array
					my @filter_up = add_tag($tag,@line_up);
					unshift(@temp_pos_up,[@filter_up]);
				} else {
					# merge arrays
					my @merged = (@output_bps,@temp_pos_up);
					@output_bps = @merged;
					# reset temp_pos_up
					splice(@temp_pos_up,0);
					# add BP call
					push (@output_bps,[@line]);
					# stop iterating
					$outdist = "1";
				}
			}# end loop here
			# downstream: just continue iterating through array
			push (@output_bps,[@line]);
			push( @temp_pos_down,[@line]);
			
		} else {
			# check if position is within dist and tag it, else write to output array
			if ($n_type eq "BP" && exists $temp_pos_down[0]){
				splice (@temp_pos_down, 0);
				push(@temp_pos_down,[@line]);
				push(@output_bps,[@line]);
			} elsif (exists $temp_pos_down[0]) { # if we are tracing BP dist
				my $first_val = "0";# array index
				my ($down_chrom,$down_pos,$down_type,$down_len,$down_vartype,$down_call) = var_temp($first_val,@temp_pos_down);
				if( $down_chrom eq $n_chrom && ($n_pos <= ($down_pos + $dist))){
					# tag and write to array
					my @filter_down = add_tag($tag,@line);
					push (@output_bps,[@filter_down]);
				} else { 
					# if we are out of dist, clear @temp_pos_down array
					splice (@temp_pos_down, 0);
					push (@output_bps,[@line]);
				}	
			} elsif (!exists $temp_pos_down[0]) {# if we have no BP nearby, just print out lines
				push (@output_bps,[@line]);
			}
			
		}# if another BP end
		
	}
	# sort array to get correct pos index
	my @output_bps_sorted = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @output_bps;
	
	return(@output_bps_sorted);
}

sub add_tag {
	my ($tag,@line) = @_;
	
	if ( $line[6] eq "."){
		$line[6] = $tag;
	} else {
		# my $filter =~ s/CV[0-9]?/CSV$tag/;
		# if tag exists, do not add more
		if ($line[6] !~ /$tag/){
			$line[6] = $line[6].",".$tag;
		}
		
	}
	return(@line);
}

# print array
sub print_array {
	my @input = @_;
	foreach my $line (@input){
		# print $line."\n";
		my @line = @{$line};
		print (join("\t",@line)."\n");
	}
}
# print header
sub print_header {
	my @input = @_;
	foreach my $line (@input){
		# print $line."\n";
		print $line."\n";
	}
}

# read in input files
sub read_vcf_files {
	my ($vcf_filename) = @_;
	my @vcf_file = open_vcf($vcf_filename);
	my $counter_line = "0";
	my @coll_lines;
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
    			$vcf_info =~ m/SVLEN=([-,0-9]*)\W/;
    			my $sv_len = $1;
    			
    			# add lines to array/array
    			push (@coll_lines,[@splitline]);
    			
    			
    		} else {
    			# print header
    			push (@header,$line);
    		}
    }
    return (@coll_lines);
}

# open and process file
sub open_vcf {
	my ($ref_name) = @_;
	my @ref_file;
	open (FILE, "<", $ref_name) || die "File not found: ".$ref_name;
	@ref_file = <FILE>;
	close FILE;
	return (@ref_file);
}

