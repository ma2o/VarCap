# subs for variant extraction for collect_variants.pl

# extract data from result output files and call sub to write to hashtable
sub extract_samtools {
    my ($infile, $repeat, $caller, $min_reads) = @_;
    # my $type;
    # open .vcf file from samtools
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line !~ /^#/ ){
    		# print $line;
    		my @splitline = split(/\t/, $line);
    		chomp(@splitline);
    		# print "length: ".@splitline." : ".$splitline[1]."\n";
    		my $type = "UNDEF";
    		my $length = "0";
    		if ($splitline[7] =~ /^INDEL/ ){
    		  $type = "IND";
    		  $length = "1";
    		} else {
    		  $type = "SNP";
    		  $length = "0";
    		}
    		$splitline[7] =~ /DP4=(\d*),(\d*),(\d*),(\d*)\;/;
    		my $reads_ref = ($1+$2);
    		my $reads_alt = ($3+$4);
    		my $ref_base = $splitline[3];
    		my $var_base = $splitline[4];
    		# $var_base =~ s/\,/\;/;
    		my $chrom = $splitline[0];
    		# feed data to hashtable
    		if ($reads_alt > $min_reads){
    			add_data($splitline[1], $caller, $repeat, $type, $length, $ref_base, $var_base, $reads_ref, $reads_alt, $chrom);
    			# print "$splitline[1], $caller, $repeat, $type, $length, $ref_base, $var_base, $reads_ref, $reads_alt, $chrom\n";
    		}
        }
    }
}

sub extract_gatk {
	my ($infile, $repeat, $caller, $min_reads) = @_;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line =~ /^\w/ ){
    		# print $line;
    		my @splitline = split(/\t/, $line);
    		my $type = "non";
    		my $length;
    		if ($splitline[7] =~ /RU=/){
    			$length = "1";
    			$type = "IND";
    		} else {
    			$length = "0";
    			$type = "SNP";
    		}
    		my @genome_data = split(/:/, $splitline[9]);
    		my @read_stats = split(/,/, $genome_data[1]);
    		my $reads_ref = $read_stats[0];
    		my $reads_alt = $read_stats[1];
    		my $ref_base = $splitline[3];
    		my $var_base = $splitline[4];
    		my $chrom = $splitline[0];
    		# feed data to hashtable
    		if ($reads_alt > $min_reads){
    			# print "reads_alt: ".$reads_alt."\n";
    			add_data($splitline[1], $caller, $repeat, $type, $length, $ref_base, $var_base, $reads_ref, $reads_alt, $chrom);
    		}
        }
    }
}

sub extract_varscan {
    my ($infile, $repeat, $caller, $min_reads) = @_;
    my $type;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line =~ /^\w/ && $line !~ /^Chrom/){
    		#print $line;
    		my @splitline = split(/\t/, $line);
    		chomp(@splitline);
    		#print "length: ".@splitline." : ".$splitline[1]."\n";
    		if ($infile =~ /snp$/ ){
    		  $type = "SNP";
    		} else {
    		  $type = "IND";
    		}
    		my $reads_ref = $splitline[4];
    		my $reads_alt = $splitline[5];
    		my $length;
    		my ($ref_base, $var_base);
    		if ($type eq "SNP"){
    			$ref_base = $splitline[2];
    		    $var_base = $splitline[18];
    		    $length = "0";
    		} else {
    			$splitline[18] =~ /^([+-]?)(\w+)/;
    			if ($1 eq "-"){
    				$ref_base = $splitline[2].$2;
    		        $var_base = $splitline[2];
    		        $length = $1."1";
    			} else {
    				$ref_base = $splitline[2];
    		        $var_base = $splitline[2].$2;
    		        $length = "1";
    			}
    		}
    		my $chrom = $splitline[0];
    		# feed data to hashtable
    		if ($reads_alt > $min_reads){
    			add_data($splitline[1], $caller, $repeat, $type, $length, $ref_base, $var_base, $reads_ref, $reads_alt, $chrom);
    		}
        }
    }
}

sub extract_lofreq {
    my ($infile, $repeat, $caller, $min_reads) = @_;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line =~ /^\w/){
    		#print $line;
    		my @splitline = split(" ", $line);
    		chomp(@splitline);
    		my $chrom = $splitline[0];
    		my $pos = $splitline[1];
    		my $type = "SNP";
    		my $length = "0";
    		# get bases
    		my ($ref_base, $var_base);
    		my $bases = $splitline[2];
    		$bases =~ m/^(\w{1})/;
    		$ref_base = $1;
    		$bases =~ m/>(\w{1})/;
    		$var_base = $1;
    		$bases =~ m/>(.*)$/;
    		my $all_var_bases = $1;
    		# get reads
    		my $reads_string = $splitline[4];
    		my %bases_reads;
    		$reads_string =~ m/basecount-A.(\w{1,5})/;
    		$bases_reads{"A"} = $1;
    		$reads_string =~ m/basecount-C.(\w{1,5})/;
    		$bases_reads{"C"} = $1;
    		$reads_string =~ m/basecount-G.(\w{1,5})/;
    		$bases_reads{"G"} = $1;
    		$reads_string =~ m/basecount-T.(\w{1,5})/;
    		$bases_reads{"T"} = $1;
    		my $reads_ref = $bases_reads{$ref_base};
    		my $reads_alt = $bases_reads{$var_base};
    		
    		# feed data to hashtable
    		if ($reads_alt > $min_reads){
    			add_data($pos, $caller, $repeat, $type, $length, $ref_base, $all_var_bases, $reads_ref, $reads_alt, $chrom);
    		}
        }
    }
}

sub extract_lofreq2 {
    my ($infile, $repeat, $caller, $min_reads) = @_;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line =~ /^\w/){
    		#print $line;
    		my @splitline = split("\t", $line);
    		chomp(@splitline);
    		my $chrom = $splitline[0];
    		my $pos = $splitline[1];
    		my $type = "SNP";
    		my $length = "0";
    		my $ref_base = $splitline[3];
    		my $var_base = $splitline[4];
    		# get reads out of info field
    		my ($reads_ref, $reads_alt);
    		my @info = split(";",$splitline[7]);
    		foreach my $info_tag (@info){
    			if($info_tag =~ m/^DP4/){
    				$info_tag =~ s/DP4=//;
    				my @reads = split(',',$info_tag);
    				$reads_ref = $reads[0] + $reads[1];
    				$reads_alt = $reads[2] + $reads[3];
    			}
    		}
    		# feed data to hashtable
    		if ($reads_alt > $min_reads){
    			add_data($pos, $caller, $repeat, $type, $length, $ref_base, $var_base, $reads_ref, $reads_alt, $chrom);
    		}
        }
    }
}

sub extract_freebayes {
    my ($infile, $repeat, $caller, $min_reads) = @_;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line =~ /^\w/){
    		#print $line;
    		my @splitline = split("\t", $line);
    		chomp(@splitline);
    		my $chrom = $splitline[0];
    		my $pos = $splitline[1];
    		# my $type = "SNP";
    		# my $length = "0";
    		my $ref_base = $splitline[3];
    		my $var_base = $splitline[4];
    		# get reads out of info field
    		my ($reads_ref, $reads_alt, $type, $length);
    		my @info = split(";",$splitline[7]);
    		foreach my $info_tag (@info){
				# ref reads
    			if($info_tag =~ m/^RO=/){
    				$info_tag =~ s/RO=//;
    				$reads_ref = $info_tag;
    			}
    			# alt reads
    			if($info_tag =~ m/^AO=/){
    				$info_tag =~ s/AO=//;
    				$reads_alt = $info_tag;
    			}
    			# type
    			if($info_tag =~ m/^TYPE/){
    				$info_tag =~ s/TYPE=//;
    				$type = uc $info_tag;
				# resolve COMPLEX type into SNP or INS
				if($type eq "COMPLEX" && ( length $ref_base == length $var_base ) ){
					$type = "SNP";
    				} elsif ($type eq "COMPLEX" && ( length $ref_base < length $var_base ) ){
					$type = "INS";
				} elsif ($type eq "COMPLEX" && ( length $ref_base > length $var_base ) ){
                                        $type = "DEL";
                                } elsif ( $type eq "MNP" ){
					$type = "SNP";
				}
			}
    			# length
    			if($info_tag =~ m/^LEN=/){
    				$info_tag =~ s/LEN=//;
    				if ( $type == "SNP" ){
						$info_tag = "0";
					}
    				$length = $info_tag;
    			}
    			
    		}
    		# feed data to hashtable
    		if ($reads_alt > $min_reads){
    			add_data($pos, $caller, $repeat, $type, $length, $ref_base, $var_base, $reads_ref, $reads_alt, $chrom);
    		}
        }
    }
}

sub extract_breakdancer {
    my ($infile, $repeat, $caller, $ref_format, $min_reads) = @_;
    my $type;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line =~ /^\w/ ){
    		#print $line;
    		my @splitline = split(/\t/, $line);
    		#print "length: ".@splitline." : ".$splitline[1]."\n";
    		$type = $splitline[6];
    		my $length = $splitline[7];
    		my $pos = $splitline[1];
    		my $ref = pos2base($pos, $ref_format);
    		my $alt = $splitline[4];
    		my $reads_ref = "nd";
    		my $reads_alt = $splitline[9];
    		my $chrom = $splitline[0];
    		#feed data to hashtable
    		if ($reads_alt > $min_reads){
    			add_data($pos, $caller, $repeat, $type, $length, $ref, $alt, $reads_ref, $reads_alt, $chrom);
    		}
        }
    }
}

sub extract_delly {
    my ($infile, $repeat, $caller, $ref_format, $min_reads) = @_;
    my $type;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line =~ /<$/ ){
    		#print $line;
    		my @splitline = split(/\t/, $line);
    		#print "length: ".@splitline." : ".$splitline[1]."\n";
    		$splitline[6] =~ /^>(...)/;
    		$type = uc ($1);#uppercase
    		my $reads_alt = $splitline[4];
    		# add filter, to add only if supporting read_nr > 3
    		if ($reads_alt > $min_reads){
    		  my $length = $splitline[3];
    		  my $pos = $splitline[1];
    		  my $ref = pos2base($pos, $ref_format);
    		  my $alt = $splitline[2];
    		  my $reads_ref = "nd";
    		  my $chrom = $splitline[0];
    		  #feed data to hashtable
    		  add_data($pos, $caller, $repeat, $type, $length, $ref, $alt, $reads_ref, $reads_alt, $chrom);
    	    }
        }
    }
}

sub extract_delly_072 {
    my ($infile, $repeat, $caller, $ref_format, $min_reads) = @_;
    my $type;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line =~ /^\w/ ){
    		#print $line;
    		my @splitline = split(/\t/, $line);
    		my @gtinfo = split(/:/, $splitline[9]);
    		my $reads_ref = $gtinfo[8];
    		my $reads_alt = $gtinfo[9];
    		# add filter, to add only if supporting read_nr > 3
    		
    		if ($reads_alt > $min_reads){
    		  $splitline[7] =~ /END=(\d+)/;
    		  my $endpos = $1;
    		  my $pos = $splitline[1];
    		  my $length = ($endpos - $pos);
    		  my $ref = pos2base($pos, $ref_format);
    		  my $alt = $splitline[4];
    		  my $chrom = $splitline[0];
    		  $splitline[7] =~ /SVTYPE=(\w+)/;
    		  my $type = $1;
    		  #feed data to hashtable
    		  add_data($pos, $caller, $repeat, $type, $length, $ref, $alt, $reads_ref, $reads_alt, $chrom);
    	    }
        }
    }
}

sub extract_pindel {
    my ($infile, $repeat, $caller, $ref_format, $min_reads) = @_;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	# DEL/INS are vcf files
    	if ( $infile =~ /vcf$/ && $line !~ /^#/ ){
    		my @splitline = split(/\t/, $line);
    		#print "length: ".@splitline." : ".$splitline[1]."\n";
    		my $pos = $splitline[1];
    		$splitline[7] =~ /SVTYPE=(\w+)/;
    		my $type = $1;
    		$splitline[7] =~ /SVLEN=(-{0,1}\d+)/;
    		my $length = $1;
    		my $reads_ref = "nd";
    		my $reads_alt = "nd";
    		if ($splitline[9] =~ /:[0-9]*\,/ ){
    			$splitline[9] =~ /:(\d+)\,(\d+)/;
    			# $reads_ref = $1;
    			$reads_alt = $2;
    		} else {
    			$splitline[9] =~ /:(\d+)/;
    			$reads_alt = $1;
    		}
    		my $chrom = $splitline[0];
    		#feed data to hashtable
    		if ($reads_alt > $min_reads){
                # print "$pos, $caller, $repeat, $type, $length, $splitline[3], $splitline[4], $reads_ref, $reads_alt, $chrom";
    			add_data($pos, $caller, $repeat, $type, $length, $splitline[3], $splitline[4], $reads_ref, $reads_alt, $chrom);
    		}
        }
        # BP file lines start with ChrID
        if ( $infile =~ /BP$/ && $line =~ /^ChrID/ ){
        	my @splitline = split(/\s+/, $line);
        	my $pos = $splitline[2];
        	my $caller_mod = $caller."_B";
        	my $length = "non";
        	my $type = "BP";
        	my $ref = pos2base($pos, $ref_format);
        	#print "ref_BP:".$ref."\n";
        	my $alt = "nd";
        	my $reads_ref = "nd";
        	my $reads_alt = $splitline[6];
        	my $chrom = $splitline[1];
        	if ($reads_alt > $min_reads){
        		add_data($pos, $caller_mod, $repeat, $type, $length, $ref, $alt, $reads_ref, $reads_alt, $chrom);
        	}
        }
        # LI file lines start with digit
        if ( ($infile =~ /LI$/)  && $line =~ /^\d/ ){
        	my @splitline = split(/\s+/, $line);
        	my $pos = $splitline[4];
        	my $caller_mod = $caller."_LI";
        	my $length = "non";
        	my $type = "LI";
        	my $ref = pos2base($pos, $ref_format);
        	#print "ref_BP:".$ref."\n";
        	my $alt = "nd";
        	my $reads_ref = "nd";
        	my $reads_alt = $splitline[6] + $splitline[9];
        	my $chrom = $splitline[3];
        	if ($reads_alt > $min_reads){
        		add_data($pos, $caller_mod, $repeat, $type, $length, $ref, $alt, $reads_ref, $reads_alt, $chrom);
        	}
        }
    }
}

sub extract_cortex {
    my ($infile, $repeat, $caller, $min_reads) = @_;
    my $type;
    open(my $handle, "<", $infile) || die "cannot open $infile: $!\n";
    while( my $line = <$handle>){
    	if ( $line =~ /^\w/ ){
    		#print $line;
    		my @splitline = split(/\t/, $line);
    		#print "length: ".@splitline." : ".$splitline[7]."\n";
    		$splitline[7] =~ /SVLEN=(-{0,1}\d+)/;
    		my $length = $1;
    		$splitline[7] =~ /SVTYPE=(\w+)/;
    		$type = $1;
    		my @genome_data = split(/:/, $splitline[9]);
    		my @read_stats = split(/,/, $genome_data[1]);
    		my $reads_ref = $read_stats[0];
    		my $reads_alt = $read_stats[1];
    		my $chrom = $splitline[0];
    		# if ($chrom eq "NC_005861.1"){
    		#	$chrom = "gi|46445634|ref|NC_005861.1|";
    		# }
    		#feed data to hashtable
    		if ($reads_alt > $min_reads){
    			add_data($splitline[1], $caller, $repeat, $type, $length, $splitline[3], $splitline[4], $reads_ref, $reads_alt, $chrom);
    		}
        }
    }
}

# add the following line containing 1 to return true
1;
