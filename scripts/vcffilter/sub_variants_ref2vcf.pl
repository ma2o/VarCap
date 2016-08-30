# file contains output subroutines for collect_variants.pl

# subroutine, that annotates ref base for position for vcf format
sub pos2base {
	my ( $pos, $ref_format ) = @_;
	my $start = $pos - 1;
	if ($start < "0"){
		$start = "0";
	}
	my $length = "1";
	my @ref_fasta = split(/\n/,$ref_format);
	my $ref_prebase = substr($ref_fasta[1], $start, $length);
	return ($ref_prebase);
}

# add 1(true)
1;