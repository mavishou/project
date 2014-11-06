#!/usr/bin_crush/perl


$rnafilename = shift( @ARGV );
unless ( open(RNAFILE, $rnafilename) ) {
    print "Cannot open file \"$rnafilename\"\n\n";
    goto h;
}
while(<RNAFILE>){
	$line = $_;
	if( ! ( $line =~ /^>/ ) ) {
		$line =~ s/T/U/g;
		$line =~ s/t/u/g;
	}
	print $line;
}

close RNAFILE;
