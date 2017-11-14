#!/usr/bin/perl

if($#ARGV<1){
    print qq{
    Program: matchGFFwithFIMO.pl (split peaks from GFF into those containing, or not containing, motif hits from fimo.gff)
    Written by: Shaun Mahony <mahony\@psu.edu>
    Usage:  perl matchGFFwithFIMO.pl [-d <number> -o <order>] cwpair.gff fimo.gff 

    Options: cwpair.gff       path to the GFF file that contains peak pairs input to FIMO
             fimo.gff       path to the FIMO GFF output file (must be fimo.gff file)
             -d <number>    motif assigned to peak if within this distance, default = 40bp
             -o <peakscore,motifscore>    order to output motif-containing peaks, default = peakscore

    Example:
      perl matchGFFwithFIMO.pl -d 40 -o peakscore peaks.gff fimo.gff
      
    Output: (Files are output to same directory as input GFF file)
      - GFF file containing peaks that have motif matches (recentered on best motif match within "d" bp of peak)
      - GFF file containing peaks without motif match
      - TXT file plotting % of peaks that have a motif match as a function of peak rank (based on original GFF file ordering)
  
  };
}else{

    #read options
    my $d = 40;
    my $order = "peakscore";
    my $rankvmotifwin=200;
    for($i=0; $i<$#ARGV; $i++){
	if($ARGV[$i] eq "-d"){
	    $d = $ARGV[$i+1];
	}
	if($ARGV[$i] eq "-o"){
	    $order = $ARGV[$i+1];
	}
    }
    my $gffFile = $ARGV[$#ARGV-1];
    my $fimoFile = $ARGV[$#ARGV];
    
    my %chrID =();
    my $numChr=0;
    my @chrCount=();
    my @mList=();

    #Read in FIMO motif hit GFF
    unless(open(FIMOFILE, $fimoFile)){
	die "Cannot open file\n";}
    my @fimolines = <FIMOFILE>;
    my @motifhits=();
    my $nummh=0;
    for($x=0; $x<=$#fimolines; $x++){
	chomp($fimolines[$x]);
	@curr = split(/\s+/, $fimolines[$x]);
	if($curr[0] !~ m/^#/){
	    $chr = $curr[0];
	    $chr =~ s/chr//g;
	    $start = $curr[3];
	    $stop = $curr[4];
	    $str = $curr[6];
	    $info = $curr[8];
	    
	    $pval = 0.0;
	    @infobits = split(/\;/, $info);
	    for($ib=0; $ib<=$#infobits; $ib++){
		@pieces = split(/\=/, $infobits[$ib]);
		if($pieces[0] eq "pvalue"){
		    $pval = $pieces[1];
		}
	    }

	    if(!defined($chrID{$chr})){$chrID{$chr}=$numChr; $chrCount[$numChr]=0;$numChr++;}
	    $thisChr = $chrID{$chr};
	    $numC=$chrCount[$thisChr];
	    $mList[$thisChr][$numC][0]=$chr;
	    $mList[$thisChr][$numC][1]=$start;
	    $mList[$thisChr][$numC][2]=$stop;
	    $mList[$thisChr][$numC][3]=$str;     
	    $mList[$thisChr][$numC][4]=$pval;
	    $mList[$thisChr][$numC][5]=$info;
	    $chrCount[$thisChr]++;
	}
    }

    my $baseName = $gffFile;
    $baseName =~ s/.gff$//g;
    my $rvmfile = $baseName."_rankvsmotif.txt";
    unless(open(RANKVMOTIF, ">$rvmfile")){
	die "Cannot open file\n";}
    print RANKVMOTIF "Quantile\tFractionWithMotif\n";


    my @peaksWithMotifs=();
    my $numPeaksWithMotifs=0;
    my @peaksNoMotifs=();
    my $numPeaksNoMotifs=0;
    my $numPeaks=0;
    my $winWithMotif=0;
    my %added=(); #List of positions that are in a list already
    #Read in GFF peaks
    unless(open(GFFFILE, $gffFile)){
	die "Cannot open file\n";}
    my @gfflines = <GFFFILE>;
    for($x=0; $x<=$#gfflines; $x++){
	@curr = split(/\s+/, $gfflines[$x]);
	if($curr[0] !~ m/^#/){
	    $chr = $curr[0];
	    $chr =~ s/chr//g;
	    $start = $curr[3];
	    $stop = $curr[4];
	    $str = $curr[6];
	    $score = $curr[5];
	    $note = $curr[8];
	    $mid = int(($start+$stop)/2);
	    $numPeaks++;

	    $motifStart=-1; $motifStop=-1; $motifStr="+"; $motifPval=1; $motifInfo="";
	    if(defined($chrID{$chr})){
		$thisChr=$chrID{$chr};
		$ccount=$chrCount[$thisChr];
		
		for($j=0; $j<$ccount; $j++){
		    if(abs($mid - $mList[$thisChr][$j][1])<=$d || abs($mid - $mList[$thisChr][$j][2])<=$d){
			if($mList[$thisChr][$j][4]<$motifPval){
			    $motifStart = $mList[$thisChr][$j][1];
			    $motifStop = $mList[$thisChr][$j][2];
			    $motifStr = $mList[$thisChr][$j][3];
			    $motifPval = $mList[$thisChr][$j][4];
			    $motifInfo = $mList[$thisChr][$j][5];
			}
		    }
		}
	    }
	    
	    if($motifStart != -1){ #Has motif
		$id = $chr.":".$motifStart.":".$motifStr;
		if(defined($added{$id})){
		    #If this motif instance is alread added, the only piece that could be changed is the peak score
		    if($order ne "motifscore" && $peaksWithMotifs[$added{$id}][4]<$score){
			$peaksWithMotifs[$added{$id}][4] = $score;
		    }
		    $peaksWithMotifs[$added{$id}][5]="$note\;$motifInfo";
		}else{
		    $peaksWithMotifs[$numPeaksWithMotifs][0]=$chr;
		    $peaksWithMotifs[$numPeaksWithMotifs][1]=$motifStart;
		    $peaksWithMotifs[$numPeaksWithMotifs][2]=$motifStop;
		    $peaksWithMotifs[$numPeaksWithMotifs][3]=$motifStr;
		    $peaksWithMotifs[$numPeaksWithMotifs][4]=$score;
		    if($order eq "motifscore"){
			$peaksWithMotifs[$numPeaksWithMotifs][4]=$motifPval;
		    }
		    $peaksWithMotifs[$numPeaksWithMotifs][5]="$note\;$motifInfo";
		    $added{$id}=$numPeaksWithMotifs;
		    $numPeaksWithMotifs++;
		    $winWithMotif++;
		}
	    }else{ #No motif: assume that all points in this set are unique
		$peaksNoMotifs[$numPeaksNoMotifs][0]=$chr;
		$peaksNoMotifs[$numPeaksNoMotifs][1]=$start;
		$peaksNoMotifs[$numPeaksNoMotifs][2]=$stop;
		$peaksNoMotifs[$numPeaksNoMotifs][3]=$str;
		$peaksNoMotifs[$numPeaksNoMotifs][4]=$score;
		$peaksNoMotifs[$numPeaksNoMotifs][5]=$note;
		$numPeaksNoMotifs++;
	    }

	    if($numPeaks>0 && $numPeaks % $rankvmotifwin==0){
		$perc = $winWithMotif/$rankvmotifwin;
		$lastWin = $numPeaks-$rankvmotifwin+1;
		print RANKVMOTIF "$lastWin-$numPeaks\t$perc\n";
		$winWithMotif=0;
	    }
	}
    }
    close(GFFFILE);
    close(FIMOFILE);
    close(RANKVMOTIF);

    #Sort peaks with motifs
    if($order eq "motifscore"){
	@sortedPeaksWithMotifs = sort { $a->[4] <=> $b->[4] } @peaksWithMotifs;
    }else{
	@sortedPeaksWithMotifs = sort { $b->[4] <=> $a->[4] } @peaksWithMotifs;
    }
    @sortedPeaksNoMotifs = sort { $b->[4] <=> $a->[4] } @peaksNoMotifs;

    #Write the output GFF files
    my $hmfile = $baseName."_withmotif.gff";
    unless(open(HASMOTIF, ">$hmfile")){
	die "Cannot open file\n";}
    my $nmfile = $baseName."_nomotif.gff";
    unless(open(NOMOTIF, ">$nmfile")){
	die "Cannot open file\n";}
    for($a=0; $a<$numPeaksWithMotifs; $a++){
	print HASMOTIF "chr$sortedPeaksWithMotifs[$a][0]\twithmotif\t.\t$sortedPeaksWithMotifs[$a][1]\t$sortedPeaksWithMotifs[$a][2]\t$sortedPeaksWithMotifs[$a][4]\t$sortedPeaksWithMotifs[$a][3]\t.\t$sortedPeaksWithMotifs[$a][5]\n";
    }
    for($a=0; $a<$numPeaksNoMotifs; $a++){
	print NOMOTIF "chr$sortedPeaksNoMotifs[$a][0]\tnomotif\t.\t$sortedPeaksNoMotifs[$a][1]\t$sortedPeaksNoMotifs[$a][2]\t$sortedPeaksNoMotifs[$a][4]\t$sortedPeaksNoMotifs[$a][3]\t.\t$sortedPeaksNoMotifs[$a][5]\n";
    }
    close(HASMOTIF);
    close(NOMOTIF);
    print "Output:\t".$baseName."_rankvsmotif.txt\n";
    print "Output:\t".$baseName."_withmotif.gff\n";
    print "Output:\t".$baseName."_nomotif.gff\n";
}
