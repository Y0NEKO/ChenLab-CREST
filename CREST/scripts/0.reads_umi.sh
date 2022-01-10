#!/usr/bin/

BEGIN {
	FS="\t";
	OFS="\t";
	sep=":"
	xcregex = "CB:Z:[ACTGN0-9]+"
	xmregex = "UB:Z:[ACTGN0-9]+"
	print "Read.ID", "Read.Seq", "Cell.BC", "UMI"
}
{
	match($0, xcregex);
	xc = substr($0, RSTART, RLENGTH);
	xc = split(xc, CBC, sep);
	xc = CBC[3];
	if(!xc){
		match($0, "CR:Z:[ACTGN0-9]+");
		xc = substr($0, RSTART, RLENGTH);
		xc = split(xc, CBC, sep)
		xc = CBC[3]
	}

	match($0, xmregex);
	xm = substr($0, RSTART, RLENGTH);
	xm = split(xm, UMI, sep);
	xm = UMI[3]
	#match($10,/TTC(([AT][CG])+)GGATCC/, celltag);
	#if( xc && xm && celltag[1] ){
		print $1, $10, xc, xm;
	#}
}
