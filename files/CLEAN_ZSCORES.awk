{
a1=toupper($2);
a2=toupper($3);
ref1 = toupper($7);
ref2 = toupper($8);

ref = ref1;
if ( ref == "A" ) flip = "T";
else if ( ref == "T" ) flip = "A";
else if ( ref == "G" ) flip = "C";
else if ( ref == "C" ) flip = "G";
flip1 = flip;

ref = ref2;
if ( ref == "A" ) flip = "T";
else if ( ref == "T" ) flip = "A";
else if ( ref == "G" ) flip = "C";
else if ( ref == "C" ) flip = "G";
flip2 = flip;

z = 0;

if( !((a1=="A" && a2=="T") || (a1=="T" && a2=="A") || (a1=="C" && a2=="G") || (a1=="G" && a2=="C")) ) {

if( a1 == ref1 && a2 == ref2 ) z = $5;
else if ( a1 == ref2 && a2 == ref1 ) z = -1*$5;
else if ( a1 == flip1 && a2 == flip2 ) z = $5;
else if ( a1 == flip2 && a2 == flip1 ) z = -1*$5;
}

if ( z != 0 ) print $1,a1,a2,z,$6
}
