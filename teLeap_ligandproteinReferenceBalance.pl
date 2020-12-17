#!/usr/bin/perl -w
#use warnings;
#use strict;
use File::Copy;
use List::Util qw(shuffle);

# specify the path to working directory for teLeap here
open(IN, "<"."paths.ctl") or die "could not find paths.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $path = @INrow[1];
	 if ($header eq "teleap_path"){$teleap_path = $path;}
}
close IN;
print "path to teLeap .exe\t"."$teleap_path\n";

###########################################################

print "control file inputs\n\n";

open(IN, "<"."MDr.ctl") or die "could not find MD.ctl control file\n";
@IN = <IN>;
for (my $c = 0; $c <= scalar @IN; $c++){
    $INrow = $IN[$c];
    @INrow = split (/\s+/, $INrow);
    $header = $INrow[0];
    $value = $INrow[1];
    print "$header\t"."$value\n";
    if ($header eq "PDB_ID") { $PDB_IDref = $value;}
    if ($header eq "Force_Field") { $Force_Field = $value;}
    if ($header eq "LIGAND_Field") { $LIGAND_Field = $value;}
    if ($header eq "Box_Size") { $Box_Size = $value;}
    if ($header eq "Number_Runs") { $Number_Runs = $value;}
    if ($header eq "Heating_Time") { $Heating_Time = $value;}
    if ($header eq "Equilibration_Time") { $Equilibration_Time = $value;}
	if ($header eq "Production_Time") { $Production_Time = $value;}
    if ($header eq "Solvation_Method") { $Solvation_Method = $value;}

}

open(IN, "<"."MDq.ctl") or die "could not find MD.ctl control file\n";
@IN = <IN>;
for (my $c = 0; $c <= scalar @IN; $c++){
    $INrow = $IN[$c];
    @INrow = split (/\s+/, $INrow);
    $header = $INrow[0];
    $value = $INrow[1];
    print "$header\t"."$value\n";
    if ($header eq "PDB_ID") { $PDB_IDquery = $value;}
    
}
close IN;


#my $protein_label = $ARGV[0];
my $protein_labelR = $PDB_IDref;
my $protein_labelQ = $PDB_IDquery;
my $method = $Solvation_Method; # "explicit" or "implicit"
my $prmtop;
my $igb;
my $ntb;
my $cut;

if ($method eq "explicit") {
	$prmtop = "wat"; # "vac" or "wat"
	$igb = 0;
	$ntb = 1;
	$cut = 8.5;
}

if ($method eq "implicit") {
	$prmtop = "vac"; # "vac" or "wat"
	$igb = 1;
	$ntb = 0;
	$cut = 999;
}

my $num_runs = $Number_Runs; # Number of repeated production runs
my $len_prod = $Production_Time; # Length of each production run in fs (nstlim value)
my $len_eq = $Equilibration_Time; # Length of equilibration run in fs
my $len_heat = $Heating_Time; # Length of heat run in fs
my $proteinfield = $Force_Field; # specify AMBER protein forcefield
my $ligandfield = $LIGAND_Field; # specify AMBER ligand forcefield

=pod

if (-e "$protein_label.pdb") { print "$protein_label.pdb found\n"; }
#print "Reducing $protein_label\n";
#system("pdb4amber -i $protein_label.pdb -o reduced_$protein_label.pdb --reduce --dry 2> $protein_label"."reduce.log");# --reduce --dry");

# PDBs further reduced manually by deleting heteroatoms
=cut

####################################################################
# find bound protein (complex) water box size
####################################################################
open(IN, "<"."wat_$protein_labelQ.inpcrd") or die "could not open input coordinate file\n";
@IN = <IN>;
for (my $a = 0; $a < scalar @IN; $a++){  # catches values from last line of .inpcrd file
    $INrow = $IN[$a];
    @INrow = split (/\s+/, $INrow);
    $Xval = $INrow[1];
    $Yval = $INrow[2];
    $Zval = $INrow[3];
    if ($Xval != '') {$Xsize = $Xval;};
    if ($Yval != '') {$Ysize = $Yval;};
    if ($Zval != '') {$Zsize = $Zval;};
    print "$Xval\t"."$Yval\t"."$Zval\n";
}
print "protein-ligand complex water box size\n";
print "$Xsize\t"."$Ysize\t"."$Zsize\n";
sleep(1);
close IN;

####################################################################
# Protein: Prepare the input file for tleap 
#################################################################

open(LEAP_PROTEIN, ">"."$protein_labelR"."adjust.bat") or die "could not open LEAP file\n";
	print LEAP_PROTEIN "source "."$teleap_path"."$proteinfield\n";
     print LEAP_PROTEIN "source "."$teleap_path"."$ligandfield\n";
	print LEAP_PROTEIN "source "."$teleap_path"."leaprc.water.tip3p\n";
     print LEAP_PROTEIN "protein$protein_labelR = loadpdb $protein_labelR"."adjust.pdb\n";
     print LEAP_PROTEIN "addions protein$protein_labelR Na+ 0\n"; # to charge or neutralize explicit solvent
     print LEAP_PROTEIN "addions protein$protein_labelR Cl- 0\n"; # to charge or neutralize explicit solvent
	print LEAP_PROTEIN "set protein$protein_labelR box {$Xsize $Ysize $Zsize}\n";
     #print LEAP_PROTEIN "setBox protein$protein_label vdw {$Box_Size $Box_Size $Box_Size}\n";
     print LEAP_PROTEIN "saveamberparm protein$protein_labelR wat"."_$protein_labelR.prmtop wat"."_$protein_labelR.inpcrd\n";
     print LEAP_PROTEIN "quit\n";
close LEAP_PROTEIN;

print "  preparing input file for teLeap\n\n";
sleep(1);

######################################################################
print "  edit teLeap setup if desired (e.g. change water model, box size/shape, or charging)\n";
print "  default is simple rigid 3 point model, charge neutralized with Na+\n";
print "  close .bat when done\n\n";
sleep(2);

system "gedit $protein_labelR"."adjust.bat\n";

######################################################################################
# Run sequence through tleap: prepare topology (prmtop) and coordinate (inpcrd) files
######################################################################################
open(TLEAP_PROTEIN, '|-', "tleap -f $protein_labelR"."adjust.bat");
	print<TLEAP_PROTEIN>;
close TLEAP_PROTEIN;


sleep(1);

######################################################################

print "teLeap on Reference structure is complete\n\n";

exit;
