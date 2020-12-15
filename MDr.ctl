PDB_ID	1yet_unboundREDUCED	# Protein Data Bank ID for MD run
Number_Chains	1	# Number of chains on structure
lengthA	213	 #end of chain designated
LIGAND_ID	1yet_ligandREDUCED	# Protein Data Bank ID for MD run
WATER_ID	1yet_waterREDUCED	# Protein Data Bank ID for MD run
Force_Field	leaprc.protein.ff14SB	# AMBER force field to use in MD runs
LIGAND_Field	leaprc.gaff2	# AMBER force field to use in MD runs
Box_Size	40	# water box size (buffer dist in angtroms)
Number_Runs	5	# number of repeated samples of MD runs
Heating_Time	100000	# length of heating run (fs)
Equilibration_Time	100000	# length of equilibration run (fs)
Production_Time	100000	# length of production run (fs)
Solvation_Method	explicit	# method of solvation (implicit or explicit)
Salt_Conc	0	# salt concentration (implicit only, PME=O)
