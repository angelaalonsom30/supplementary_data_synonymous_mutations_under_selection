#!/usr/bin/perl
use strict; 
use warnings;


## GOAL: To count the number of substitutions of each type in each stochastic mapped tree
## 
## The file format we're working with: (taken from the phylobayes manual)
######
## The substitution mappings are written into one single (large) file,
## using an extension of the newick format. Here is an example of a substitution mapping for
## a subtree of 3 taxa, Tax1, Tax2 and Tax3:
## ((Tax1_I:0.4:I:0.6:M:0.5:L,Tax2_L:0.7:L)L:0.6:L,Tax3_L:1.8:L)L;
## The substitution mapping is such that the ancestor is in state L, and 2 substitutions have
## occurred along the branch leading to Tax1, successively, a transition from L to M, and a
## transition from M to I:
## 
## |------------------ (L) Tax3
## |(L)
## |      |------- (L) Tax2
## |------|(L)
##        |-----/------/---- (I) Tax1
##             L->M   M->I
######
##
## **Note that the order of substitutions in the newick format is backwards in time.

## Run the script in a command terminal using the following command:
## perl tree_parser.pl infile.map > outfile.csv
##

my ($file) = @ARGV;
open my $in, $file or die "Could not open $file: $!";

#print "Tip_count, ";
print "A_T, A_C, A_G, T_A, T_C, T_G, C_A, C_T, C_G, G_A, G_T, G_C \n";

while( my $line = <$in>)  {  
       chomp $line;

	   next unless length($line); 	   
	   
	   my @spl = split(m[(?:[',',\),\(])+], $line);
	   
	   my $tip_count = 0;
	   my $orig_state;
	   
	   my $A_T = 0;
	   my $A_C = 0;
	   my $A_G = 0;
	   
	   my $T_A = 0;
	   my $T_C = 0;
	   my $T_G = 0;
	   
	   my $C_A = 0;
	   my $C_T = 0;
	   my $C_G = 0;
	   
	   my $G_A = 0;
	   my $G_T = 0;
	   my $G_C = 0;
	   
	   
	   foreach my $tips (@spl) {
		   
		   if(!index($tips, "_") == 0 ) {

			   $tip_count++;
		   }
		   $tips =~ s/\:\d+\.\d+\:/,/g;

		   my @states = split( "," , $tips);
		   shift @states for 1;
		   
		   my $nt0 = $states[0];
		   shift @states for 1;		   
		   foreach my $mut (@states){
			   
			   if( $mut eq "A" && $nt0 eq "T" ) { $A_T++; }
			   if( $mut eq "A" && $nt0 eq "C" ) { $A_C++; }
			   if( $mut eq "A" && $nt0 eq "G" ) { $A_G++; }
			   
			   if( $mut eq "T" && $nt0 eq "A" ) { $T_A++; }
			   if( $mut eq "T" && $nt0 eq "C" ) { $T_C++; }
			   if( $mut eq "T" && $nt0 eq "G" ) { $T_G++; }
			   
			   if( $mut eq "C" && $nt0 eq "A" ) { $C_A++; }
			   if( $mut eq "C" && $nt0 eq "T" ) { $C_T++; }
			   if( $mut eq "C" && $nt0 eq "G" ) { $C_G++; }
			   
			   if( $mut eq "G" && $nt0 eq "A" ) { $G_A++; }
			   if( $mut eq "G" && $nt0 eq "T" ) { $G_T++; }
			   if( $mut eq "G" && $nt0 eq "C" ) { $G_C++; }
			   
			   $nt0 = $mut;

		   }

	   }
#	   print "$tip_count, ";
	   print "$A_T, $A_C, $A_G, $T_A, $T_C, $T_G, $C_A, $C_T, $C_G, $G_A, $G_T, $G_C \n";	   
	
}

close $in;