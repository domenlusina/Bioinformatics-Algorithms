                     GenCompress
-----------------------------------------------------------------------------------------------

 The attached is the linux version of GenCompress program.
 
 GenCompress can effectively compress DNA sequences.
 If a file is composed of only four lower-case characters
 of {a, c, g, t}, it will be considered as a DNA sequence file by GenCompress program. 
 
 Usage :  GenCompress original.dat  [-c  reference.dat]  [-o  k  b]
 
 Result : A compressed file called original.GEN will be obtained.
 
 a> If one want to compress a DNA sequence file, run "GenCompress original.dat";
 
 b> If one want to compress a DNA sequence file conditional to another DNA sequence file, 
   i.e., Compress(original_file|reference_file), run "GenCompress original.dat -c reference.dat";
 
 c> [-o k b] is an optional parameter described in our paper. Its default value is "-o 12 3".

-------------------------------------------------------------------------------------------------


Xin Chen  chxin@cs.ucsb.edu