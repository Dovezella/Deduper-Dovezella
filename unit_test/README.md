## Deduper Unit Test Read ME

For this assignment I wanted to try to cover as many cases as I could think of before hand. So I created a fake same file called unit_tes_input.sam located in the unit_test folder of Deduper-Dovezella.

The headers were taken from an example and shouldn't change.
Starting on the first read line:
25: This represents a unique read and should be written to output sam file. 
26: Unique - write out
27: Unique - write out
28: Despite the same start position and other parameters as line 27, it is reversed and the adjusted position would make it unique. write it out with fixed position and rev_comped umi
29: unmapped don't write out
30: secondary alignment, don't write out
31: secondary, reverse strand, don't write out
32: unique, but reverse so fix start position and rev_comp the umi, write out
33: unique, reverse strand, but soft clipped. fix position and write out
34: unique, reverse strand, soft clipped at the 3' end after adjustment so no need to worry about the soft clip. write out
35: unique, forward strand, soft clipped at 5' end. fix position and write out
36: reverse strand that matches 33 but w/o soft clipping. unique. fix position and write out with rev_comped umi
37: matches 26, but doesn't align as well because of insert and different umi. don't write out
38: duplicate of 26, don't write out
39: duplicate of 26 with different umi, don't write out
40: duplicate of 28, don't write out
41: unique, reversed, with indent; adjust position considering indent and write out with rev_comped umi
42: unique, reversed, with deletion; adjust position considering deletion and write out with rev_comped umi
43: unique, reversed, with large skip from splice; adjust position considering skip and write out with rev_comped umi
44: unique, has skipped splice, write out
45: duplicate of 44, different umi, don't write out
46: unique, but reverse so fix start position and rev_comp the umi, write out
47: unknown umi, don't write out
48: error in umi, don't write out
49: unmapped, reverse, don't write out