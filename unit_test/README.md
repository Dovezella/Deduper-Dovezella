## Deduper Unit Test Read ME

For this assignment I wanted to try to cover as many cases as I could think of before hand. So I created a fake same file called unit_tes_input.sam located in the unit_test folder of Deduper-Dovezella.

The headers were taken from an example and shouldn't change.
Starting on the first read line:
28: unique write out 
29: Unique - write out
30: Unique - write out
31: write to erunk (N in umi)
32: unique write out
33: unique write out
34: write to erunk (unknown umi)
35: unique write out
36: unique write out
37: unique write out
38: unique write out
39: write to dupe (dupe of 37)
40: unique write out (similar to 39 but diff umi)
41: unique write out
42: unique write out
43: unique write out
44: unique write out (similar to 43 but reverse)
45: write to dupe (dupe of 44) 
46: unique write out (forward, left soft clip) 
47: unique write out (reverse, right soft clip)
48: write to dupe (dupe of 47)
49: unique write out (Reverse, left soft clip)
50: unique write out (insertion)
51: dupe of line 50 after adjust_pos (deletion), write dupe