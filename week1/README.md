Step 1.1:
( 1,000,000 * 3 ) / 100 = 30,000 reads for 3x coverage

Step 1.4:
About 50,000 positions in our genome have no coverage, which is ~ 5%.
Inspecting the graph shows the poisson distribution fits the data very well. The normal distribution does not fit as well, looking slightly right shifted.

Step 1.5:
106 positions in our genome have no coverage, which is ~ 0.01%.
This simulation is also nicely fit by the poisson distribution. Again, the normal distribution does not fit as well, looking right shifted again. However it looks closer than 3x coverage.

Step 1.6:
8 positions in our genome have no coverage, which is 0.0008%.
This simulation is also nicely fit by the poisson distribution. The normal distribution also looks like it fits well here.

Step 2.4:
dot -Tpng ex2_edges.dot -o ex2_digraph.png

Step 2.5:
ATTCATTGATTGATTCTTATTT

Step 2.6:
The kmers you use to make a directed graph should not repeat. That way there are no ambigous walks through the graph. Longer reads will also help.



