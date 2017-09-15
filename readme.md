## Overview

The main reference is [hep-lat/9503028](https://arxiv.org/pdf/hep-lat/9503028.pdf). It Section 3.2 therein, it says a _clusters algorithm_ is used to undate the spin states configuration, which is, reputedly, more effective.

The details of the _clusters algorithm_ can be found in Ref. _[R. H. Swendsen and J.-S. Wang, Phys. Rev. Lett. 58 (1987) 86]_. Honestly, I have not followed everything in this paper, and I do not know how the detailed balence is assured. But the numerical procedurr is clear. Another similar paper is [hep-lat/9206004](https://arxiv.org/pdf/hep-lat/9206004.pdf). However, I think the first two are enough.

As mentioned by the main reference, the clusters can be labelled by the [Hoshen-Kopelman algorithm](https://en.wikipedia.org/wiki/Hoshen–Kopelman_algorithm). With this tool in hand, together with the steps introduced in the PRL paper and the probabilities of linking given in the main reference, the code can be implemented relatively easily. 

## TODO
- ~~Finish~~ and test the Hoshen-Kopelman method for 4-dim case.
- Implemente the rest part of the code.
 - Maybe start from one field case? 
- Calculate the critical point as a test.
