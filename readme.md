The main reference is [hep-lat/9503028](https://arxiv.org/pdf/hep-lat/9503028.pdf). In Section 3.2 therein, it says a _clusters algorithm_ to undate the spin states, which is, reputedly, more effective.

The details of the _clusters algorithm_ can be found in Ref. _[R. H. Swendsen and J.-S. Wang, Phys. Rev. Lett. 58 (1987) 86]_. Another similar paper is [hep-lat/9206004](https://arxiv.org/pdf/hep-lat/9206004.pdf). However, I think the first two are enough.

As mentioned by the main reference, the clusters can be labelled by the [Hoshen-Kopelman algorithm](https://en.wikipedia.org/wiki/Hoshen–Kopelman_algorithm). With this tool in hand, together with the steps introduced in the PRL paper and the probabilities of the linking given in the main reference, the code can be implemented relatively easily.

Please note the definition of the effective $\kappa$,
$$\kappa_{eff}(x,\mu)=\kappa_{\phi}-\frac{1}{2}g(\rho_x+\rho_{x+\mu}).$$