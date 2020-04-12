# Network-TDA

This repository contains the code for the algorithms to construct and update **community trees** to summarize the topological structures in a network. The concept of  community tree, a tree structure established based on clique communities from the clique percolation method, is proposed in the following paper
```
@article{chen2017note,
  title={A note on community trees in networks},
  author={Chen, Ruqian and Chen, Yen-Chi and Guo, Wei and Banerjee, Ashis G},
  journal={arXiv preprint arXiv:1710.03924},
  year={2017}
}
```

With the information revealed by community trees and the corresponding persistence diagrams, our proposed approach is able to detect clique communities and keep track of the major structural changes during their evolution given a stability threshold (TSN bound).

<p align="center">
    <img src="https://github.com/w-guo/Network-TDA/blob/master/images/sms_evolution.png" width="900"> <br />
    <em> Network evolution for the short-msg dataset</em>
</p>
