# Network-TDA

This repository contains the code for the algorithms to construct and update **community trees** to summarize the topological structures in a network. Please cite the following paper if you make use of the code.
```
@article{guo2022efficient,
  title={Efficient Community Detection in Large-Scale Dynamic Networks Using Topological Data Analysis},
  author={Guo, Wei and Chen, Ruqian and Chen, Yen-Chi and Banerjee, Ashis G},
  journal={arXiv preprint arXiv:2204.03191},
  year={2022}
}
```

The concept of  community tree, a tree structure established based on clique communities from the clique percolation method, is proposed by
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
    <img src="https://github.com/w-guo/wguo/blob/master/content/publication/network-tda/sms_evolution.png" width="840"> <br />
    <em> Network evolution for the short message correspondence dataset</em>
</p>
