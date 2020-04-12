/*
 * EulerTourTree.h
 *
 */

#ifndef EULERTOURTREE_H_
#define EULERTOURTREE_H_

#include <iostream>
#include <limits>

using namespace std;

const int numDeltaWs = 2;

struct ETNode;
typedef TPair<ETNode *, ETNode *> ETArcPair;
typedef TVec<ETArcPair> ETArcPairs;

struct CGNode;
typedef TVec<CGNode *> CGNodeArr;

struct CGEdge;
typedef TVec<CGEdge *> AdjacencyList;

struct CGNode {

  CGNode(int m, int size_)
      : cgId(m), size(size_), adjTreeCGEdges(), loopNode(NULL) {}

  ~CGNode() {}

public:
  // the ID of this CG node
  int cgId;

  // this size of the CG node (i.e., the size of MC represented by the CG node)
  int size;

  // adjacency list of edges incident to this CGNode that are in an ET tree
  AdjacencyList adjTreeCGEdges;

  ETNode *loopNode;
};

// edges between CG nodes
struct CGEdge {

  CGEdge(CGNodeArr *cga1_, CGNodeArr *cga2_, int weight_)
      : cga1(cga1_), cga2(cga2_) {
    arcs.Gen(weight_);
  }

  ~CGEdge() {
    // note: arcs should already have been deleted when the tree was
    // cleaned up, so don't re-delete here (that memory may have already
    // been reallocated so can't do checks here either)
  }

public:
  void clear() {
    // note that we're assuming the tree gets destroyed elsewhere
    arcs.Clr();
  };

public:
  // the corresponding CGNode arrays for this edge
  CGNodeArr *cga1, *cga2;

  // list of pairs of arcs that correspond to this edge (indexed by level,
  // for all levels <= this->level)
  ETArcPairs arcs;
};

typedef TVec<ETNode *> ETNodeStack;

struct ETNode {

  ETNode(CGNode *a, CGNode *b)
      : cg1(a), cg2(b), parent(NULL), left(NULL), right(NULL) {
    // IAssert(cg1 != NULL && cg2 != NULL);
    // cout << "created ETNode: ";
    // this->print();
    // cout << endl;

    if (a->cgId == b->cgId) {
      dw = a->cgId; // dRepID
    } else {
      dw = 0;
    }
  }

  ~ETNode() {
    // cout << "deleting ETNode: ";
    // this->print();
    // cout << endl;
    // IAssert(isDisconnected());
  }

  void set(CGNode *a, CGNode *b) {
    cg1 = a;
    cg2 = b;
  }

public:
  // return true if this is a loop node
  bool isLoop() const { return (cg1 == cg2); }

  bool operator==(const ETNode &rhs) const {
    return (cg1 == rhs.cg1 && cg2 == rhs.cg2);
  }

  // return the root of this node's tree
  const ETNode *findRoot() const {
    const ETNode *z = this;
    while (z->parent != NULL) {
      z = z->parent;
    }
    return z;
  }

  ETNode *findRoot() {
    ETNode *z = this;
    while (z->parent) {
      z = z->parent;
    }
    return z;
  }

  void print(ostream &os = cout) const {
    os << "(";

    if (isLoop())
      os << cg1->cgId;
    else
      os << cg1->cgId << "," << cg2->cgId;

    os << ")";
  }

  void clear() {
    // IAssert(isDisconnected());
    cg1 = NULL;
    cg2 = NULL;
  }

  void disconnect() {
    // NOTE: This doesn't change the nodes connected to this node!
    parent = left = right = NULL;
  }

  bool isDisconnected() const {
    return (parent == NULL && left == NULL && right == NULL);
  }

public:
  // split to the left of this node and return the new tree (which was
  // the left subtree of this node)
  ETNode *splitLeft() {
    ETNode *z(left);
    left = NULL;
    if (z) {
      z->parent = NULL;
      z->dw += dw;
    }
    return z;
  }

  // split to the right of this node and return the new tree (which was
  // the right subtree of this node)
  ETNode *splitRight() {
    ETNode *z(right);
    right = NULL;
    if (z) {
      z->parent = NULL;
      z->dw += dw;
    }
    return z;
  }

  void joinLeft(ETNode *leftRoot) {
    // IAssert(!left && !leftRoot->parent);
    left = leftRoot;
    if (left) {
      left->parent = this;
      left->dw -= dw;
    }
  }

  void joinRight(ETNode *rightRoot) {
    // IAssert(!right && !rightRoot->parent);
    right = rightRoot;
    if (right) {
      right->parent = this;
      right->dw -= dw;
    }
  }

  // return the next node following this one (or null if there aren't any)
  // in an in-order traversal of the tree
  ETNode *next() const {
    ETNode *nxt(NULL);

    if (right) {
      nxt = right->subtree_minimum();
    } else if (parent) {
      if (parent->left == this) {
        nxt = parent;
      } else {
        // IAssert(parent->right == this);
        for (ETNode *n = parent; n->parent; n = n->parent) {
          if (n == n->parent->left) {
            nxt = n->parent;
            break;
          }
        }
      }
    }

    return nxt;
  }

  // return the node that precedes the specified node (or null) in an
  // in-order traversal of the tree
  ETNode *prev() const {

    ETNode *prv(NULL);

    if (left) {
      prv = left->subtree_maximum();
    } else if (parent) {
      if (parent->right == this) {
        prv = parent;
      } else {
        // IAssert(parent->left == this);
        for (ETNode *n = parent; n->parent; n = n->parent) {
          if (n == n->parent->right) {
            prv = n->parent;
            break;
          }
        }
      }
    }
    return prv;
  }

  // descend the left branch until we find the minimum node
  const ETNode *subtree_minimum() const {
    const ETNode *u(this);
    while (u->left) {
      u = u->left;
    }
    return u;
  }

  ETNode *subtree_minimum() {
    ETNode *u(this);
    while (u->left) {
      u = u->left;
    }
    return u;
  }

  // descend the right branch until we find the maximum node
  const ETNode *subtree_maximum() const {
    const ETNode *u(this);
    while (u->right) {
      u = u->right;
    }
    return u;
  }

  ETNode *subtree_maximum() {
    ETNode *u(this);
    while (u->right) {
      u = u->right;
    }
    return u;
  }

  // splay this node to the top of the tree that it's currently in and
  // return this node (for convenience when stringing together node ops.)
  // NOTE: This is not a very clean implementation, as splaying a node
  // makes it the new root of its tree, but this can't actually update an
  // instance of EulerTourTree if it exists, so that will no longer point to
  // the correct root if a node is splayed without that being updated, so
  // use this carefully.
  ETNode *splay() {
    ETNode *x = this;
    while (x->parent) {
      if (!x->parent->parent) {
        if (x->parent->left == x)
          x->parent->right_rotate();
        else
          x->parent->left_rotate();
      } else if (x->parent->left == x && x->parent->parent->left == x->parent) {
        x->parent->parent->right_rotate();
        x->parent->right_rotate();
      } else if (x->parent->right == x &&
                 x->parent->parent->right == x->parent) {
        x->parent->parent->left_rotate();
        x->parent->left_rotate();
      } else if (x->parent->left == x &&
                 x->parent->parent->right == x->parent) {
        x->parent->right_rotate();
        x->parent->left_rotate();
      } else {
        x->parent->left_rotate();
        x->parent->right_rotate();
      }
    }
    // IAssert(this->parent == NULL);
    return this;
  }

protected:
  void left_rotate() {
    ETNode *x = this;
    ETNode *y = x->right;

    int dwx, dwy;
    dwx = x->dw;
    dwy = y->dw;

    x->right = y->left;
    if (y->left) {
      y->left->parent = x;
      x->right->dw += dwy;
    }

    y->dw += dwx;
    x->dw = -dwy;

    y->parent = x->parent;
    if (x->parent) {
      if (x == x->parent->left)
        x->parent->left = y;
      else
        x->parent->right = y;
    }

    y->left = x;
    x->parent = y;
  }

  void right_rotate() {
    ETNode *x = this;
    ETNode *y = x->left;

    int dwx, dwy;
    dwx = x->dw;
    dwy = y->dw;

    x->left = y->right;
    if (y->right) {
      y->right->parent = x;
      x->left->dw += dwy;
    }

    y->dw += dwx;
    x->dw = -dwy;

    y->parent = x->parent;

    if (x->parent) {
      if (x == x->parent->left)
        x->parent->left = y;
      else
        x->parent->right = y;
    }

    y->right = x;
    x->parent = y;
  }

public:
  // pointers to the two CG nodes that are at the ends of the arc (directed
  // edge) that this node represents in the Euler tour
  CGNode *cg1, *cg2;

  // pointers to the attached nodes in the tree
  ETNode *parent, *left, *right;

  int dw;
}; // struct ETNode

class EulerTourTree {

public:
  EulerTourTree() : root(NULL) {}

  ~EulerTourTree() {}

  static void createTour(CGNode *v) {
    // IAssert( v != NULL && v->loopNode == NULL );
    v->loopNode = new ETNode(v, v);
  }

  // delete all of the nodes in this tree
  static void destroyTour(ETNode *root) {
    EulerTourTree::postOrderDelete(root);
    ;
  }

  static void postOrderDelete(ETNode *root) {
    if (!root)
      return;
    ETNode *cur = root;
    ETNode *prev = NULL;

    while (cur) {
      if (prev == cur->parent || prev == NULL) { // traversing down the tree
        prev = cur;
        if (cur->left)
          cur = cur->left;
        else if (cur->right)
          cur = cur->right;
        else
          cur = cur->parent;          // reached the bottom
      } else if (prev == cur->left) { // traversing up from the left
        // perform the visit of the node we are coming from (left)
        //				( this->*(onVisit) )( prev );
        prev->disconnect();
        delete prev;

        prev = cur;
        if (cur->right)
          cur = cur->right;
        else
          cur = cur->parent;
      } else if (prev == cur->right) { // traversing up from the right
        // perform the visit of the node we are coming from (right)
        prev->disconnect();
        delete prev;

        prev = cur;
        cur = cur->parent;
      } else {
        // IAssert(false);
      }
    }

    // finally, visit the root
    // IAssert(prev == root);
    prev->disconnect();
    delete prev;
  }

  static void cut(ETNode *n1, ETNode *n2, ETNode **r1, ETNode **r2) {
    // IAssert(n1->findRoot() == n2->findRoot());
    // IAssert(!n1->isLoop() && !n2->isLoop());
    // IAssert(n1->cg1 == n2->cg2 && n1->cg2 == n2->cg1);

    // cut out the two edges in the tour (represented by n1 and n2)

    n1->splay();
    // IAssert(n1->parent == NULL);
    ETNode *l1 = (ETNode *)n1->splitLeft();
    ETNode *l2 = (ETNode *)n1->splitRight();

    // get the actual previous/next node (not just the root of the pre-split
    // subtree)
    if (l1 != NULL)
      l1 = (ETNode *)l1->subtree_maximum();
    if (l2 != NULL)
      l2 = (ETNode *)l2->subtree_minimum();

    n2->splay();
    // IAssert(n2->parent == NULL);
    ETNode *l3 = (ETNode *)n2->splitLeft();
    ETNode *l4 = (ETNode *)n2->splitRight();

    // get the actual previous/next node (not just the root of the pre-split
    // subtree)
    if (l3 != NULL)
      l3 = (ETNode *)l3->subtree_maximum();
    if (l4 != NULL)
      l4 = (ETNode *)l4->subtree_minimum();

    // check whether n1 and n2 were in the correct order or reversed, such
    // that their order in the tour was actually ..., n2, ..., n1, ...
    if (!l2 || !l3 || l2->findRoot() != l3->findRoot()) {
      // IAssert(l1->findRoot() == l4->findRoot());
      Swap(l1, l3);
      Swap(l2, l4);
    }

    // join l1 and l4 (l2 and l3 should already be the same subtree)
    l1 = join(l1, l4);

    // delete these nodes that are no longer needed
    delete n1;
    delete n2;

    // return the roots of the two trees
    *r1 = l1->findRoot();
    *r2 = l3->findRoot();
  }

  static ETNode *link(CGNode *u, CGNode *v, ETNode **neuv, ETNode **nevu) {
    ETNode *nluu = u->loopNode;
    ETNode *nlvv = v->loopNode;

    // IAssert(nluu != NULL && nlvv != NULL);
    // IAssert(!connected(nluu, nlvv));
    // IAssert(*neuv == NULL && *nevu == NULL);
    // split each tour (represented as a list embedded in a BST), after each
    // loop node, to get L1 -> (L11, L12) and L2 -> (L21, L22)

    nluu->splay();
    ETNode *l11 = nluu;
    ETNode *l12 = (ETNode *)nluu->splitRight();

    nlvv->splay();
    ETNode *l21 = nlvv;
    ETNode *l22 = (ETNode *)nlvv->splitRight();

    // create the new nodes for the two arcs representing the new edge
    *neuv = new ETNode(u, v);
    *nevu = new ETNode(v, u);

    // concat as: ( l12, l11, neuv, l22, l21, nevu ) -- takes 5 joins
    ETNode *tmp1 = join(l11, *neuv);
    ETNode *tmp2 = join(l21, *nevu);

    ETNode *rl = join(l12, tmp1);
    ETNode *rr = join(l22, tmp2);

    return join(rl, rr);
  }

  static bool connected(const ETNode *u, const ETNode *v) {
    return (u->findRoot() == v->findRoot());
  }

  static bool connected(const CGNode *u, const CGNode *v) {
    return connected(u->loopNode, v->loopNode);
  }

  static void print(const ETNode *n) {

    int count(0);
    n = (ETNode *)n->findRoot()->subtree_minimum();
    // IAssert(n != NULL);

    while (n != NULL) {
      cout << (count > 0 ? ", " : "");
      n->print(cout);
      // IAssert(n->checkSizeAndWeight());
      ++count;
      n = (ETNode *)n->next();
    }
  }

  static ETNode *join(ETNode *p, ETNode *q) {
    if (p == NULL)
      return q->findRoot();
    if (q == NULL)
      return p->findRoot();
    p = (ETNode *)p->findRoot()->subtree_maximum()->splay();
    // IAssert(p->right == NULL);
    p->joinRight(q->findRoot());
    return p;
  }

  static int findRepID(ETNode *cur) {
    ETNode *n = cur;
    int sum = 0;
    while (n->parent != NULL) {
      sum += n->dw;
      n = n->parent;
    }
    sum += n->dw;
    return sum;
  }

protected:
  // the root node of this tree
  ETNode *root;
};

#endif // EULERTOURTREE_H_