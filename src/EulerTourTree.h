/*
 * EulerTourTree.h
 *
 */

#ifndef EULERTOURTREE_H_
#define EULERTOURTREE_H_

#include <iostream>
#include <vector>
#include <map>
using namespace std;

const int numDeltaWs = 3;
const int numDeltaMinWs = 2;

struct ETNode;
typedef pair<ETNode *, ETNode *> ETArcPair;
typedef vector<ETArcPair> ETArcPairs;


struct CGNode;
typedef vector<CGNode *> CGNodeArr;

struct CGEdge;
typedef vector<CGEdge *> AdjacencyList;

struct CGNode {

	CGNode(int m, int size_)
	: cgid(m)
	, size(size_)
	, adjTreeCGEdges()
	, adjNontreeCGEdges()
	, loopNode(NULL)
	{}

	~CGNode() {}

public:

	static CGNodeArr * initialize(int m, int nlevels);

public:
	// the ID of this CG node
	int cgid;

	// this size of the CG node (i.e., the size of MC represented by the CG node)
	int size;

	// adjacency list of edges incident to this CGNode that are in an ET tree
	AdjacencyList adjTreeCGEdges;

	// adjacency list containing the incident edges that are not in an ET tree
	AdjacencyList adjNontreeCGEdges;

	ETNode * loopNode;
};


// edges between CG nodes
struct CGEdge {

	CGEdge(CGNodeArr * cga1_, CGNodeArr * cga2_, int weight_)
		: cga1(cga1_)
		, cga2(cga2_)
		, weight(weight_)
	{
		arcs.reserve(1);
	}

	~CGEdge() {
		// note: arcs should already have been deleted when the tree was
		// cleaned up, so don't re-delete here (that memory may have already
		// been reallocated so can't do checks here either)
	}

public:

	void clear() {
		weight = -1;

		// note that we're assuming the tree gets destroyed elsewhere
		arcs.clear();
	};

	/*
		public:
			bool isequal(int uu, int vv) const {
				return ((cgid1 == uu && cgid2 == vv) ||
						(cgid1 == vv && cgid2 == uu));
			}

			bool operator==(const CGEdge & rhs) const {
				bool res = isequal( rhs.cgid1, rhs.cgid2 );
				return res;
			}
	 */
public:
	// the corresponding CGNode arrays for this edge
	CGNodeArr * cga1, * cga2;

	// the weight of this edge
	int weight;

	// list of pairs of arcs that correspond to this edge (indexed by level,
	// for all levels <= this->level)
	ETArcPairs arcs;
};

typedef vector<ETNode *> ETNodeStack;

struct ETNode {

	ETNode(CGNode * a, CGNode * b)
		: cg1(a)
		, cg2(b)
		, parent(NULL)
		, left(NULL)
		, right(NULL)
		, stsize(1)
		, weight(0)
		, stweight(weight)
	{
		assert( cg1 != NULL && cg2 != NULL );
		cout << "created ETNode: "; this->print(); cout << endl;

		for ( int i = 0; i < numDeltaMinWs; ++i ) {
			dminw[i] = 0; //dMinID: dminw[0]; dMinV: dminw[1]
		}

		if (a->cgid == b->cgid){
			dw[0] = a->cgid;  //dID
			dw[1] = -a->size; //dV
			dw[2] = a->cgid;  //dRepID
		}
		else{
			dw[0] = numeric_limits<int>::max();
			dw[1] = 1;
			dw[2] = 0;
		}
	}

	~ETNode() {
		cout << "deleting ETNode: "; this->print(); cout << endl;
		assert( isDisconnected() );
	}

	void set(CGNode * a, CGNode * b) {
		cg1 = a;
		cg2 = b;
	}

public:

	// return true if this is a loop node
	bool isLoop() const { return ( cg1 == cg2 ); }

	bool operator==( const ETNode & rhs ) const {
		return ( cg1 == rhs.cg1 && cg2 == rhs.cg2 );
	}

	// return the root of this node's tree
	const ETNode * findRoot() const {
		const ETNode * z = this;
		while ( z->parent != NULL ) { z = z->parent; }
		return z;
	}

	ETNode * findRoot() {
		ETNode * z = this;
		while ( z->parent ) { z = z->parent; }
		return z;
	}

	void print( ostream & os = cout ) const {
		os << "(";

		if (isLoop()) os << cg1->cgid;
		else os << cg1->cgid << "," << cg2->cgid;

		os << ")";
	}

	void clear() {
		assert( isDisconnected() );
		cg1 = NULL;
		cg2 = NULL;
	}

	void disconnect() {
		// NOTE: This doesn't change the nodes connected to this node!
		parent = left = right = NULL;
	}

	bool isDisconnected() const {
		return ( parent == NULL && left == NULL && right == NULL );
	}

public:
	// split to the left of this node and return the new tree (which was
	// the left subtree of this node)
	ETNode * splitLeft() {
		ETNode * z( left );
		left = NULL;
		if ( z ) {
			z->parent = NULL;
			stsize -= z->stsize;
			stweight -= z->stweight;
			for ( int i = 0; i < numDeltaWs; ++i ) {
				z->dw[i] += dw[i];
			}
			updateDeltaMinWeight(this);
		}
		return z;
	}

	// split to the right of this node and return the new tree (which was
	// the right subtree of this node)
	ETNode * splitRight() {
		ETNode * z( right );
		right = NULL;
		if ( z ) {
			z->parent = NULL;
			stsize -= z->stsize;
			stweight -= z->stweight;
			for ( int i = 0; i < numDeltaWs; ++i ) {
				z->dw[i] += dw[i];
			}
			updateDeltaMinWeight(this);
		}
		return z;
	}

	void joinLeft( ETNode * leftRoot ) {
		assert( !left && !leftRoot->parent );
		left = leftRoot;
		if ( left ) {
			stsize += left->stsize;
			stweight += left->stweight;
			left->parent = this;
			for ( int i = 0; i < numDeltaWs; ++i ) {
				left->dw[i] -= dw[i];
			}
			updateDeltaMinWeight(this);
		}
	}

	void joinRight( ETNode * rightRoot ) {
		assert( !right && !rightRoot->parent );
		right = rightRoot;
		if ( right ) {
			stsize += right->stsize;
			stweight += right->stweight;
			right->parent = this;
			for ( int i = 0; i < numDeltaWs; ++i ) {
				right->dw[i] -= this->dw[i];
			}
			updateDeltaMinWeight(this);
		}
	}

	// return the next node following this one (or null if there aren't any)
	// in an in-order traversal of the tree
	ETNode * next() const {
		ETNode * nxt( NULL );

		if ( right ) {
			nxt = right->subtree_minimum();
		} else if ( parent ) {
			if ( parent->left == this ) {
				nxt = parent;
			} else {
				assert( parent->right == this );
				for ( ETNode * n = parent; n->parent; n = n->parent ) {
					if ( n == n->parent->left ) {
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
	ETNode * prev() const {

		ETNode * prv( NULL );

		if ( left ) {
			prv = left->subtree_maximum();
		} else if ( parent ) {
			if ( parent->right == this ) {
				prv = parent;
			} else {
				assert( parent->left == this );
				for ( ETNode * n = parent; n->parent; n = n->parent ) {
					if ( n == n->parent->right ) {
						prv = n->parent;
						break;
					}
				}
			}
		}

		return prv;
	}

	// return the next node in a pre-order traversal of this node's tree,
	// using the stack to maintain state across calls to this function (can
	// also use to traverse only those nodes with positive weights)
	ETNode * nextPreOrder(ETNodeStack & stack, bool onlyWeighted = false) const {

		ETNode * nxt( NULL );

		// remember the right node for later
		if ( right != NULL && ( !onlyWeighted ||
				right->stweight > 0 ) ) {
			stack.push_back( right );
		}

		// try to move to the left node now, or else get one from the stack
		if ( left != NULL &&
				( !onlyWeighted || left->stweight > 0 ) ) {
			nxt = left;

		} else if ( !stack.empty() ) {
			nxt = stack.back();
			stack.pop_back();
		}

		assert( nxt != NULL || stack.empty() );

		return nxt;
	}

	// descend the left branch until we find the minimum node
	const ETNode * subtree_minimum() const {
		const ETNode * u( this );
		while ( u->left ) { u = u->left; }
		return u;
	}

	ETNode * subtree_minimum() {
		ETNode * u( this );
		while ( u->left ) { u = u->left; }
		return u;
	}

	// descend the right branch until we find the maximum node
	const ETNode * subtree_maximum() const {
		const ETNode * u( this );
		while ( u->right ) { u = u->right; }
		return u;
	}

	ETNode * subtree_maximum() {
		ETNode * u( this );
		while ( u->right ) { u = u->right; }
		return u;
	}

	void changeWeight( int wnew ) {
		int diff( wnew - weight );
		weight = wnew;
		stweight += diff;

		assert( weight >= 0 && stweight >= 0 );
		assert( checkSizeAndWeight() );

		// walk up the tree and adjust all subtree weights accordingly
		ETNode * n( this );
		while ( ( n = n->parent ) ) {
			n->stweight += diff;
			n->checkSizeAndWeight();
		}
	}

	void addWeight( int wdelta ) {
		changeWeight( weight + wdelta );
	}

	inline bool checkSizeAndWeight() const {
#ifdef DEBUG
		bool sz( stsize == ( 1 + ( left ? left->stsize : 0 ) +
				( right ? right->stsize : 0 ) ) );

		bool wt( true );

		wt = wt && ( stweight[i] == ( weight[i] +
				( left ? left->stweight : 0 ) +
				( right ? right->stweight : 0 ) ) );

		assert( weight >= 0 && weight <= 1000000 );


		return ( sz && wt );
#else
		return true;
#endif // DEBUG
	}

	// splay this node to the top of the tree that it's currently in and
	// return this node (for convenience when stringing together node ops.)
	// NOTE: This is not a very clean implementation, as splaying a node
	// makes it the new root of its tree, but this can't actually update an
	// instance of EulerTourTree if it exists, so that will no longer point to
	// the correct root if a node is splayed without that being updated, so
	// use this carefully.
	ETNode * splay() {
		ETNode * x = this;
		while (x->parent) {
			if (!x->parent->parent) {
				if (x->parent->left == x) x->parent->right_rotate();
				else x->parent->left_rotate();

			} else if (x->parent->left == x
					&& x->parent->parent->left == x->parent) {
				x->parent->parent->right_rotate();
				x->parent->right_rotate();

			} else if (x->parent->right == x
					&& x->parent->parent->right == x->parent) {
				x->parent->parent->left_rotate();
				x->parent->left_rotate();

			} else if (x->parent->left == x
					&& x->parent->parent->right == x->parent) {
				x->parent->right_rotate();
				x->parent->left_rotate();

			} else {
				x->parent->left_rotate();
				x->parent->right_rotate();
			}
		}

		assert( checkSizeAndWeight() );
		assert( this->parent == NULL );
		return this;
	}

	void updateDeltaMinWeight(ETNode * u) {
		for ( int i = 0; i < numDeltaMinWs; ++i ) {
			u->dminw[i] = max(0, u->left ? u->left->dminw[i] - u->left->dw[i] : 0);
			u->dminw[i] = max(u->dminw[i], u->right ? u->right->dminw[i] - u->right->dw[i] : 0);
		}
	}

protected:
	void left_rotate() {
		ETNode *x = this;
		ETNode *y = x->right;

		int dwx[numDeltaWs], dwy[numDeltaWs];
		for ( int i = 0; i < numDeltaWs; ++i ) {
			dwx[i] = x->dw[i];
			dwy[i] = y->dw[i];
		}

		x->right = y->left;
		if (y->left) {
			y->left->parent = x;
			for ( int i = 0; i < numDeltaWs; ++i ) {
				x->right->dw[i] += dwy[i];
			}
		}

		for ( int i = 0; i < numDeltaWs; ++i ) {
			y->dw[i] += dwx[i];
			x->dw[i] = -dwy[i];
		}

		y->parent = x->parent;
		if ( x->parent ) {
			if (x == x->parent->left) x->parent->left = y;
			else x->parent->right = y;
		}

		// update sizes before setting y->left
		x->stsize -= y->stsize - ( y->left ? y->left->stsize : 0 );
		y->stsize += 1 + ( x->left ? x->left->stsize : 0 );

		y->left = x;
		x->parent = y;

		// update weights
		x->stweight = x->weight +
				( x->left ? x->left->stweight : 0 ) +
				( x->right ? x->right->stweight : 0 );
		y->stweight = y->weight +
				( y->left ? y->left->stweight : 0 ) +
				( y->right ? y->right->stweight : 0 );

		assert( !x || x->checkSizeAndWeight() );
		assert( !y || y->checkSizeAndWeight() );
		assert( !y->parent || y->parent->checkSizeAndWeight() );

		updateDeltaMinWeight(x);
		updateDeltaMinWeight(y);
	}

	void right_rotate() {
		ETNode *x = this;
		ETNode *y = x->left;

		int dwx[numDeltaWs], dwy[numDeltaWs];
		for ( int i = 0; i < numDeltaWs; ++i ) {
			dwx[i] = x->dw[i];
			dwy[i] = y->dw[i];
		}

		x->left = y->right;
		if (y->right) {
			y->right->parent = x;
			for ( int i = 0; i < numDeltaWs; ++i ) {
				x->left->dw[i] += dwy[i];
			}

		}

		for ( int i = 0; i < numDeltaWs; ++i ) {
			y->dw[i] += dwx[i];
			x->dw[i] = -dwy[i];
		}

		y->parent = x->parent;

		if ( x->parent ) {
			if (x == x->parent->left) x->parent->left = y;
			else x->parent->right = y;
		}

		// update the sizes before setting y->right
		x->stsize -= y->stsize - ( y->right ? y->right->stsize : 0 );
		y->stsize += 1 + ( x->right ? x->right->stsize : 0 );

		y->right = x;
		x->parent = y;

		// update weights after rotation
		x->stweight = x->weight +
				( x->left ? x->left->stweight : 0 ) +
				( x->right ? x->right->stweight : 0 );
		y->stweight = y->weight +
				( y->left ? y->left->stweight : 0 ) +
				( y->right ? y->right->stweight : 0 );

		assert( !x || x->checkSizeAndWeight() );
		assert( !y || y->checkSizeAndWeight() );
		assert( !y->parent || y->parent->checkSizeAndWeight() );

		updateDeltaMinWeight(x);
		updateDeltaMinWeight(y);
	}

public:
	// pointers to the two CG nodes that are at the ends of the arc (directed
	// edge) that this node represents in the Euler tour
	CGNode * cg1, * cg2;

	// pointers to the attached nodes in the tree
	ETNode * parent, * left, * right;

	// the size of the subtree rooted at this node (includes this node)
	size_t stsize;

	// the weight(s) of this node
	int weight;

	// the weight(s) of the subtree rooted at this node (includes this node)
	int stweight;

	int dw[numDeltaWs];
	int dminw[numDeltaMinWs];

}; // struct ETNode


class EulerTourTree {

public:
	EulerTourTree()
		: root(NULL)
	{}

	~EulerTourTree() {}

	static void createTour(CGNode * v) {
		//assert( v != NULL && v->loopNode == NULL );
		v->loopNode = new ETNode(v, v);
	}

	// delete all of the nodes in this tree
	static void destroyTour(ETNode * root) {
		EulerTourTree::postOrderDelete( root );;
	}

	// join <rightTree> to this tree by splaying the max node and then inserting
	// <rightTree> as the right subtree of the splayed node
//	void join( ETNode * rightRoot ) {
//		ETNode * z = root->subtree_maximum();
//		z->splay();
//		root = z;
//
//		assert( !rightRoot->parent );
//		assert( !z->right );
//		z->joinRight( rightRoot );
//	}


	static void postOrderDelete(ETNode * root) {
		if (!root) return;
		ETNode * cur = root;
		ETNode * prev = NULL;

		while (cur) {
			if ( prev == cur->parent || prev == NULL ) { // traversing down the tree
				prev = cur;
				if ( cur->left )
					cur = cur->left;
				else if ( cur->right )
					cur = cur->right;
				else
					cur = cur->parent; // reached the bottom

			} else if ( prev == cur->left ) { // traversing up from the left
				// perform the visit of the node we are coming from (left)
//				( this->*(onVisit) )( prev );
				prev->disconnect();
				delete prev;

				prev = cur;
				if ( cur->right ) cur = cur->right;
				else cur = cur->parent;

			} else if ( prev == cur->right ) { // traversing up from the right
				// perform the visit of the node we are coming from (right)
//				( this->*(onVisit) )( prev );
				prev->disconnect();
				delete prev;

				prev = cur;
				cur = cur->parent;
			} else {
				assert( false );
			}
		}

		// finally, visit the root
		assert( prev == root );
		prev->disconnect();
		delete prev;
	}


	static void cut(ETNode * n1, ETNode * n2, ETNode ** r1, ETNode ** r2) {
		assert( n1->findRoot() == n2->findRoot() );
		assert( !n1->isLoop() && !n2->isLoop() );
		assert( n1->cg1 == n2->cg2 && n1->cg2 == n2->cg1 );

		// cut out the two edges in the tour (represented by n1 and n2)

		n1->splay();
		assert( n1->parent == NULL );
	//	cout << "splitting 1 "; n1->print(); cout << " from ";
	//	if ( n1->left ) n1->left->print();
	//	else cout << "null";
	//	cout << " and ";
	//	if ( n1->right ) n1->right->print();
	//	else cout << "null";
	//	cout << endl;
		ETNode * l1 = (ETNode *) n1->splitLeft();
		ETNode * l2 = (ETNode *) n1->splitRight();

		// get the actual previous/next node (not just the root of the pre-split subtree)
		if ( l1 != NULL ) l1 = (ETNode *) l1->subtree_maximum();
		if ( l2 != NULL ) l2 = (ETNode *) l2->subtree_minimum();

		n2->splay();
		assert( n2->parent == NULL );
	//	cout << "splitting 2 "; n2->print(); cout << " from ";
	//	if ( n2->left ) n2->left->print();
	//	else cout << "null";
	//	cout << " and ";
	//	if ( n2->right ) n2->right->print();
	//	else cout << "null";
	//	cout << endl;
		ETNode * l3 = (ETNode *) n2->splitLeft();
		ETNode * l4 = (ETNode *) n2->splitRight();

		// get the actual previous/next node (not just the root of the pre-split subtree)
		if (l3 != NULL) l3 = (ETNode *) l3->subtree_maximum();
		if (l4 != NULL) l4 = (ETNode *) l4->subtree_minimum();

		// check whether n1 and n2 were in the correct order or reversed, such that
		// their order in the tour was actually ..., n2, ..., n1, ...
		if ( !l2 || !l3 || l2->findRoot() != l3->findRoot() ) {
			assert( l1->findRoot() == l4->findRoot() );
			swap(l1, l3);
			swap(l2, l4);
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


	static ETNode * link(CGNode * u, CGNode * v, ETNode ** neuv, ETNode ** nevu ) {
		ETNode * nluu = u->loopNode;
		ETNode * nlvv = v->loopNode;

		assert( nluu != NULL && nlvv != NULL );
		assert( !connected( nluu, nlvv ) );
		assert( *neuv == NULL && *nevu == NULL );
		// split each tour (represented as a list embedded in a BST), after each
		// loop node, to get L1 -> (L11, L12) and L2 -> (L21, L22)

		nluu->splay();
		ETNode * l11 = nluu;
		ETNode * l12 = (ETNode *) nluu->splitRight();

		nlvv->splay();
		ETNode * l21 = nlvv;
		ETNode * l22 = (ETNode *) nlvv->splitRight();

		// create the new nodes for the two arcs representing the new edge
		*neuv = new ETNode(u, v);
		*nevu = new ETNode(v, u);

		// concat as: ( l12, l11, neuv, l22, l21, nevu ) -- takes 5 joins

		ETNode * tmp1 = join( l11, *neuv );
		ETNode * tmp2 = join( l21, *nevu );

		ETNode * rl = join( l12, tmp1 );
		ETNode * rr = join( l22, tmp2 );

		return join(rl, rr);
	}

	static bool connected(const ETNode * u, const ETNode * v) {
		return (u->findRoot() == v->findRoot());
	}

	static bool connected(const CGNode * u, const CGNode * v) {
		return connected( u->loopNode, v->loopNode);
	}

	static void print(const ETNode * n) {

		int count( 0 );
		n = (ETNode *) n->findRoot()->subtree_minimum();

		assert(n != NULL);

		while (n != NULL) {
			cout << ( count > 0 ? ", " : "" );
			n->print( cout );

			assert(n->checkSizeAndWeight());

			++count;
			n = (ETNode *) n->next();
		}
	}


	static ETNode * join(ETNode * p, ETNode * q) {
		if (p == NULL) return q->findRoot();
		if (q == NULL) return p->findRoot();
		p = (ETNode *) p->findRoot()->subtree_maximum()->splay();
		assert(p->right == NULL);
		p->joinRight(q->findRoot());
		return p;
	}


	static ETNode * findRep(ETNode * root) {
		ETNode * cur = root;
		while (cur->dminw[1] != 0) {
			if (cur->left->dminw[1] - cur->left->dw[1] > cur->right->dminw[1] - cur->right->dw[1]) {
				cur = cur->left;
			}
			else if (cur->left->dminw[1] - cur->left->dw[1] < cur->right->dminw[1] - cur->right->dw[1]) {
				cur = cur->right;
			}
			else if (cur->left->dminw[0] - cur->left->dw[0] > cur->right->dminw[0] - cur->right->dw[0]) {
				cur = cur->left;
			}
			else {
				cur = cur->right;
			}
		}
		cur->splay();
		return cur;
	}

protected:
	// the root node of this tree
	ETNode * root;
};


CGNodeArr * CGNode::initialize(int m, int nlevels) {
	CGNodeArr * cga = (CGNodeArr *)malloc((nlevels-1)*sizeof(CGNode));
	for (int i = 0; i < nlevels-1; i++) {
		CGNode * cg = new CGNode(m, nlevels);
		EulerTourTree::createTour(cg);
		//assert(cg == cg->loopNode->cg1);
		cga->emplace_back(cg);
		//cout << "cga address: " << cga << endl;
	}
	return cga;
};

#endif // EULERTOURTREE_H_
