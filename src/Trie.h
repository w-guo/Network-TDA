/*
 * trie.h
 *
 *  Created on: Jan 10, 2019
 *      Author: robustor
 */

#ifndef TRIE_H_
#define TRIE_H_

#include <iostream>
#include <vector>
#include <map>
#include <list>

using namespace std;
typedef vector<int> MCVec;
typedef map<int, MCVec> MCPair;

// define a global variable to save the multiple outputs from the recursive function
MCPair * MCP = new MCPair();

class TrieNode {

public:
	TrieNode(TrieNode * parent, int vid)
		: key(vid)
		, value(-1)
		, par(parent)
		, chi()
	{}

	~TrieNode() {
//		if (!chi.empty()) {
//			chi.erase(chi.begin(), chi.end());
//		}
		par = NULL;
		chi.clear();
	}

	int key, value;
	TrieNode * par;
	map<int, TrieNode *> chi;
};

class Trie {

public:
	Trie(int vid){
		root = new TrieNode(NULL, vid);
		NodeList[vid].push_back(root);
	}

	~Trie() {			// delete all MCs in the trie
		destroy(root);
		NodeList.clear();
	}

	void destroy(TrieNode * cur);
	void addMC(MCVec * MC, int m);
	void removeMC(MCVec * MC);
	 int getID(MCVec * MC);
	void searchDown(TrieNode * cur, MCVec * MC);
	MCPair* outputMC(TrieNode * cur);

public:
	TrieNode * root;
	map<int, list<TrieNode *>> NodeList;
};

//Delete the trie following post-order traversal
void Trie::destroy(TrieNode * cur) {
	if (cur == NULL)
		return;
	if (!cur->chi.empty()) {
		for (auto it : cur->chi) {
			destroy(it.second);
		}
	}
	delete cur;
}


void Trie::addMC(MCVec * MC, int m) {
	// start from root node
	TrieNode * cur = root;

	for (auto it = MC->begin()+1; it != MC->end(); it++) {
		// create a new node if path doesn't exists
		if (cur->chi.find(*it) == cur->chi.end()) {
			cur->chi[*it] = new TrieNode(cur, *it);
			// add this new node to the NodeList
			NodeList[*it].push_back(cur->chi[*it]);
		}
		// go to next node
		cur = cur->chi[*it];
	}
	// current node is now the leaf node
	cur->value = m;
}


void Trie::removeMC(MCVec * MC) {
	TrieNode * cur = root;
	// move down to the leaf node of this MC
	for (auto it = MC->begin()+1; it != MC->end(); it++) {
		if (cur->chi.find(*it) == cur->chi.end())
			return;
		cur = cur->chi[*it];
	}

	// walk up looking for nodes to remove
	while (cur->chi.empty()) {
		int cKey = cur->key;
		// remove this node from the NodeList
		NodeList[cKey].remove(cur);
		if (NodeList[cKey].empty()) {
			NodeList.erase(cKey);
		}
		// remove this node from this trie
		if (cur->par != NULL) {
			cur = cur->par;
			delete cur->chi[cKey];
			cur->chi[cKey] = NULL;
			cur->chi.erase(cKey);
		}
		else {
			delete cur;
			cur = NULL;
			return;
		}
	}

}


int Trie::getID(MCVec * MC) {
	TrieNode * cur = root;
	// move down to the leaf node of this MC
	for (auto it = MC->begin()+1; it != MC->end(); it++) {
		if (cur->chi.find(*it) == cur->chi.end())
			return -1;
		cur = cur->chi[*it];
	}
	return cur->value;
}


void Trie::searchDown(TrieNode * cur, MCVec * MC) {
	if (cur->chi.empty()) {
		MCP->insert(make_pair(cur->value, *MC));
		return;
		}
	else {
		for (auto it : cur->chi) {
			MC->push_back(it.first);
			searchDown(it.second, MC);
			MC->pop_back();
		}
	}
}


MCPair * Trie::outputMC(TrieNode * cur) {
	MCP->clear();
	MCVec * MC = new MCVec();
	TrieNode * pre = cur;
	while (pre != NULL) {
		MC->push_back(pre->key);
		pre = pre->par;
	}
	reverse(MC->begin(), MC->end());
	searchDown(cur, MC);

	return MCP;
}

#endif /* TRIE_H_ */
