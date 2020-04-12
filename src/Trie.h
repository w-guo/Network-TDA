/*
 * trie.h
 *
 *  Created on: Jan 10, 2019
 *      Author: robustor
 */

#ifndef TRIE_H_
#define TRIE_H_

#include <iostream>

using namespace std;

// define a global variable to save the multiple outputs from the recursive
// function
TIntIntVH *MCVH = new TIntIntVH();

class TrieNode {

public:
  TrieNode(TrieNode *parent, int vId)
      : key(vId), value(-1), par(parent), next(NULL), chi() {}

  ~TrieNode() {
    par = NULL;
    next = NULL;
    chi.Clr();
  }

  int key, value;
  TrieNode *par, *next;
  THash<TInt, TrieNode *> chi;
};

class Trie {

public:
  Trie(int vId) {
    root = new TrieNode(NULL, vId);
    NodeList.AddDat(vId, root);
  }

  ~Trie() { // delete all MCs in the trie
    destroy(root);
    NodeList.Clr();
  }

  void destroy(TrieNode *cur);
  void addMC(TIntV &MCV, int m);
  void removeMC(TIntV &MCV);
  int getID(TIntV &MCV);
  void searchDown(TrieNode *cur, TIntV &MCV);
  TIntIntVH *outputMC(TrieNode *cur);

public:
  TrieNode *root;
  THash<TInt, TrieNode *> NodeList;
};

// Delete the trie following post-order traversal
void Trie::destroy(TrieNode *cur) {
  if (cur == NULL)
    return;
  if (!cur->chi.Empty()) {
    for (THashKeyDatI<TInt, TrieNode *> it = cur->chi.BegI(); !it.IsEnd();
         it++) {
      destroy(it.GetDat());
    }
  }
  delete cur;
}

void Trie::addMC(TIntV &MCV, int m) {
  // start from root node
  TrieNode *cur = root;

  for (int i = 1; i < MCV.Len(); i++) {
    int vId = MCV[i];
    // create a new node if path doesn't exists
    if (!cur->chi.IsKey(vId)) {
      TrieNode *node = new TrieNode(cur, vId);
      cur->chi.AddDat(vId, node);
      // add this new node to the NodeList
      if (!NodeList.IsKey(vId)) {
        NodeList.AddDat(vId, node);
      } else {
        node->next = NodeList.GetDat(vId);
        NodeList.GetDat(vId) = node;
      }
    }
    // go to next node
    cur = cur->chi.GetDat(vId);
  }
  // current node is now the leaf node
  cur->value = m;
}

void Trie::removeMC(TIntV &MCV) {
  TrieNode *cur = root;
  // move down to the leaf node of this MCV
  for (int i = 1; i < MCV.Len(); i++) {
    int vId = MCV[i];
    if (!cur->chi.IsKey(vId))
      return;
    cur = cur->chi.GetDat(vId);
  }
  cur->value = -1;
  // walk up looking for nodes to remove
  while (cur->chi.Empty()) {
    int cKey = cur->key;
    // remove this node from the NodeList
    TrieNode *node = NodeList.GetDat(cKey);
    if (node == cur) {
      if (node->next == NULL) {
        NodeList.DelKey(cKey);
      } else {
        NodeList.GetDat(cKey) = node->next;
      }
    }
    while (node != NULL) {
      if (node->next == cur) {
        node->next = cur->next;
        break;
      }
      node = node->next;
    }
    // remove this node from this trie
    if (cur->par != NULL) {
      cur = cur->par;
      delete cur->chi.GetDat(cKey);
      cur->chi.GetDat(cKey) = NULL;
      cur->chi.DelKey(cKey);
    } else {
      delete cur;
      cur = NULL;
      return;
    }
  }
}

int Trie::getID(TIntV &MCV) {
  TrieNode *cur = root;
  // move down to the leaf node of this MCV
  for (int i = 1; i < MCV.Len(); i++) {
    int vId = MCV[i];
    if (!cur->chi.IsKey(vId))
      return -1;
    cur = cur->chi.GetDat(vId);
  }
  return cur->value;
}

void Trie::searchDown(TrieNode *cur, TIntV &MCV) {
  if (cur->chi.Empty()) {
    MCVH->AddDat(cur->value, MCV);
    return;
  } else {
    for (THashKeyDatI<TInt, TrieNode *> it = cur->chi.BegI(); !it.IsEnd();
         it++) {
      MCV.Add(it.GetKey());
      searchDown(it.GetDat(), MCV);
      MCV.DelLast();
    }
  }
}

TIntIntVH *Trie::outputMC(TrieNode *cur) {
  MCVH->Clr();
  TIntV MCV;
  TrieNode *pre = cur;
  while (pre != NULL) {
    MCV.Add(pre->key);
    pre = pre->par;
  }
  MCV.Reverse();
  searchDown(cur, MCV);
  return MCVH;
}

#endif /* TRIE_H_ */