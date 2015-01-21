#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <ctime>
#include <string.h>
#include <string>
#include <sys/time.h>
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <sys/stat.h>

bool outputTree = true;
bool relativeDist = false;
bool branchDist = false;
bool branchDistPop = false;
double distParameter = 1.0;
int minBr = 2, maxBr = 0;
int lineNum = 0;
int bDPsample = 1, lineCounter = 10000;
int numOfSites = -1;
time_t t1, t2;
double t3=0.0;


struct Cursor{
	int M; //Number of haplotypes
	int N; //sites
	int lf; //look forward, size of buffer
	int bufP; //buffer pointer
	int numBrRec, numSiteRec;
	int treeSize;
	std::vector<int> countBrRec;
	std::vector< std::vector<int> > buf; //buffer of sites
	std::vector<int> y; //Alleles at the site
	std::vector<int> a; //PBWT permutation
	std::vector<int> br;
	std::vector<std::pair<int, int> > FindBlocks();
	Cursor(){
		N = 0;
		bufP=0;
		lf = 100;
		numBrRec = 0;
		numSiteRec = 0;
		treeSize = 0;
//		maxBrRec = 0;
	}
} cur;

struct Branch{
	int start; //Beginning ob the block of ones
	int length; //Length of block of ones
	int parent; //Parent branch, 0 if top branch
	int tmp;
	int fSeen, lSeen;
	std::set<int> ch;
	std::list<std::pair<int,int> > unaf; //Positions of unaffected (by recombination) individuals relatively to branch
	char mode; //'c' for normal coalescent branch, 'r' for recombined branch with unknown position, 's' for suspected as recombined
	Branch(){
		mode = 'c';
	}
};

class Tree{
	public:
		std::vector<Branch> branch; //Tree structure
		int size; //Number of branches in the tree
		std::vector<std::vector<int> > recombs;
		std::vector<std::pair<int, std::vector<std::pair<int, int> > > > cr;
		//tree[0] is used to keep top branches
		Tree(){
			Branch nb;
			branch.push_back(nb);
			size = 1;
		}
		int lab;
		std::vector<int> delBrs; //indices of logically deleted branches
		void InitTree(Cursor* );
		std::pair<int,int> LevelTest(std::vector<std::pair<int, int> > , int , Cursor* );
		void AddNewBranch(int , int , int , std::vector<int> , Cursor* );
		int UpdateTree(Cursor*, std::vector<std::pair<int, int> > );
		void ParseRecomb(int , std::vector<std::pair<int, int> > );
		void PrintTree(bool , bool);
		int PositionRecomb(std::vector<int> , Cursor* );
		int ApplyPBWT(Cursor* , int , bool );
		void AddRecomb(int , std::vector<int>, Cursor* , bool );
		int ApplyRecomb(std::vector<int> , Cursor* );
		void NewickTree(Cursor*);
		void Distance(double , std::vector<double>&, std::vector<double>&, Cursor* );
		void BranchDist(std::vector< std::vector<unsigned long long> >& , std::vector<int>& , Cursor* );
		void BranchDistPop(std::vector< std::vector<unsigned long long> >& , std::vector<int>& , std::vector<std::pair<int, int> > , int , Cursor* );
		void TreeSize(Cursor* );
		int CommonAncestor(int , int , int*);
		void BranchDistPop1(std::vector< std::vector<unsigned long long> >& , std::vector<int>& , std::vector<std::pair<int, int> > , int , Cursor* );
		
		void NewickTreeBuf(Cursor*);
} tree;

void Tree::ParseRecomb(int parent, std::vector<std::pair<int, int> > blocks){
	int i,j,k=0;
	int startBranch, endBranch;
	int startBlock, endBlock;
	std::pair<int, std::vector<std::pair<int, int> > > nbr;
	std::set<int>::iterator it;
	char c;

	Tree::cr.pop_back();
	for (it = tree.branch[parent].ch.begin(); it != tree.branch[parent].ch.end(); it++){
		nbr.first = *it;
		startBranch = tree.branch[*it].start;
		endBranch = tree.branch[*it].start + tree.branch[*it].length - 1;
		for (i = 0; i < blocks.size(); i++){
			startBlock = blocks[i].first;
			endBlock = blocks[i].first + blocks[i].second - 1;
			if (endBlock < startBranch || endBranch < startBlock )
				continue;
			if (startBranch >= startBlock && endBranch <= endBlock)
				continue;
			if (startBranch <= startBlock && endBranch >= endBlock){
				nbr.second.push_back(blocks[i]);
				blocks.erase(blocks.begin()+i);
				i--;
				continue;
			}
			if (startBranch < startBlock){
				nbr.second.push_back(std::make_pair(blocks[i].first, tree.branch[*it].start + tree.branch[*it].length - blocks[i].first));
				blocks[i].first = tree.branch[*it].start + tree.branch[*it].length;
				blocks[i].second = endBlock - blocks[i].first + 1;
			}
			if (startBranch > startBlock){
				nbr.second.push_back(std::make_pair(tree.branch[*it].start, blocks[i].first + blocks[i].second - tree.branch[*it].start));
				blocks[i].second = tree.branch[*it].start - blocks[i].first;
			}
		}
		if (nbr.second.size() > 0)
			Tree::cr.push_back(nbr);
		nbr.second.clear();
		k++;
	}
	if (blocks.size() > 0){
		nbr.first = parent;
		nbr.second.swap(blocks);
		Tree::cr.push_back(nbr);
	}
}

int Tree::ApplyRecomb(std::vector<int> cr, Cursor* cur){
	int i, ml = 0, mb;
	std::vector<int> y;
	cur->numBrRec += cr.size();
	cur->numSiteRec++;
	cur->countBrRec[cr.size()]++;
	for (i = 0; i < cr.size(); i++)
		if (Tree::branch[cr[i]].length > ml){
			mb = i;
			ml = Tree::branch[cr[i]].length;
		}
	Tree::AddRecomb(mb, cr, cur, false);
	return 1;
}

int Tree::PositionRecomb(std::vector<int> nr, Cursor* cur){
	int i, j, k, pos, p, p1, posF, sb, mode;
	int u1 = 0, u2 = 0, v1 = 0, v2 = 0, w = 0;
	bool f = false;
	for (i = 0; i < nr.size(); i++){
		for (j = 0; j < Tree::recombs.size(); j++){
			for (k = 0; k < Tree::recombs[j].size(); k++){
				if ( nr[i] != Tree::recombs[j][k] && !(Tree::branch[nr[i]].start > Tree::branch[Tree::recombs[j][k]].start + Tree::branch[Tree::recombs[j][k]].length - 1 || Tree::branch[Tree::recombs[j][k]].start > Tree::branch[nr[i]].start + Tree::branch[nr[i]].length - 1)){
					pos = Tree::branch[nr[i]].start;
					p = Tree::branch[nr[i]].parent;
					f = true;
					break;
				}
			}
			if (f)
				break;
		}
		if (f)
			break;
	}
	if (!f)
		return 0;
	sb = nr[i];
	Tree::AddRecomb(k, Tree::recombs[j], cur, false);
	Tree::recombs.erase(Tree::recombs.begin()+j);
	for (i = 0; i < cur->M; i++)
		cur->y[i] = 0;
	for (i = 0; i < nr.size(); i++)
		for (j = Tree::branch[nr[i]].start; j < Tree::branch[nr[i]].start + Tree::branch[nr[i]].length; j++)
			cur->y[j] = 1;
	mode = Tree::UpdateTree(cur, cur->FindBlocks());
	if (mode == 0){
		i = 0;
		while (Tree::recombs.back()[i] != sb)
			i++;
		Tree::AddRecomb(i, Tree::recombs.back(), cur, false);
		Tree::recombs.erase(Tree::recombs.end());
	}
	return 1;
}

void Tree::AddRecomb(int sb, std::vector<int> nr, Cursor* cur, bool unafF = false){
	int i, j, k, pos, p, p1, posF, ch1, p2 = -1, par, newpar;
	int u1 = 0, u2 = 0, v1 = 0, v2 = 0, w = 0, is;
	int startBranch, endBranch;
	bool f = false, f1 = true, f2;
	std::vector<int> y;
	std::vector<int> b1, b2, c1, c2, perm;
	std::vector<Branch>::iterator it;
	std::set<int>::iterator jt;
	std::list<std::pair<int,int> >::iterator it1;
	pos = Tree::branch[nr[sb]].start;
	char symb;
	
	if (cur->N >= -1 && cur->N <= -10000){
		Tree::PrintTree(true, false);
		std::cout << "Recombination:";
		for (i = 0; i < nr.size(); i++)
			std::cout << " " << nr[i];
		std::cout << std::endl << "stable branch " << nr[sb] << std::endl;
	}

	if (unafF){//Update unaffected guys
		par = 0;
		f2 = true;
		while (f2){
			for (jt = tree.branch[par].ch.begin(); jt != tree.branch[par].ch.end(); jt++){
				startBranch = tree.branch[*jt].start;
				endBranch = tree.branch[*jt].start + tree.branch[*jt].length - 1;
				if (Tree::branch[nr[0]].start >= startBranch && Tree::branch[nr[0]].start + Tree::branch[nr[0]].length - 1 <= endBranch){
					newpar = *jt;
					break;
				}
			}
			for (i = 1; i < nr.size(); i++)
				if (!(Tree::branch[nr[i]].start >= startBranch && Tree::branch[nr[i]].start + Tree::branch[nr[i]].length - 1 <= endBranch)){
					f2 = false;
					break;
				}
			if (f2)
				par = newpar;
		}
		for (i = 0; i < nr.size(); i++){
			startBranch = Tree::branch[nr[i]].start;
			endBranch = Tree::branch[nr[i]].start + Tree::branch[nr[i]].length - 1;
			p1 = Tree::branch[nr[i]].parent;
			while (p1 != par){
				for (it1 = tree.branch[par].unaf.begin(); it1 != tree.branch[par].unaf.end(); it1++){
					if (it1->first + it1->second - 1 < startBranch)
						continue;
					if (it1->first > endBranch)
						break;
					if (it1->first >= startBranch && it1->first + it1->second - 1 <= endBranch)
						tree.branch[par].unaf.erase(it1);
					else if (it1->first >= startBranch){// && it1->first + it1->second - 1 > endBranch)
						it1->second = it1->first + it1->second - Tree::branch[p1].start - Tree::branch[p1].length;
						it1->first = Tree::branch[p1].start + Tree::branch[p1].length;
					}
					else if (it1->first < startBranch)// && it1->first + it1->second - 1 <= endBranch)
						it1->second = startBranch - it1->first;
					else{
						tree.branch[par].unaf.insert(it1, std::make_pair(it1->first, startBranch - it1->first));
						it1->second = it1->first + it1->second - Tree::branch[p1].start - Tree::branch[p1].length;
						it1->first = Tree::branch[p1].start + Tree::branch[p1].length;
					}
				}
				p1 = Tree::branch[p1].parent;
			}
		}
	}

	for (i = 0; i < cur->M; i++)
		y.push_back(0);
	for (i = 0; i < nr.size(); i++){
		for (j = Tree::branch[nr[i]].start; j < Tree::branch[nr[i]].start+Tree::branch[nr[i]].length; j++){
			y[j] = 1;
		}
	}
	for (i = 0; i < Tree::branch.size(); i++)
		Tree::branch[i].tmp = 0;
	for (i = 0; i < nr.size(); i++){
		if (i == sb)
			continue;
		p1 = Tree::branch[nr[i]].parent;
		Tree::branch[p1].ch.erase(nr[i]);

		while (p1 != 0){
			Tree::branch[p1].length = Tree::branch[p1].length - Tree::branch[nr[i]].length;
			if (Tree::branch[p1].start == Tree::branch[nr[i]].start)
				Tree::branch[p1].start = Tree::branch[p1].start + Tree::branch[nr[i]].length;
			if (f1){
				if (Tree::branch[p1].ch.size() == 1){
					ch1 = *Tree::branch[p1].ch.begin();
					if (Tree::branch[p1].start == Tree::branch[ch1].start && Tree::branch[p1].length == Tree::branch[ch1].length){
						Tree::branch[ch1].parent = Tree::branch[p1].parent;
						Tree::branch[p1].length = 0;
						Tree::branch[p1].ch.clear();
						Tree::branch[Tree::branch[p1].parent].ch.insert(ch1);
						Tree::branch[Tree::branch[p1].parent].ch.erase(p1);
						Tree::delBrs.push_back(p1);
					}
				}
				f1 = false;
			}
			p1 = Tree::branch[p1].parent;
		}
		f1 = true;
		p1 = Tree::branch[nr[sb]].parent;
		while(p1 != 0){
			Tree::branch[p1].length = Tree::branch[p1].length + Tree::branch[nr[i]].length;
			if (Tree::branch[p1].start > Tree::branch[nr[i]].start)
				Tree::branch[p1].tmp += Tree::branch[nr[i]].length;
			Tree::branch[p1].mode = 'r';
			p2 = p1;
			p1 = Tree::branch[p1].parent;
		}
		if (p2 > 0){
			if (Tree::branch[p1].start == Tree::branch[p2].start && Tree::branch[p1].length == Tree::branch[p2].length){
				Tree::branch[p1].ch.erase(p2);
				for (jt = Tree::branch[p2].ch.begin(); jt != Tree::branch[p2].ch.end(); jt++){
					Tree::branch[p1].ch.insert(*jt);
					Tree::branch[*jt].parent = p1;
				}
				Tree::branch[p2].ch.clear();
				Tree::branch[p2].length = 0;
				Tree::delBrs.push_back(p2);
			}
		}
		Tree::branch[nr[i]].parent = nr[sb];
		Tree::branch[nr[sb]].length += Tree::branch[nr[i]].length;
		Tree::branch[nr[sb]].ch.insert(nr[i]);
	}
	for (i = 0; i < Tree::branch.size(); i++)
		Tree::branch[i].start -= Tree::branch[i].tmp;
	for (i = 0; i < pos; i++){
		if (y[i] == 1)
			w++;
	}
	posF = pos - w;
	Tree::branch[nr[sb]].start = posF;
	Tree::branch[nr[sb]].mode = 'r';
	perm.reserve(cur->M);
	for (i = 0; i < pos; i++){
		if (y[i] == 1){
			b1.push_back(cur->a[i]);
			perm[i] = u1++;
		}
		else{
			c1.push_back(cur->a[i]);
			perm[i] = v1++;
		}
	}
	for (i = pos; i < cur->M; i++){
		if (y[i] == 1){
			b2.push_back(cur->a[i]);
			perm[i] = u2++;
		}
		else{
			c2.push_back(cur->a[i]);
			perm[i] = v2++;
		}
	}
	for (i = 0; i < v1; i++)
		cur->a[i] = c1[i];
	for (i = v1; i < v1 + u1; i++)
		cur->a[i] = b1[i - v1];
	for (i = v1 + u1; i < v1 + u1 + u2; i++)
		cur->a[i] = b2[i - v1 - u1];
	for (i = v1 + u1 + u2; i < v1 + u1 + u2 + v2; i++)
		cur->a[i] = c2[i - v1 - u1 - u2];

	for ( it = Tree::branch.begin()+1; it != Tree::branch.end(); it++){//Update branch positions
		if (it->mode == 'c'){
			if (it->start < pos){
				if (y[it->start] == 1)
					it->start = v1 + perm[it->start];
				else
					it->start = perm[it->start];
			}
			else{
				if (y[it->start] == 1)
						it->start = u1 + v1 + perm[it->start];
				else
					it->start = u1 + v1 + u2 + perm[it->start];
			}
		}
		else
			it->mode = 'c';
	}
}

void Tree::InitTree(Cursor* cur){
	Tree::branch[0].start = 0;
	Tree::branch[0].length = cur->M;
	Tree::branch[0].parent = 0;
}

std::pair<int,int> Tree::LevelTest(std::vector<std::pair<int, int> > blocks, int parent, Cursor* cur){ //parent = 0 if we are at the top, return -1 for recombination, -2 for placing the branch at the given level
//TODO function does not capture identical columns
	int i, j;
	std::set<int>::iterator it;
	int newParent = -1, st, len, f=0;
	char mode = 'c';//"c" for coalescent (can change to "p" if a new parent is found), "r" for recombination,
	cur->br.clear();
	if (blocks.size() == 0)
		return std::make_pair(-3,0);
	if (blocks.size() == 1){
		if (blocks[0].first == Tree::branch[parent].start && blocks[0].second == Tree::branch[parent].length)
			return std::make_pair(-4, parent);
		for (it = Tree::branch[parent].ch.begin(); it != Tree::branch[parent].ch.end(); it++){
			st = Tree::branch[*it].start;
			len = Tree::branch[*it].length;
			if (blocks[0].first == st && blocks[0].first + blocks[0].second - 1 == st + len - 1)
				return std::make_pair(-4,*it);
		}
	}
		
	for (it = Tree::branch[parent].ch.begin(); it != Tree::branch[parent].ch.end(); it++){//Here we check if we need to go to another level
		st = Tree::branch[*it].start;
		len = Tree::branch[*it].length;
			if ((blocks[0].first >= st && blocks[0].first + blocks[0].second - 1 < st + len - 1) || (blocks[0].first > st && blocks[0].first + blocks[0].second - 1 <= st + len - 1)){
			newParent = *it;
			mode = 'p';
			for (i = 0; i < blocks.size(); i++){
				if (blocks[i].first < st || blocks[i].first + blocks[i].second - 1 > st + len - 1){
					mode = 'r';
					break;
				}
			}
		}
		if (mode == 'p' || mode == 'r')
			break;
	}
	if (mode == 'c'){//If we do not go to the next level, we look for children   TODO optimize for(){for()...}
		for (i = 0; i < blocks.size(); i++){
			for (it = Tree::branch[parent].ch.begin(); it != Tree::branch[parent].ch.end(); it++){
				st = Tree::branch[*it].start;
				len = Tree::branch[*it].length;
				if (blocks[i].first <= st && blocks[i].first + blocks[i].second - 1 >= st + len - 1){
					cur->br.push_back (*it);
				}
				else if (!(blocks[i].first + blocks[i].second - 1 < st || blocks[i].first > st + len - 1)){
					mode = 'r';
					break;
				}
			}
			if (mode == 'r')
				break;
		}
	}
	if (mode == 'p')
		return std::make_pair(0,newParent);
	else if (mode == 'r')
		return std::make_pair(-1,0);
	else
		return std::make_pair(-2,parent);
}

void Tree::AddNewBranch(int start, int len, int parent, std::vector<int> ch, Cursor* cur){
	std::vector<int>::iterator it;
	Branch nb;
	int ind;
	if (Tree::delBrs.size() > 0){
		ind = Tree::delBrs.back();
		Tree::branch[ind].start = start;
		Tree::branch[ind].length = len;
		Tree::branch[ind].parent = parent;
		Tree::branch[ind].fSeen = cur->N;
		Tree::branch[ind].lSeen = cur->N;
		Tree::branch[ind].unaf.clear();
		Tree::branch[ind].unaf.push_back(std::make_pair(start, len));
		Tree::branch[ind].ch.clear();
		for ( it = ch.begin(); it != ch.end(); it++){
			Tree::branch[ind].ch.insert(*it);
			Tree::branch[parent].ch.erase(*it);
			Tree::branch[*it].parent = ind;
		}
		Tree::branch[parent].ch.insert(ind);
		Tree::delBrs.pop_back();
		Tree::lab = ind;
	}
	else{
		nb.start = start;
		nb.length = len;
		nb.parent = parent;
		nb.fSeen = cur->N;
		nb.lSeen = cur->N;
		nb.unaf.push_back(std::make_pair(start, len));
		for ( it = ch.begin(); it != ch.end(); it++){
			nb.ch.insert(*it);
			Tree::branch[parent].ch.erase(*it);
			Tree::branch[*it].parent = Tree::branch.size();
		}
		Tree::branch.push_back(nb);
		Tree::branch[parent].ch.insert(Tree::branch.size()-1);
		Tree::lab = Tree::branch.size()-1;
		Tree::size++;
	}
}

std::vector<std::pair<int, int> > Cursor::FindBlocks(){
	int i, t1;
	std::vector<std::pair<int, int> > blocks;
	
	for (i = 0; i < Cursor::M; i++){//Find all blocks of ones in the column, give their start positions and lenghts
		if (Cursor::y[i] == 1){
			t1 = i;
			while(Cursor::y[i] == 1 && i < Cursor::M)
				i++;
			blocks.push_back(std::make_pair(t1, i-t1));
		}
	}
	return blocks;
}

int Tree::ApplyPBWT(Cursor* cur, int p = 0, bool updateTree = false){
		int posF = 0, u = 0, v = 0;
		int i,j;
		std::vector<int> b, c, perm;
		std::vector<Branch>::iterator it;

		while (cur->y[posF] == 0)//Apply PBWT
			posF++;
		perm.reserve(cur->M-posF);
		for (i = posF; i < cur->M; i++){
			if (cur->y[i] == 1){
				b.push_back(cur->a[i]);
				perm[i-posF] = u++;
			}
			else{
				c.push_back(cur->a[i]);
				perm[i-posF] = v++;
			}
		}
		//TODO faster operation, analogue of memcpy?
		for (i = 0; i < u; i++)
			cur->a[posF+i] = b[i];
		for (i = 0; i < v; i++)
			cur->a[posF+u+i] = c[i];
		if (updateTree){
			for ( it = Tree::branch.begin()+1; it != Tree::branch.end(); it++){//Update branch positions
				if (it->start >= posF){
					if (cur->y[it->start] == 1){
						it->start = posF + perm[it->start-posF];
					}
					else{
						it->start = posF + u + perm[it->start-posF];
					}
				}
			}
			for (i = 0; i < Tree::cr.size(); i++){
				for (j = 0; j < Tree::cr[i].second.size(); j++)
					if (Tree::cr[i].second[j].first >= posF){
						if (cur->y[Tree::cr[i].second[j].first] == 1){
							Tree::cr[i].second[j].first = posF + perm[Tree::cr[i].second[j].first-posF];
						}
						else{
							Tree::cr[i].second[j].first = posF + u + perm[Tree::cr[i].second[j].first-posF];
						}
					}
			}
			Tree::AddNewBranch(posF, b.size(), p, cur->br, cur);
		}
		return 0;
}

int Tree::UpdateTree(Cursor* cur, std::vector<std::pair<int, int> > blocks){//Returns -1 for a new branch, 0 for a new recombination, -2 otherwise
	int i, j, posF = 0;
	int t1, t2;
	int u = 0, v = 0;
	int p = 0, tmp1 = 0;
	char symb;
	std::pair<int, int> f = std::make_pair(0,0);
	std::vector<int> b, c, perm;
	std::vector<Branch>::iterator it;
	std::pair<int, std::vector<std::pair<int, int> > > tb;
	std::vector<int> nr;
	
	while (f.first == 0)
		f = Tree::LevelTest(blocks, f.second, cur);
	p = f.second;
	if (f.first == -2){
		Tree::ApplyPBWT(cur, p, true);
		return -1;
	}
	if (f.first == -4){
		Tree::branch[p].lSeen = cur->N;
	}
	if (f.first == -1){
		Tree::cr.push_back(std::make_pair(p, blocks));
		Tree::recombs.reserve(Tree::recombs.size()+1);
		while (Tree::cr.size() > 0){
			tmp1++;
			if (Tree::cr.back().second.size() == 0)
				Tree::cr.pop_back();
			tb = Tree::cr.back();
			f = Tree::LevelTest(tb.second, tb.first, cur);
			switch(f.first){
				case -1:
					Tree::ParseRecomb(tb.first, tb.second);
					break;
				case -2:
					for (i = 0; i < cur->M; i++){
						cur->y[i] = 0;
					}
					for (j = 0; j < tb.second.size(); j++){
						for (i = tb.second[j].first; i < tb.second[j].first + tb.second[j].second; i++)
							cur->y[i] = 1;
					}
					Tree::ApplyPBWT(cur, f.second, true);
					nr.push_back(Tree::lab);
					Tree::cr.pop_back();
					break;
				case -3:
					std::cout << "Tree::UpdateTree: Unexpected LevelTest returned value.\n" << std::endl;
					exit (EXIT_FAILURE);
					break;
				case -4:
					nr.push_back(f.second);
					Tree::cr.pop_back();
					break;
				case 0:
					Tree::cr[Tree::cr.size()-1].first = f.second;
					break;
				default:
					std::cout << "Tree::UpdateTree: Unknown LevelTest exeption.\n" << std::endl;
					exit (EXIT_FAILURE);		
					break;
			}
		}
		if (!Tree::ApplyRecomb(nr, cur))
			Tree::recombs.push_back(nr);
		return 0;
	}
	return -2;
}

void Tree::PrintTree(bool tr = true, bool rec = true){
	int i, j, k =0;
	std::vector<Branch>::iterator it, kt;
	std::set<int>::iterator jt;
	if (!outputTree)
		return;
	i = 0;
	for (it = Tree::branch.begin(); it != Tree::branch.end(); it++){
		for (kt = Tree::branch.begin(); kt != Tree::branch.end(); kt++){
			k++;
			if (i == k-1)
				continue;
			if (kt->start > it->start + it->length - 1)
				continue;
			if (it->start > kt->start + kt->length - 1)
				continue;
			if (it->start >= kt->start && it->start + it->length - 1 <= kt->start + kt->length - 1)
				continue;
			if (it->start <= kt->start && it->start + it->length - 1 >= kt->start + kt->length - 1)
				continue;
			std::cout << "INVALID TREE: branches " << i << " " << k-1 << std::endl;
		}
		i++;
		k=0;
	}
	i = 0;
	if (tr){
		std::cout << "\t\tTREE" << std::endl;
		for (it = Tree::branch.begin(); it != Tree::branch.end(); it++){
//			if (it->length == 0 || it->length == 1){
//				i++;
//				continue;
//			}
			std::cout << "\t\t\tBranch " << i << ", parent " << it->parent << ", begin at " << it->start << ", length " << it->length << std::endl;
			std::cout << "\t\t\t\tChildren:";
			for (jt = it->ch.begin(); jt != it->ch.end(); jt++)
				if (Tree::branch[*jt].length != 0 && Tree::branch[*jt].length != 1)
					std::cout << *jt << " ";
			std::cout << std::endl;
			i++;
		}
	}
	if (rec){
		std::cout << "\t\tRECOMBINATIONS" << std::endl;
		for (i = 0; i < recombs.size(); i++){
			std::cout << "\t\t\tRecombination " << i << ":";
			for (j = 0; j < recombs[i].size(); j++){
				std::cout << " " << recombs[i][j];
			}
			std::cout << std::endl;
		}
		std::cout << "\t\t--------------" << std::endl;
	}
}

void Tree::Distance(double s, std::vector<double>& dist, std::vector<double>& cd, Cursor* cur){
	int i, j;
	int alpha = 0;
	std::vector<Branch>::iterator it;
	char c;
//	if (cur->N%1000 == 0)
//	std::cout << "Tree size " << Tree::branch.size() << std::endl;		
//	std::cout << "Vector size " << dist.size() << std::endl;
//	std::cin >> c;
	for (i = 0; i < cur->M*cur->M; i++)
		cd[i] = 0;
	if (s == 1.0)
		s = distParameter;
	for (it = Tree::branch.begin(); it != Tree::branch.end(); it++){
		if (it->length < 2 || it->start < 0)
			continue;
		for (i = it->start; i < it->start+it->length-1; i++){
			for (j = i+1; j < it->start+it->length; j++){
				if (cd[cur->M*cur->a[i]+cur->a[j]] < exp(-s*(it->length-minBr)+alpha*(2*s-1)) && it->length >= minBr && (maxBr == 0 || it->length <= maxBr) ){
					cd[cur->M*cur->a[i]+cur->a[j]] = exp(-s*(it->length-minBr)+alpha*(2*s-1));
					cd[cur->M*cur->a[j]+cur->a[i]] = exp(-s*(it->length-minBr)+alpha*(2*s-1));
				}
			}
		}
	}
	for (i = 0; i < cur->M*cur->M; i++)
		dist[i] += cd[i];
//	if (cur->N == 1000)
//		std::cout << dist[1] << std::endl;
}

/*
void Tree::BranchDist(std::vector< std::vector<int> >& brDist, std::vector<int>& cd, Cursor* cur){
	int i, j;
	int alpha = 0;
	std::vector<Branch>::iterator it;
	char c;
//	if (cur->N%1000 == 0)
//	std::cout << "Tree size " << Tree::branch.size() << std::endl;		
//	std::cout << "Vector size " << dist.size() << std::endl;
//	std::cin >> c;
	for (i = 0; i < cur->M*cur->M; i++)
		cd[i] = cur->M;
	for (it = Tree::branch.begin(); it != Tree::branch.end(); it++){
		if (it->length < 2 || it->start < 0)
			continue;
		for (i = it->start; i < it->start+it->length-1; i++){
			for (j = i+1; j < it->start+it->length; j++){
				if (cd[cur->M*cur->a[i]+cur->a[j]] > it->length){
					cd[cur->M*cur->a[i]+cur->a[j]] = it->length;
					cd[cur->M*cur->a[j]+cur->a[i]] = it->length;
				}
			}
		}
	}
	for (i = 0; i < cur->M*cur->M; i++)
		brDist[i][cd[i]-1] ++;
//	if (cur->N == 1000)
//		std::cout << dist[1] << std::endl;
}
*/

int Tree::CommonAncestor(int h1, int h2, int *iter){
	int brID = 0;
	int st, end;
	std::set<int>::iterator it;
	char c;
	bool flag = true;
	int control = 0;
//	if (cur->N%1000 == 0)
//	std::cout << "Tree size " << Tree::branch.size() << std::endl;		
//	std::cout << "Vector size " << dist.size() << std::endl;
//	std::cin >> c;
	if (Tree::branch[brID].ch.size() == 0)
		return brID;
	while(Tree::branch[brID].ch.size() > 0 && flag){
		control++;
		flag = false;
		for (it = Tree::branch[brID].ch.begin(); it != Tree::branch[brID].ch.end(); it++){
			if (Tree::branch[*it].length < 2)
				continue;
			*iter = *iter + 1;
			st = Tree::branch[*it].start;
			end = st+Tree::branch[*it].length;
			if (h1 >= st && h1 < end && h2 >= st && h2 < end){
				brID = *it;
				flag = true;
				break;
			}
		}
	}
	return brID;
}

void Tree::BranchDistPop1(std::vector< std::vector<unsigned long long> >& brDistPop, std::vector<int>& ar, std::vector<std::pair<int, int> > hapGroups, int numPop, Cursor* cur){
	int i, j, k, l;
	std::vector<Branch>::iterator it;
	char c;
	int ca, ind;
	int counter = 0, count2 = 0;
	int iter = 0;

	for (i = 0; i < cur->M; i++) //reversed permutation cur->a
		ar[cur->a[i]] = i;

	for (k = 0; k < numPop; k++){//Within population
		ind = (k+1)*(k+2)/2-1;
		for (i = hapGroups[k].first; i < hapGroups[k].second-1; i++){
			for (j = i+1; j < hapGroups[k].second; j++){
//				time(&t1);
				count2++;
				ca = Tree::CommonAncestor(ar[i], ar[j], &iter);
//				time(&t2);
//				t3 += difftime(t2,t1);
				brDistPop[ind][Tree::branch[ca].length-1]++;
			}
		}
	}
	for (k = 0; k < numPop-1; k++){//Between populations
		for (l = k+1; l < numPop; l++){
			ind = (l+1)*l/2+k;
			for (i = hapGroups[k].first; i < hapGroups[k].second; i++){
				for (j = hapGroups[l].first; j < hapGroups[l].second; j++){
//					time(&t1);
					count2++;
					ca = Tree::CommonAncestor(ar[i], ar[j], &iter);
//					time(&t2);
//					t3 += difftime(t2,t1);
					brDistPop[ind][Tree::branch[ca].length-1]++;
	//				counter+=4;
				}
			}
		}
	}
	Tree::TreeSize(cur);
	std::cout << "iter = " << iter << ", couples= " << count2 << ", average operations = " << (double)iter/(double)count2 << std::endl;
	std::cout << "tree size = " << cur->treeSize << std::endl;
}

void Tree::BranchDistPop(std::vector< std::vector<unsigned long long> >& brDistPop, std::vector<int>& cd, std::vector<std::pair<int, int> > hapGroups, int numPop, Cursor* cur){
	int i, j, k, l;
	std::vector<Branch>::iterator it;
	char c;
	int h1, h2, ind, ind1;
	int counter = 0;
//	if (cur->N%1000 == 0)
//	std::cout << "Tree size " << Tree::branch.size() << std::endl;		
//	std::cout << "Vector size " << dist.size() << std::endl;
//	std::cin >> c;
	for (i = 0; i < cur->M*(cur->M-1)/2; i++)
		cd[i] = cur->M;
	for (it = Tree::branch.begin()+1; it != Tree::branch.end(); it++){
		if (it->length < 2 || it->start < 0)
			continue;
		for (i = it->start; i < it->start+it->length-1; i++){
			for (j = i+1; j < it->start+it->length; j++){
				counter+=2;
				h1 = cur->a[i]<cur->a[j]?cur->a[i]:cur->a[j];
				h2 = cur->a[i]+cur->a[j]-h1;
				ind = (h2-1)*h2/2+h1;
				if (cd[ind] > it->length)
					cd[ind] = it->length;
			}
		}
	}
	for (k = 0; k < numPop; k++){//Within population
		ind1 = (k+1)*(k+2)/2-1;
		for (i = hapGroups[k].first; i < hapGroups[k].second; i++){
			for (j = i+1; j < hapGroups[k].second; j++){
				ind = j*(j-1)/2+i;
				brDistPop[ind1][cd[ind]-1]++;
	//			counter+=4;
			}
		}
	}
	for (k = 0; k < numPop-1; k++){//Between populations
		for (l = k+1; l < numPop; l++){
			ind1 = (l+1)*l/2+k;
			for (i = hapGroups[k].first; i < hapGroups[k].second; i++){
				for (j = hapGroups[l].first; j < hapGroups[l].second; j++){
					ind = j*(j-1)/2+i;
					brDistPop[ind1][cd[ind]-1]++;
	//				counter+=4;
				}
			}
		}
	}
	//std::cout << "Number of operations " << counter << std::endl;
//	if (cur->N == 1000)
//		std::cout << dist[1] << std::endl;
}


void Tree::BranchDist(std::vector< std::vector<unsigned long long> >& brDist, std::vector<int>& cd, Cursor* cur){
	int i, j, k, l;
	std::vector<Branch>::iterator it;
	int h1, h2;
	char c;
	int counter = 0;
//	if (cur->N%1000 == 0)
//	std::cout << "Tree size " << Tree::branch.size() << std::endl;		
//	std::cout << "Vector size " << dist.size() << std::endl;
//	std::cin >> c;
	for (i = 0; i < cur->M; i++)
		cd[i] = cur->M;
	for (it = Tree::branch.begin()+1; it != Tree::branch.end(); it++){
		if (it->length < 2 || it->start < 0)
			continue;
		for (i = it->start; i < it->start+it->length-1; i++){
			for (j = i+1; j < it->start+it->length; j++){
				h1 = cur->a[i];
				h2 = cur->a[j];
				if (cd[h1+cur->M*h2] > it->length){
					cd[h1+cur->M*h2] = it->length;
					cd[h2+cur->M*h1] = it->length;
				}
			}
		}
	}
	for (i = 0; i < cur->M*cur->M; i++){
		brDist[i][cd[i]-1]++;
	}
}


void Tree::TreeSize(Cursor* cur){
	std::vector<Branch>::iterator it;
//	if (cur->N < cur->M)
//		return;
	for (it = Tree::branch.begin(); it != Tree::branch.end(); it++){
		if (it->length > 1)
			cur->treeSize++;
	}
}

void ReadPureBinary(char* fname, Cursor* cur){
	time_t timer1, timer2;
	FILE *fp;
	char str[100000];
	char fname1[255], fname2[255];
	char dn[255];
	int i, j, k;
	int sum;
	int ltp;
	int countBrRec10 = 0;
	int tmp;
//	std::string str;
	std::vector<int> x;
	std::vector<int> ar;
	std::vector<double> dist, cd;
	std::vector< std::vector<unsigned long long> > brDist;
	std::vector<unsigned long long> init;
	int sampleNum = 0;
	double norm;
//	int popSize[] = {96, 61, 99, 113, 99, 85, 108, 94, 64, 85, 104, 93, 103, 105, 104, 99, 99, 99, 91, 107, 107, 86, 103, 102, 96, 102};
//	std::string popName[] = {"ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU"};
//	int numPop = 26;
	int popSize[] = {500};
	std::string popName[] = {"SIMMACS"};
	int numPop = 1;
	std::vector<std::pair<int, int> > hapGroups;
	std::vector<int> cd1;
//	std::istream fs (fname);
	char c, c1;
	tmp = 0;
	for (i = 0; i < numPop; i++){
		hapGroups.push_back(std::make_pair(tmp, tmp+2*popSize[i]));
		tmp+=(2*popSize[i]);
	}
	fp = fopen(fname, "r");
	if (!fp){
		std::cout << "Unable to open file " << fname << "." << std::endl;
		exit(EXIT_FAILURE);
	}
	time(&timer1);
	fgets(str, 100000, fp);
//	if (strncmp(str, "#BINARY", 7) != 0){
//		std::cout << "Corrupted binary file " << fname << ", header is corrupted." << std::endl;
//		exit(EXIT_FAILURE);
//	}
	cur->buf.resize(cur->lf);
	//Read first line of data, find cur->M, initialize tree
	do{
		c = fgetc(fp);
		if (c == '0' || c == '1'){
			if (c == '0')
				x.push_back(0);
			else
				x.push_back(1);
		}
		else if (c != '\n' && c != EOF){
			std::cout << "Fatal error: unexpected charachter \"" << c << "\"." << std::endl;
			exit(EXIT_FAILURE);
		}
	}while(c != '\n' && c!= EOF);
	if (x.size() < 1 || c == EOF){
		std::cout << "Too few of data." << std::endl;
		exit(EXIT_FAILURE);
	}
	cur->M = x.size();
	cur->buf[0].swap(x);
	tree.InitTree(cur);
	//Allocate memory, initialize vectors
	for (i = 0; i < cur->M; i++){
		cur->a.push_back(i);
		cur->y.push_back(0);
		cur->countBrRec.push_back(0);
	}
	if (relativeDist)
		for (i = 0; i < cur->M*cur->M; i++){
			dist.push_back(0);
			cd.push_back(0);
		}
	if (branchDist && branchDistPop){
		std::cout << "Simultaneous computation of haplotype to haplotype and population to population distances is not supported." << std::endl;
		exit(EXIT_FAILURE);
	}
	if (branchDist)
		for (i = 0; i < cur->M; i++)
			init.push_back(0);
		for (i = 0; i < cur->M*cur->M; i++){
			brDist.push_back(init);
			cd1.push_back(0);
		}
		init.clear();
	if (branchDistPop){
		for (i = 0; i < cur->M; i++){
			ar.push_back(0);
			init.push_back(0);
		}
		for (i = 0; i < cur->M*(cur->M-1)/2; i++)
			cd1.push_back(0);
		for (i = 0; i < numPop*(numPop+1)/2; i++)
			brDist.push_back(init);
	}
	x.clear();
	
	//Fill buffer with first lines
	ltp = 1;
	while (ltp < cur->lf && c != EOF){
		do{
			c = fgetc(fp);
			if (c == '0' || c == '1'){
				if (c == '0')
					x.push_back(0);
				else
					x.push_back(1);
			}
			else if (c != '\n' && c != EOF){
				std::cout << "Fatal error: unexpected charachter \"" << c << "\"." << std::endl;
				exit(EXIT_FAILURE);
			}
		}while(c != '\n' && c!= EOF);
		if (x.size() == 0)
			break;
		if (cur->M != x.size()){
			std::cout << "Fatal error: wrong number of haplotypes, line " << i << "." << std::endl;
			exit(EXIT_FAILURE);
		}
		cur->buf[ltp].swap(x);
		x.clear();
		ltp++;
	}
	//Start tcPBWT processing
	while(ltp > 0 && (cur->N < lineNum || lineNum == 0)  && (cur->N < numOfSites || numOfSites == -1)){
		for (i = 0; i < cur->M; i++)
			cur->y[i] = cur->buf[cur->bufP][cur->a[i]];
//		if (cur->N >= 0 && cur->N <= 10000000000)
//			tree.PrintTree(true, false);
		tree.NewickTree(cur);

		if ( (cur->N+1)%lineCounter == 0)
			std::cout << cur->N+1 << " lines analyzed." << std::endl;
		
//		if (cur->N > 3000000)
//			tree.NewickTreeBuf(cur);
		
		tree.UpdateTree(cur, cur->FindBlocks());
//		tree.TreeSize(cur);
		if (relativeDist)
			tree.Distance(1.0, dist, cd, cur);
		if (branchDist && cur->N%bDPsample==0 && cur->N > cur->M){
			tree.BranchDist(brDist, cd1, cur);
			sampleNum++;
		}
		if (branchDistPop && cur->N%bDPsample==0 && cur->N > cur->M){
			tree.BranchDistPop(brDist, cd1, hapGroups, numPop, cur);
			sampleNum++;
		}

//		for (i = 0; i < cur->lf; i++){
//			for (j = 0; j < cur->M; j++)
//				std::cout << cur->buf[(i+cur->bufP)%cur->lf][j];
//			std::cout << std::endl;
//		}
		if (c != EOF){
			x.clear();
			c = '2';
			do{
				c = fgetc(fp);
				switch(c){
					case '0':
						x.push_back(0);
						break;
					case '1':
						x.push_back(1);
						break;
					case '\n': case EOF:
						break;
					default:
						std::cout << "Fatal error: unexpected charachter \"" << c << "\"." << std::endl;
						exit(EXIT_FAILURE);
						break;
				}
			}while(c != '\n' && c!= EOF);
			if (x.size() == 0)
				break;
			if (cur->M != x.size()){
				std::cout << "Fatal error: wrong number of haplotypes, line " << cur->N+cur->lf << "." << std::endl;
				exit(EXIT_FAILURE);
			}
			cur->buf[cur->bufP].swap(x);
		}
		else
			ltp--;
		
		cur->N++;
		cur->bufP = (cur->bufP+1)%cur->lf;
	}
	time(&timer2);
	std::cout << "Binary file is read successfuly." << std::endl;
	std::cout << "M = " << cur->M << ", N = " << cur->N << ", average complexity of recombinations " << (double)cur->numBrRec/(double)cur->numSiteRec << ", number of sites with recombination " << cur->numSiteRec << std::endl;
	std::cout << "Total number of recombinations " << cur->numBrRec-cur->numSiteRec << ", average number of branches " << (double)cur->treeSize/(double)(cur->N-cur->M) << std::endl;
//	for (i = 0; i < cur->M; i++){
//		if (i > 9)
//			countBrRec10 += cur->countBrRec[i];
//		if (cur->countBrRec[i] > 0)
//			std::cout << i << "\t" << cur->countBrRec[i] << std::endl;
//	}
//	std::cout << ">9\t" << countBrRec10 << std::endl;
	std::cout << "Algorithm time " << difftime(timer2,timer1) << "s" << std::endl;
	fclose(fp);
	if (relativeDist){
		strcpy (fname1, fname);
		strcat (fname1, ".reldist");
		fp = fopen(fname1, "w");
		sum = 0;
		for (i = 0; i < cur->M*cur->M; i++){
//			dist[i] = log (dist[i]);
			sum += dist[i];
		}
		for (i = 0; i < cur->M*cur->M; i++){
			fprintf (fp, "%f", 2*cur->M*cur->M*dist[i]/sum);
	//		std::cout << std::fixed << 2*cur->M*cur->M*dist[i]/sum;
			if (i%cur->M == cur->M - 1)
				fprintf (fp, "\n");
	//			std::cout << std::endl;
			else
				fprintf (fp, "\t");
	//			std::cout << "\t";
		}
		fclose(fp);
	}
	if (relativeDist){
		strcpy (fname1, fname);
		strcat (fname1, ".reldist");
		fp = fopen(fname1, "w");
		sum = 0;
		for (i = 0; i < cur->M*cur->M; i++){
//			dist[i] = log (dist[i]);
			sum += dist[i];
		}
		for (i = 0; i < cur->M*cur->M; i++){
			fprintf (fp, "%f", 2*cur->M*cur->M*dist[i]/sum);
	//		std::cout << std::fixed << 2*cur->M*cur->M*dist[i]/sum;
			if (i%cur->M == cur->M - 1)
				fprintf (fp, "\n");
	//			std::cout << std::endl;
			else
				fprintf (fp, "\t");
	//			std::cout << "\t";
		}
		fclose(fp);
	}
	if (branchDistPop){
		strcpy (fname2, fname);
		strcat (fname2, ".branch.pop");
		fp = fopen(fname2, "w");
		fprintf (fp, "Within population\n");
		for (i = 0; i < numPop; i++){
			fprintf (fp, "%s\t", popName[i].c_str());
			tmp = popSize[i]*(2*popSize[i]-1);
			norm = 0;
			for (k = 0; k < cur->M; k++){
				fprintf (fp, "%f", (double)brDist[i*(i+3)/2][k]/(double)tmp);
				norm+=(double)brDist[i*(i+3)/2][k]/(double)tmp;
				if (k < cur->M-1)
					fprintf (fp, " ");
			}
			if (norm - sampleNum > 1 || norm - sampleNum < -1){
				printf("Norm is not equal to the number of sampling times.");
			}
			fprintf (fp, "\n");
		}

		fprintf (fp, "\n\nBetween populations\n");
		for (i = 0; i < numPop-1; i++){
			for (j = i+1; j < numPop; j++){
				fprintf (fp, "%s\t%s\t", popName[i].c_str(), popName[j].c_str());
				tmp = 4*popSize[i]*popSize[j];
				norm = 0;
				for (k = 0; k < cur->M; k++){
					fprintf (fp, "%f", (double)brDist[(j+1)*j/2+i][k]/(double)tmp);
					norm+=(double)brDist[i*(i+3)/2][k]/(double)tmp;
					if (k < cur->M-1)
						fprintf (fp, " ");
				}
				if (norm - sampleNum > 1 || norm - sampleNum < -1){
					printf("Norm is not equal to the number of sapling times.");
				}
				fprintf (fp, "\n");
			}
		}
		fclose(fp);
	}
	if (branchDist){
		fprintf (fp, "Within population\n");
		mkdir("PBWTout", 0775);
		for (i = 0; i < cur->M; i++){
			for (j = 0; j < cur->M; j++)
				brDist[i*(cur->M+1)][j] = 0;
		}
		for (i = 0; i < cur->M; i++){
			if (i%200 == 0){
				sprintf(str, "%d", i/200);
				strcpy(dn, "PBWTout/datap");
				strcat(dn, str);
				mkdir(dn, 0775);
			}
			sprintf(str, "%03d", i);
			strcpy(fname2,dn);
			strcat(fname2,"/denOneToOne_");
			strcat(fname2, str);
			strcat(fname2, ".csv");
			fp = fopen(fname2, "w");
			for (j = 0; j < cur->M; j++){
				fprintf (fp, "\"h%d_h%d\"", i, j);
				if (j < cur->M-1)
					fprintf (fp, ",");
				else
					fprintf (fp, "\n");
			}
			for (k = 1; k < cur->M; k++)
				for (j = 0; j < cur->M; j++){
					fprintf (fp, "%f", (double)brDist[j+i*cur->M][k]/(double)sampleNum);
					if (j < cur->M-1)
						fprintf (fp, ",");
					else
						fprintf (fp, "\n");
				}
			fclose(fp);
		}
	}
}

void Tree::NewickTree(Cursor *cur){
	if (!outputTree)
		return;
	
	std::vector<int> lbr, rbr;
	std::vector<Branch>::iterator it;
	int i, left, right;
	
	for (i = 0; i < cur->M; i++){
		lbr.push_back(0);
		rbr.push_back(0);
	}
	
	for (it = Tree::branch.begin(); it != Tree::branch.end(); it++){
		if (it->length == 0 || it->length == 1){
			i++;
			continue;
		}
		left = it->start;
		right = it->start + it->length - 1;
		lbr[left]++;
		rbr[right]++;
	}
	for (i = 0; i < cur->M; i++){
		while (lbr[i] > 0){
			std::cout << "(";
			lbr[i]--;
		}
		std::cout << cur->a[i];
		while (rbr[i] > 0){
			std::cout << ")";
			rbr[i]--;
		}
		if (i < cur->M - 1)
			std::cout << ",";
		else
			std::cout << "\n";
	}
}

void Tree::NewickTreeBuf(Cursor *cur){	
	std::vector<int> lbr, rbr, lbr1, rbr1;
	std::vector<Branch>::iterator it;
	int i, left, right;
	int j, count;
	char c1;
	
	for (i = 0; i < cur->M; i++){
		lbr.push_back(0);
		rbr.push_back(0);
	}
	
	for (it = Tree::branch.begin(); it != Tree::branch.end(); it++){
		if (it->length == 0 || it->length == 1){
			i++;
			continue;
		}
		left = it->start;
		right = it->start + it->length - 1;
		lbr[left]++;
		rbr[right]++;
	}
	rbr1 = rbr;
	lbr1 = lbr;
	for (i = 0; i < cur->M; i++)
		std::cout << cur->a[i] << ", ";
	std::cout << std::endl;
	for (j = 0; j < cur->lf; j++){
		rbr = rbr1;
		lbr = lbr1;
		count = 0;
		for (i = 0; i < cur->M && count < 2; i++)
			if (cur->buf[(j+cur->bufP)%cur->lf][i] == 1)
				count++;
		if (count < 2){
			if (j == 0)
				return;
			else
				continue;
		}
		for (i = 0; i < cur->M; i++){
			while (lbr[i] > 0){
				std::cout << "<";
				lbr[i]--;
			}
			if (cur->buf[(j+cur->bufP)%cur->lf][cur->a[i]] == 0)
				std::cout << ".";
			else
				std::cout << "#";
			while (rbr[i] > 0){
				std::cout << ">";
				rbr[i]--;
			}
			if (i < cur->M - 1)
				std::cout << "";
			else
				std::cout << "\n";
		}
	}
	std::cin >> c1;
}

void Help(){
	std::cout << "Help page for ARGentum." << std::endl;
	std::cout << "Compile simply with g++ -o argentum argentum.cpp" << std::endl;
	std::cout << "To get help run ./argentum without parameters." << std::endl;
	std::cout << "Run the program with the command ./argentum -r [input file name] ...optional flags." << std::endl;
	std::cout << "-not - cancel the output of local trees. By default local trees are outputed in Newick format to the stadndard output." << std::endl;
	std::cout << "-dist - calculate the relative distance between haplotypes based on a sum of exponents of minimal branch size distribution. You can also set the metric parameter with -param [FLOAT]. -minBr [INTEGER] and -maxBr [INTEGER] allows to set tresholds for minimal branch size wich are included in the summary." << std::endl;
	std::cout << "-brDist - calculate the pairwise minimal branch size distribution. Output in PBWTout folder. Ignores first M sites." << std::endl;
	std::cout << "-brDistPop - calculate the average of -brDist within and between populations. Population parameters should be pregiven in the code... TODO" << std::endl;
	std::cout << "-bDPsample [INTEGER] - allows to sample -brDist or -brDistPop, chooses every 1 in INTEGER trees." << std::endl;
	std::cout << "-lc [INTEGER] - set the line counter (reports every INTEGER line) to track the progress." << std::endl;
	std::cout << "-lf [INTEGER] - set the size of buffer (for example for look forward purpose)." << std::endl;
	std::cout << "-ll [INTEGER] - set the maxmimal number of sites to process."  << std::endl;
	std::cout << "-ns [INTEGER] - seems to be the same as -ll... TODO" << std::endl;
}

int main(int argc, char *argv[]){
	int i;
	char path[255];
	int lf;
//	strcpy (path,"data/sim_data_");
//	char* index = getenv ("LSB_JOBINDEX");
//	strcat (path, index);
//	char fname[]="../PBWT/data/pbwt2bin.txt";
//	char fname[]="PBWT-Data/20.macs";
//	char fname[]="30-notrees.txt";
	//char fname[]="PBWT-Data/1.1k.macs";
	std::cout << "Welcome to tcPBWT ARG." << std::endl;
//	FirstTest();
//	ReadMacs(fname, &cur);
	if (argc < 3){
		Help();
		return 0;
	}
	if (strcmp(argv[1], "-r") == 0){
		if (strlen(argv[2]) <= 255)
			strcpy (path, argv[2]);
		else{
			std::cout << "Too long filename, 255 symbols maximum." << std::endl;
			return 0;
		}
	}
	else{
		std::cout << "Type -r [filename] for input first." << std::endl;
		return 0;
	}
	i = 3;
	while (i < argc){
		if (strcmp(argv[i], "-not") == 0){
			outputTree = false;
			i++;
		}
		else if (strcmp(argv[i], "-dist") == 0){
			relativeDist = true;
			i++;
		}
		else if (strcmp(argv[i], "-param") == 0){
			distParameter = strtol(argv[i+1], NULL, 10);
			i+=2;
			if (!relativeDist){
				std::cout << "WARNING: You have set the flag -param, but you have not set the flag -dist, distance will not be computed." << std::endl;
			}
		}
		else if (strcmp(argv[i], "-min") == 0){
			minBr = strtol(argv[i+1], NULL, 10);
			i+=2;
			if (!relativeDist){
				std::cout << "WARNING: You have set the flag -min, but you have not set the flag -dist, distance will not be computed." << std::endl;
			}
		}
		else if (strcmp(argv[i], "-max") == 0){
			maxBr = strtol(argv[i+1], NULL, 10);
			i+=2;
			if (!relativeDist){
				std::cout << "WARNING: You have set the flag -max, but you have not set the flag -dist, distance will not be computed." << std::endl;
			}
		}
		else if (strcmp(argv[i], "-brDist") == 0){
			branchDist = true;
			i++;
		}
		else if (strcmp(argv[i], "-brDistPop") == 0){
			branchDistPop = true;
			i++;
		}
		else if (strcmp(argv[i], "-bDPsample") == 0){
			bDPsample = strtol(argv[i+1], NULL, 10);
			i+=2;
			if (!branchDistPop){
				std::cout << "WARNING: You have set the flag -bDPsample, but you have not set -brDist or -brDistPop flag, distance will not be computed." << std::endl;
			}
		}
		else if (strcmp(argv[i], "-lc") == 0){
			lineCounter = strtol(argv[i+1], NULL, 10);
			i+=2;
		}
		else if (strcmp(argv[i], "-lf") == 0) { //look forward
			cur.lf = strtol(argv[i+1], NULL, 10);
			i+=2;
		}
		else if (strcmp(argv[i], "-ll") == 0) { //line limit
			lineNum = strtol(argv[i+1], NULL, 10);
			i+=2;
		}
		else if (strcmp(argv[i], "-ns") == 0) { //number of sites to process
			numOfSites = strtol(argv[i+1], NULL, 10);
			i+=2;
		}
		else{
			std::cout << "Unknown flag." << std::endl;
			return 0;
		}
	}
	std::cout << "Input file: " << path << ".\n";
	ReadPureBinary(path, &cur);
//	tree.PrintTree(true, false);
//	for (i = 0; i < cur.M; i++)
//		std::cout << cur.a[i] << " ";
//	std::cout << std::endl;
	std::cout << "M = " << cur.M << ", N = " << cur.N << std::endl;
	return 0;
}

/*
branchDist (brDist) - one-to-one distribution of minimal branch sizes
branchDistPop - mean of the branchDist over populations (both within and between populations)
*/

/*Input function: get new column, put it in $cur->y[i]$ in the order given by the permutation $cur->a[i]$.
Call function Tree::UpdateTree(Cursor, Cursor::FindBlocks())
Run Tree::LevelTest until one of events happens:
- recombination detected
- branch already exists
- PBWT applied, new branch added


*/
