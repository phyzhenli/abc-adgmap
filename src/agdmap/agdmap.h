/**Header-File***************************************tab=2**************

  FileName    [agdmap.h]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Agdmap.]

  Synopsis    [External declarations.]

  Author      [Longfei Fan, Zhiyong Zhang, Chang Wu]

  Affiliation [Fudan University]

  Date        [August 1, 2023]

  Vision      [ 1.0 ]

***********************************************************************/

#ifndef _AGDMAP_H_
#define _AGDMAP_H_
#pragma once

#include "base/main/main.h" // the basic header file of abc
#include "utility.h"

namespace AgdMap {

class Node;
class Cut;
class Para;
class Solution;
class SimpleGate;
class Decompose;

// DO NOT FIX THESE TYPEDEF !
typedef float                   Level   ;
typedef float                   Depth   ;
typedef float                   Area    ;
typedef unsigned long           Sign    ;    
typedef Abc_Obj_t               AbcNode ;
typedef std::shared_ptr<Cut>    pCut    ;
typedef std::vector<pCut>       CutList ;
typedef std::vector<Node*>      NodeVec ;
typedef CutList::iterator       CutItr  ;
typedef utility::truth_table      Function;
typedef std::shared_ptr<SimpleGate>  pGate;
typedef std::shared_ptr<Decompose>   pDec;

/*=================== MACRO ===================*/
#define App   std::ios::out | std::ios::app
#define Trunc std::ios::out | std::ios::trunc
#define SSINF 1000000.0
#define Encode(node_id) ( 1UL << node_id % 63 )


//===----------------------------------------------------------------------===//
//                               Parameter class
//===----------------------------------------------------------------------===//
class Para {
public:
  Para(int k, int g, bool area, bool v) : lut_size_(k), gate_size_(g), area_oriented_(area), verbose_(v) {}
	Para& operator=(const Para& rhs) {
		lut_size_ = rhs.lut_size_;
		gate_size_ = rhs.gate_size_;
		area_oriented_ = rhs.area_oriented_;
		verbose_ = rhs.verbose_;
		return *this;
	}	
	int getLutSize() const {return lut_size_;}
	int getGateSize() const {return gate_size_;}
	bool isAreaOriented() const {return area_oriented_;}
	bool verbose() const {return verbose_;}
private:
  int  lut_size_;
  int  gate_size_;
  bool area_oriented_;
  bool verbose_;
};

//===----------------------------------------------------------------------===//
//                               NodeSol class
//===----------------------------------------------------------------------===//
class NodeSol {
public:
  typedef enum {
  Lut     = 0,
} sol_type;

NodeSol(sol_type type, Area area) : type(type), area(area), ll(0) { nodes.reserve(7); }
NodeSol(Area area) : area(area) {}

NodeSol& operator+=(const NodeSol&);

Level getLevel();

public:
  sol_type type = Lut;
  Area     area = 0.0;
  std::vector<Node*> nodes;
  pCut lut_sol;
  Level ll;
};


//===----------------------------------------------------------------------===//
//                               Node class
//===----------------------------------------------------------------------===//
class Node {
public:
typedef enum {
  Pi   = 0,
  And  = 1, // and
  Cho  = 2, // choice and
  Po   = 3,
  Cns  = 4  // const 1
} NodeType;

  Node(NodeType type, int id, const std::string& name, int fanout, bool phase, int cuts_num) : type_(type), fanout_num_(fanout), id_(id), name_(name), area_(0), flow_(0), min_dep_area_(nullptr),
                        simple_gate_(nullptr), mark_(-1), virtual_(false), phase_(phase), choice_root_(false)
  {
    fanins_.resize(2);
    choice.reserve(2);
    comple_[0] = false;
    comple_[1] = false;
    cuts_.reserve(cuts_num);
  }

  Node (int id) : Node(Node::And, id, "", 1, false, 1)
  {
    name_ = "n" + std::to_string(id);
    virtual_ = true;
  }
  
  ~Node() {}

public:
  Node *fanin (int idx) { return fanins_[idx];  }

  int faninNum()  const { return fanins_.size(); }
  int fanoutNum() const { return fanout_num_;    }

  bool complement(int idx){ return comple_[idx]; }

  NodeSol minAreaIncSol(Solution *sol, bool flow, Level required);

  bool isAnd()       const  { return type_ == And;}
  bool isPi()		     const  { return type_ == Pi; }
  bool isPo()        const  { return type_ == Po; }
  bool isConst()     const  { return type_ == Cns;}
  bool isChoice()    const  { return type_ == Cho;}
  bool isVirtual()   const  { return virtual_;    }
  bool phase()       const  { return phase_;      }
  int getId()        const  { return id_;         }
  int cutListNum()   const  { return cuts_.size();}

  pCut  minDepth()    const  { return min_dep_area_ == nullptr ? minArea() : min_dep_area_; }
  pCut  minArea()     const  { return cuts_.front();    }
  pCut  trivCut()     const  { return cuts_.back();     }
  Area  area()        const  { return area_;            }
  Area  flow()        const  { return flow_;            }
  pGate sGate()       const  { return simple_gate_;     }
  int   mark()        const  { return mark_;            }
  auto  begin()              { return cuts_.begin();    }
  auto  end()                { return cuts_.end();      }

  void setMark(int m)       { mark_ = m;    }
  void setArea(Area area)   { area_ = area; }
  void setFlow(Area flow)   { flow_ = flow; }

  void setChoice()          { type_ = Node::Cho; }

  void setTimingCut(const pCut& cut){ min_dep_area_ = cut; }

  const std::string& name() const { return name_; }

  Depth depth();
  
  void  cutEnum(int lut_size, bool area_oriented);
  void  p(std::string name); // for debug

  void  setSgate(const pGate& gate){ simple_gate_ = gate; }

  void addCut(const pCut& cut){ cuts_.emplace_back(cut); }
  void setAig(Node *fanin0, Node *fanin1, bool comple0, bool comple1);

  const pCut& getCut(int idx){ return cuts_[idx]; }

  pCut trivCutGen(bool area_mode);

  /**
   * @brief deference the nodes in AIG
   */
  void dereference();

  // void choiceMap(int fanout);

  void setFanoutNum(int fanout) { fanout_num_ = fanout; }

public:
  std::vector<Node*> choice; // the choice nodes of this node

private:
  NodeSol lutSol(Solution *sol, bool if_div, bool flow, Level required);

  pCut mergeCut(const pCut& left, const pCut& right, int lut_size, bool area_oriented);

  friend class Cut;
  friend class SimpleGate;
  friend class Decompose;

private:
  NodeType type_;
  NodeVec  fanins_;
  int      fanout_num_;
  bool     comple_[2];
  int      id_;
  std::string name_;
  CutList  cuts_;
  Area     area_; // effective area
  Area     flow_; // another estimation
  pCut     min_dep_area_;
  pGate    simple_gate_;
  int      mark_; // the multipurpose mark
  bool     virtual_;
  bool     phase_;
  bool     choice_root_;
};

struct comp_pnode_id
{
 bool operator()(const Node* lhs, const Node* rhs) { return lhs->getId() > rhs->getId(); }
};

struct pnode_id
{
 int operator()(const Node* node) { return node->getId(); }
};
//===----------------------------------------------------------------------===//
//                               Cut class
//===----------------------------------------------------------------------===//
class Cut {
public:
typedef enum {
  General = 0,
  Cut7    = 1,
  Cut8    = 2,
  Cut9    = 3,
} cut_type;

public:
  Cut(Node *root, Area af, Depth dep, Sign sign) 
  : area_(af), 
    depth_(dep), 
    root_(root), 
    sign_(sign),
    //func_(root->name())
    func_("n" + std::to_string(root_->getId()))
  { 
    inputs_.emplace_back(root);
  }
  Cut(Cut *cut) : area_(cut->area_), depth_(cut->depth_){ inputs_ = cut->inputs_; }
  Cut(const pGate& gate);
  Cut(Node *node) { inputs_.emplace_back(node->fanin(0)); inputs_.emplace_back(node->fanin(1)); }
  Cut(std::vector<Node*>& inputs, Area area, Depth depth, Node *root, Function& func, Sign sign) 
  : area_(area), depth_(depth), root_(root), sign_(sign) 
  {
    inputs_.swap(inputs);
    func_.swap(func);
  }

  Cut(const Cut& cut) = default;

public:
  Area      getArea()  const { return area_;           }
  cut_type  getType()  const { return type_;           }
  size_t    cutsize()  const { return inputs_.size() ; }
  auto      begin()    const { return inputs_.begin(); }
  auto      end()      const { return inputs_.end()  ; }
  
  bool operator==(const Cut& cut);
  pCut operator+ (const Cut& cut);
  
  Area localArea(bool if_div);
  Area localFlow();
  Area baseArea();

  Node *root() const { return root_; }
  void  setRoot(Node *root) { root_ = root; }

  float getDelay();

  std::string p(bool simple);

  const Function& func() const { return func_; }
  const Sign&     sign() const { return sign_; }

  const std::vector<Node*>& inputs() { return inputs_; }

  Depth depth() const { return depth_; }

  void invert() { func_ = ~func_; }

  friend class Node;

private:
  Area      area_  = 0;
  Depth     depth_ = 0;
  cut_type  type_  = General;
  Node     *root_  = nullptr;
  Sign      sign_  = 0;  
  Function  func_;
  std::vector<Node*> inputs_;
};

struct comp_pcut_area {
  bool operator()(const pCut& lhs, const pCut& rhs) const 
  {
    return lhs->getArea() < rhs->getArea();
  }
};

struct comp_pcut_delay {
  bool operator()(const pCut& lhs, const pCut& rhs) const 
  {
    if (lhs->depth() == rhs->depth())
    {
      return lhs->getArea() < rhs->getArea();
    } 
    return lhs->depth() < rhs->depth();
  }
};

struct sizer_pcut {
  size_t operator()(const pCut& cut) const { return cut->cutsize(); }
};

struct pcut_area {
  Area operator()(const pCut& cut) const { return cut->getArea(); }
};

static bool comp_ppcut_size_greater_func(const Cut* lhs, const Cut* rhs)
{ 
  if (lhs->cutsize() == rhs->cutsize())
  {
    return lhs->root()->getId() < rhs->root()->getId();
  }
  return lhs->cutsize() > rhs->cutsize();
}

static bool comp_ppcut_delay_func(const Cut* lhs, const Cut* rhs)
{
  if (lhs->depth() != rhs->depth())
  {
    return lhs->depth() < rhs->depth();
  }
  if (lhs->cutsize() != rhs->cutsize())
  {
    return lhs->cutsize() > rhs->cutsize();
  }
  return lhs->root()->getId() < rhs->root()->getId();
}
//===----------------------------------------------------------------------===//
//                               Mapper class
//===----------------------------------------------------------------------===//
class Mapper {
public:
  Mapper(Abc_Ntk_t *pNtk, Para& para) : pNtk_(pNtk), para_(para), best_sol_(nullptr), ideal_depth_(0)
  { 
    nodes_.reserve(pNtk->vObjs->nSize);
    Po_.reserve(Abc_NtkPoNum(pNtk));
    Pi_.reserve(Abc_NtkPiNum(pNtk));
    simple_gates_num_.resize(32);
  }

  ~Mapper();
public:
  int map();

  Abc_Ntk_t   *getNtk()     const { return pNtk_;                      }
  std::string  getNtkName() const { return std::string(pNtk_->pName);  }

  auto begin()  const { return nodes_.begin();  }
  auto end() 	  const { return nodes_.end(); 	  }
  auto rbegin() const { return nodes_.rbegin(); }
  auto rend()   const { return nodes_.rend();	  }

  int   poNum()             const { return Po_.size();  }
  Node *getPo(size_t idx)   const { return Po_[idx];    }

private:
  void cutSel(Solution *sol, bool area_oriented, bool flow);
  void itrSel(Solution*& initial_sol, int max_itr_num, bool area_oriented, bool flow);

  void updateFlow(Solution *sol);

  void printRel(Solution* sol);

  Depth getIdealDepth();
  
  void simpleGateGen(int gate_size_upper);

  void init();

  void dumpTempNetwork(); // absorb the buffer or inverter into po
  Abc_Ntk_t * agdmapToAbcLogic();
  std::string expToTrueValue(const std::vector<Node*>& inputs, const std::string & exp);

private:
  Abc_Ntk_t *pNtk_;
  Para para_;
  Solution  *best_sol_;
  Node      *const_1_ = nullptr;

  Depth ideal_depth_; // the min depth of this ntk can obtain
  
  std::vector<int> simple_gates_num_;

  std::vector<std::vector<Node*>> diff_level_nodes_;

  std::vector<Node*> nodes_; // store all nodes(include PI and exclude PO) with a topological order
  std::vector<Node*> nodes_with_virtual_;
  std::vector<Node*> Po_;
  std::vector<Node*> Pi_;
  std::map<Node*, Level> required_; // the required time for timing mode
  std::map<Node*, Node*> lut_pair_;

	bool isAreaOriented() const {return para_.isAreaOriented();}
	int getLutSize() const {return para_.getLutSize();}
	int getGateSize() const {return para_.getGateSize();}
	bool verbose() const {return para_.verbose();}
};



//===----------------------------------------------------------------------===//
//                               Solution class
//===----------------------------------------------------------------------===//
class Solution {
public:
  Solution(Mapper* mapper);
 ~Solution() {}

public:
  Level getLevel(); 

  int lutNum(int lut_size){ return lut_num_[lut_size - 1]; }

  int getArea() const { return lut_;   }
  
  void add(Node *node, const pCut& rep);

  bool isRoot(Node *node){ return sol_.find(node) != sol_.end();    }
  pCut repCut(Node *node){ assert(isRoot(node)); return sol_[node]; }

  bool shared(Node *node){ return shared_prohibited_.find(node) == shared_prohibited_.end(); }

  int   fanout(Node *node)    { return weight_.find(node) == weight_.end() ? 1 : weight_[node]; }
  int   reference(Node *node) { return weight_.find(node) == weight_.end() ? 0 : weight_[node]; }
  void  setRef(Node* node)    { assert(weight_.find(node) == weight_.end()); weight_[node] = 1; }

  Mapper* mapper;
  int rtl_f7_, rtl_f8_, rtl_f9_;

  friend class Mapper;

private:
  int lut_;
  Level logic_level_;

  std::vector<int> lut_num_; // the num of different size luts
  std::map<Node*, pCut> sol_;// roots and their representative cuts
  std::map<Node*, int> weight_;
  std::set<Node*> shared_prohibited_;
};


//===----------------------------------------------------------------------===//
//                               Bin class
//===----------------------------------------------------------------------===//
class Bin {
public:
  Bin(Cut *cut) 
  : root_(nullptr), 
    sign_(cut->sign()),
    depth_(0),
    inputs_(cut->inputs())
  {
    cuts_.reserve(3);
    cuts_.push_back(cut);
  }

  Bin() 
  : root_(nullptr),
    sign_(0),
    depth_(0)
  {}

  Bin(const Bin& bin) = default;

  Bin(Bin *bin)
  : root_(nullptr), 
    sign_(bin->sign_),
    depth_(bin->depth_),
    inputs_(bin->inputs_),
    cuts_(bin->cuts_)
  {
    connect_.reserve(bin->connect_.size());
  }
  
  ~Bin() {}

public:
  int cutNum()  const { return cuts_.size();   }
  int  size()   const { return inputs_.size(); }

  int cutsize() const { return inputs_.size() + connect_.size(); }

  int  connectSize() const { return connect_.size();  }
  void addConnect(Bin *bin){ connect_.push_back(bin); }

  const std::vector<Bin*>& connect() const { return connect_; }

  pCut cutGen(const std::set<Node*>& complemented, std::vector<Node*>& vir_nodes);

  int depth() const { return depth_; }
  void setDepth(int depth) { depth_ = depth; }

  int calculate_delay();

  Node *cutEnum(std::vector<Node*>& choice_nodes, const std::set<Node*>& complement, int lut_size, bool area_oriented);

  Node *cutsMerge(Node *fanin1, Node *fanin2, const std::set<Node*>& complement, int lut_size, bool area_oriented);

  friend class Decompose;

private:
  Node *root_;
  pCut  cut_;
  Sign  sign_;
  int   depth_;
  std::vector<Node*> inputs_;  
  std::vector<Cut*>  cuts_;    // the cuts(boxes)
  std::vector<Bin*>  connect_; // multi-level tree fanin
};

struct comp_pbin_size {
  bool operator()(const Bin *lhs, const Bin *rhs) const { return lhs->size() > rhs->size(); }
};

static bool comp_pbin_size_less(const Bin* lhs, const Bin* rhs){ return lhs->size() < rhs->size(); }

//===----------------------------------------------------------------------===//
//                               Decompose class
//===----------------------------------------------------------------------===//
class Decompose {
public:
  Decompose(const pCut& cut, int k, int gate_size)
  : cut_size_(0),
    lut_size_(k),
    area_((cut->getArea() - 1.0) / (Area)cut->root()->fanoutNum()),
    depth_(cut->depth()),
    root_bin_(nullptr),
    inputs_(cut->inputs())
  {
    cuts_.reserve(gate_size);
    cuts_.push_back(cut.get());
  }

  Decompose(Decompose *decomp, std::vector<Bin*>& bins) 
  : cut_size_(0), 
    lut_size_(decomp->lut_size_),
    area_(decomp->area_), 
    depth_(decomp->depth_),
    root_bin_(nullptr)
  {
    bins_.swap(bins);
  }

  Decompose(const Decompose& decom)
  : cut_size_(decom.cut_size_),
    lut_size_(decom.lut_size_),
    area_(decom.area_),
    depth_(decom.depth_),
    inputs_(decom.inputs_), 
    cuts_(decom.cuts_)
  {
    cuts_.reserve(cuts_.size() * 2);
  }

  Decompose(std::vector<Cut*>& cuts, int k)
  : cut_size_(0),
    lut_size_(k),
    area_(0),
    depth_(0),
    root_bin_(nullptr)
  {
    cuts_.reserve(cuts.size());
    for (Cut *cut : cuts)
    {
      area_ += (cut->getArea() - 1.0) / cut->root()->fanoutNum();
      if (depth_ < cut->depth())
      {
        depth_ = cut->depth();
      }
      cuts_.push_back(cut);
    }
  }

  ~Decompose() 
  {
    for (Bin *bin : bins_)
    {
      delete bin;
    }
  }

  /**
   * @brief combine the subcut 
   */
  void combine(const pCut& sub_cut)
  {
    area_ += (sub_cut->getArea() - 1.0) / (Area)sub_cut->root()->fanoutNum();
    std::vector<Node*> inputs;
    utility::ordered_merge<Node*, comp_pnode_id>(this->inputs_, sub_cut->inputs(), inputs);
    this->inputs_.swap(inputs);
    cuts_.emplace_back(sub_cut.get());    
    depth_ = depth_ > sub_cut->depth() ? depth_ : sub_cut->depth();
  }

  /**
   * @brief unique node inputs size
   */
  int input_size() const { return inputs_.size(); }

  /**
   * @brief the cone area of this wide cut
   */
  Area area() const { return area_; }

  int depth() const { return depth_; }

  int lut_size() const { return lut_size_; }

  /**
   * @brief the equivalent cut size of this decomposition's kcut
   * 
   * @return int cutsize 
   */
  int cutsize() const { return cut_size_; }

  /**
   * @brief delay-oriented bin-packing, method of a balanced strategy
   * 
   * @return pDec another decomposition which has a cutsize-area trad-off with this decomposition
   */
  pDec delay_binpack();
  void area_binpack();
  
  /**
   * @brief area-driven multi-level decompose
   */
  void multilevel_decompose();

  /**
   * @brief kcut generation
   */
  pCut kcut(Node *root, const std::set<Node*>& complemented, bool area_oriented);

  /**
   * @brief generating the cut of bins recursively
   * 
   * @param phase simple gate input phase
   * @return pCut kcut
   */
  pCut bins2cuts(const std::set<Node*>& complemented) { return root_bin_->cutGen(complemented, vir_nodes_); }

  /**
   * @brief print this Decompose to std::cout
   */
  void print();

  void cutEnum(std::vector<pCut>& cuts, const std::set<Node*>& complement, int lut_size, bool area_oriented);

  friend class SimpleGate;

private:
  // kcut cutsize
  int cut_size_;
  // the number of LUT inputs
  int lut_size_;
  // the cuts' cone area
  Area area_;
  // the depth of this decomposition
  int depth_;
  // root bin
  Bin *root_bin_;
  // inputs
  std::vector<Node*> inputs_;
  // the inputs cut solution
  std::vector<Cut*> cuts_;
  // all bins
  std::vector<Bin*> bins_;
  // the virtual internal nodes
  std::vector<Node*> vir_nodes_;
};

// delay-first decomposition comparision
struct comp_pdec_delay {
  bool operator()(const pDec& lhs, const pDec& rhs) const
  {
    if (lhs->depth() == rhs->depth())
    {
      return lhs->area() < rhs->area(); 
    }
    return lhs->depth() < rhs->depth();
  }
};

// area-first decomposition comparision
struct comp_pdec_area {
  bool operator()(const pDec& lhs, const pDec& rhs) const
  {
    return lhs->area() < rhs->area(); 
  }
};

static bool comp_pdec_area_fn(const pDec& lhs, const pDec& rhs) { return lhs->area() < rhs->area(); }

struct sizer_pdec {
  int operator()(const pDec& pdec) const { return pdec->input_size(); }
};

struct sizer_pdec_cutsize {
  int operator()(const pDec& pdec) const { return pdec->cutsize(); }
};

struct pdec_area {
  Area operator()(const pDec& pdec) const { return pdec->area(); }
};


//===----------------------------------------------------------------------===//
//                               SimpleGate class
//===----------------------------------------------------------------------===//
class SimpleGate {
public:
  SimpleGate(Node *root) : root_(root), trivial_(false) {}
  ~SimpleGate() {}

public:
  /**
   * @brief expand this simple gate to gate_size as possible
   * 
   * @param gate_size
   */
  void expand(int gate_size);

  // inputs
  const std::vector<Node*>& inputs() const { return inputs_; }

  int   size()    const { return inputs_.size(); }
  Node *root()    const { return root_;          }
  bool  trivial() const { return trivial_;       }

  /**
   * @brief simple gate based cut enumeration
   * 
   * @param k the max size of the general cut
   * @param vir_nodes virtual nodes 
   */
  void SimpleGateEnum(int k, bool area_oriented, std::vector<Node*>& vir_nodes);

  /**
   * @brief combine the inputs cuts with area-size pruning
   */
  void combine(int k, bool area_oriented);

  void print();

private:
  bool eat(Node *node, int idx);

private:
  Node *root_;
  // if all inputs of a simple gate are pi, then it is  trivial
  bool trivial_;
  // the decom solutions
  std::vector<pDec> decompositions_;
  // simple gate inputs
  std::vector<Node*> inputs_;
  // the interal nodes of this gate, not include gate root
  std::vector<Node*> internal_; 
  // the complemented inputs
  std::set<Node*> complemented_;
};

//===----------------------------------------------------------------------===//
//                               API for abc
//===----------------------------------------------------------------------===//


static int countOnes64(unsigned long uWord) 
{
  uWord = (uWord & 0x5555555555555555) + ((uWord>>1)  & 0x5555555555555555);
  uWord = (uWord & 0x3333333333333333) + ((uWord>>2)  & 0x3333333333333333);
  uWord = (uWord & 0x0F0F0F0F0F0F0F0F) + ((uWord>>4)  & 0x0F0F0F0F0F0F0F0F);
  uWord = (uWord & 0x00FF00FF00FF00FF) + ((uWord>>8)  & 0x00FF00FF00FF00FF);
  uWord = (uWord & 0x0000FFFF0000FFFF) + ((uWord>>16) & 0x0000FFFF0000FFFF);
  return  (uWord & 0x00000000FFFFFFFF) + (uWord>>32);
}

} // end namespace

#endif //_AGDMAP_H_
