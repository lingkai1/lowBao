#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"
//#include "Graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <string>
#include "float.h"
#include <fstream> //文件操作
#include <iterator>
#include <cstdlib>
#include <algorithm>

using namespace std;
class FindServer;

#define MAX_64 (0x7fffffffffffffffLL)
#define MAX_32 (0x7fffffff)

#define PRICEUP_UP_BOUND MAX_64
#define UNFEASIBLE          -2
#define ALLOCATION_FAULT    -3
#define PRICE_OFL           -4
#define UPDT_FREQ      0.4
#define UPDT_FREQ_S    30
#define SCALE_DEFAULT  12.0
#define PRICE_OUT_START  1
#define CUT_OFF_POWER    0.44
#define CUT_OFF_COEF     1.5
#define CUT_OFF_POWER2   0.75
#define CUT_OFF_COEF2    1
#define CUT_OFF_GAP      0.8
#define CUT_OFF_MIN      12
#define CUT_OFF_INCREASE 4
#define TIME_FOR_PRICE_IN1    2
#define TIME_FOR_PRICE_IN2    4
#define TIME_FOR_PRICE_IN3    6
#define MAX_CYCLES_CANCELLED 0
#define START_CYCLE_CANCEL   100
#define MAX( x, y ) ( ( (x) > (y) ) ?  x : y )
#define MIN( x, y ) ( ( (x) < (y) ) ? x : y )
#define ABS( x ) ( (x) >= 0 ) ? (x) : -(x)
#define NODE_NUMBER( i ) ( ( (i) == NULL ) ? -1 : ( (i) - nodeArray + nodeMin ) )
#define ARC_NUMBER( a ) ( ( (a) == NULL )? -1 : (a) - arcArray )
#define INT_TO_POS_NODE(i) (nodeArray + (i))

using namespace std;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned char uint8;

struct outputEdgeType
{
  int startNode;
  int endNode;
  int flow;
  int cost;
  bool operator == (const int &in)
  {

      if((startNode == in)&& (flow > 0))// 找到了具有该节点为起点的边，并且该边的flow不为0
          return 1;
      else
          return 0;
  }
};
struct serverDataType
{
	int capacity;
	int cost;
};
struct DirectedEdge
{
    friend bool operator < (DirectedEdge n1, DirectedEdge n2)
    {
//        降序
        return n1.cost > n2.cost;
    }

    int startNode;
    int endNode;
    int capacity;
    int residualcapacity;
    int cost;
    DirectedEdge()
    {
       startNode = 0;
       endNode = 0;
       capacity = 0.0;
       residualcapacity = 0;
       cost = 0;
    }
};
typedef struct
{
    int nodeID;
    int demand;
}ConsumerDatatype;

class graph
{
private:
    int v;              // numbel of node
    int e;              // number of arc

public:
    vector<vector<DirectedEdge>> adj; // 定义邻接
    graph(int v)
    {
        this->v = v;
        this->e = 0;
        vector<DirectedEdge> temp;
        for(int i = 0; i < v ; i++)
            adj.push_back(temp);

    }
    friend class PrimMST;
    friend class Kruskal;
    friend class Dijkstra;
    inline void addEdge(DirectedEdge edge)
    {
        adj[edge.startNode].push_back(edge);
        e++;
    }
    int V(){return v;}
    int E(){return e;}

    void read(vector<ConsumerDatatype> consumerVector);

    ~graph()
    {}
};

class FindServer
{
 public:
        typedef long long int excessType;
        typedef long long int priceType;

        class NODE;

        class ARC {
        public:
                long residualCapcity; // residual capacity;
                priceType _cost; // cost of arc;
                NODE *headNode; // head node;
                ARC *sisterArc; // opposite arc;
                int tailInt;
                int headInt;
        public:
                ARC() {}
                ~ARC() {}

                void setResidualCapacity( long acessResidualCapcity) { residualCapcity = acessResidualCapcity; }
                void decreaseResidualCapacity( long delta) { residualCapcity -= delta; }
                void increaseResidualCapacity( long delta) { residualCapcity += delta; }
                void setCost( priceType cost) { _cost = cost; }
                void multiplyCost( priceType mult) { _cost *= mult; }
                void setHead( NODE *head) { headNode = head; }
                void setSisterArc( ARC *sister) { sisterArc = sister; }
                long acessResidualCapcity() { return residualCapcity; }
                priceType cost() { return _cost; }
                NODE *head() { return headNode; }
                ARC *sister() { return sisterArc; }
        };

        class NODE {
        public:
                excessType EXCESS; // excess of the node;
                priceType _price; // distance from a sink;
                ARC *firstOutArc; //保存了当前第一条边在 arc链表里的位置
                ARC *currentOutArc; // 当前使用的边
                //
                ARC *suspendArc;
                NODE *queueNext; // next node in push queue
                NODE *bucketNext; // next node in bucket-list;
                NODE *bucketPre; // previous node in bucket-list;
                long _rank; // bucket number;
                long _inp; // auxilary field;
                long flewToSuperSink;
                double averageExcessCost;// 平均
                long long totalCost;
                long long overTotalCost;
                double averageRouteCost;
                long long totalCostMax;
                long long excessMax;
                int serverId;

        public:
                NODE() {}
                ~NODE() {}

                void setExcess( excessType excess) { EXCESS = excess; }
                void decreaseExcess( long delta) { EXCESS -= delta; }
                void increaseExcess( long delta) { EXCESS += delta; }
                void setPrice( priceType price) { _price = price; }
                void decreasePrice( long delta) { _price -= delta; }
                void setFirst( ARC *first) { firstOutArc = first; }
                void setCurrent( ARC *current) { currentOutArc = current; }
                void increaseCurrent() { currentOutArc ++; }
                void setSuspended( ARC *suspended) { suspendArc = suspended; }
                void setQueueNext( NODE *q_next) { queueNext = q_next; }
                void setBucketNext( NODE *b_next) { bucketNext = b_next; }
                void setBucketPre( NODE *b_prev) { bucketPre = b_prev; }
                void setRank( long rank) { _rank = rank; }
                void setInp( long inp) { _inp = inp; }
                excessType acessExs() { return EXCESS; }
                priceType price() { return _price; }
                ARC *first() { return firstOutArc; }
                void decreaseFirst() { firstOutArc --; }
                void increaseFirst() { firstOutArc ++; }
                ARC *current() { return currentOutArc; }
                ARC *suspended() { return suspendArc; }
                NODE *acessQNext() { return queueNext; }
                NODE *acessBNext() { return bucketNext; }
                NODE *acessBPre() { return bucketPre; }
                long rank() { return _rank; }
                long inp() { return _inp; }
        };

        class BUCKET {
        private:

                NODE *firstInBucket;
        public:
        BUCKET( NODE *ptrFirst) : firstInBucket(ptrFirst) {}
                BUCKET() {}
                ~BUCKET() {}

                void setPositiveFirst( NODE *ptrFirst) { firstInBucket = ptrFirst; }
                NODE *positiveFirt() { return firstInBucket; }
        };

 public:
        long _n; // number of nodes
        long _m; // number of arcs

        long *capArray;
        NODE *nodeArray;
        NODE *sentinelNode;
        NODE *excessQueueFirst;
        NODE *excessQueueLast;
        ARC *arcArray;
        ARC *sentinelArc;

        BUCKET *bucketArray;
        BUCKET *lastBucket;
        long numberOfBucket;
        int timeForPriceIn;

        priceType epsilon;
        priceType _dn;
        priceType PricelowBound;
        priceType multiMaxCost;
        double scaleFactor;
        double cutOffFactor;
        double cutOn; // the bound for returning suspended arcs
        double cutOff; // the bound for suspending arcs
        excessType totalExcess;


        int priceFlag;

        int updateFlag;

        int maxCycleCancel;


        ARC tempUseArc;
        NODE tempUseNode;
        NODE *AddressOftmpNode;
        NODE *dNode;

        long nRel;
        long nRef;
        long numberOfNodeWithExcess;
        long numberPush;
        long numberRelabel;
        long numberDiscarge;
        long numberRefine;
        long numberUpdate;
        long numberScan;
        long numberPreScan;
        long numberPreScan1;
        long numberPreScan2;
        long numberBadPriceIn;
        long numberBadRelabel;
        long numberPRefine;

        bool numberZeroCycles;
        bool checkSolution;
        bool computeDuals;
        bool costRestart;
        bool printAnswer;
        long long int *nodeBalance;

        // sketch variables used during reading in arcs;
        long nodeMin;
        long nodeMax;
    long *arcFirst;


        long *arcTail;
        long posCurrent;
        ARC *arcCurrent;
        ARC *arcNew;
        ARC *arcTmp;
        priceType maxCost; // maximum cost
        excessType _total_supply; // total supply
        excessType _total_demand; // total demand
        // pointers to the node structure
        NODE *iNode;
        NODE *jNode;

 public:
        FindServer( long numNodes, long numArcs) {
                _n = numNodes;
                _m = numArcs;

                priceFlag = 0;
                updateFlag = 0;
                numberPush = 0;
                numberRelabel = 0;
                numberDiscarge = 0;
                numberRefine = 0;
                numberUpdate = 0;
                numberScan = 0;
                numberPreScan = 0;
                numberPreScan1 = 0;
                numberPreScan2 = 0;
                numberBadPriceIn = 0;
                numberBadRelabel = 0;
                numberPRefine = 0;
                numberZeroCycles = false;
                checkSolution = false;
                computeDuals = false;
                costRestart = false;
                printAnswer = true;
                // allocate arrays and prepare for "receiving" arcs;
                // will also reset posCurrent, etc.;
                allocateMemory();
        }
        ~FindServer() {}

        void faultHandeler( int cc);
        void allocateMemory();
        void deleteMemory();
        void setArc( long tail_node_id, long head_node_id,
                                  long low_bound, long up_bound, priceType cost);
        void setDemandOfNode( long id, long excess);
        void preProscessing();
        int initialization();
        void nodeScan( NODE *i);
        void priceUpdate();
        int relabel( NODE *i);
        int discharge( NODE *i);
        int priceIn();
        int refine();
        int priceRefine();
        void computePrices();
        void priceOut();
        int updateEpsilon();
        int checkFeasible();
        int checkCs();
        int checkEpsOpt();
        void initSolution();
        void csCostReinit();
        int cs2CostRestart( double *objective_cost);
        void printSolution();
        void printGraph();
        void finishup( double *objective_cost);
        int cs2( double *objective_cost);
        int runSolition();

        // shared utils;
        void increase_flow( NODE *i, NODE *j, ARC *a, long df);
        bool timeForUpdate() {
                return ( nRel > _n * UPDT_FREQ + numberOfNodeWithExcess * UPDT_FREQ_S);
        }
        // utils for excess queue;
        void rstExcessQueue() {
                for ( ; excessQueueFirst != NULL; excessQueueFirst = excessQueueLast ) {
                        excessQueueLast = excessQueueFirst->acessQNext();
                        excessQueueFirst->setQueueNext( sentinelNode);
                }
        }
//        如果它的next是_sentinel节点，就说明它不在excess 队列中
        bool notInExcessQueue( NODE *i) { return ( i->acessQNext() == sentinelNode); }
        bool isEmptyExcessQ() { return ( excessQueueFirst == NULL); }
        bool isNotEmptyExcessQ() { return ( excessQueueFirst != NULL); }
//        FIFO
        void insertToExcessQueue( NODE *i) {
                if ( isNotEmptyExcessQ() ) {
                        excessQueueLast->setQueueNext( i);
                } else {
                        excessQueueFirst = i;
                }
                i->setQueueNext( NULL);
                excessQueueLast = i;
        }
        void insertToFrontExcessQueue( NODE *i) {
                if ( isEmptyExcessQ() ) {
                        excessQueueLast = i;
                }
                i->setQueueNext( excessQueueFirst);
                excessQueueFirst = i;
        }
        void rmFromExcessQ( NODE *i) {
                i = excessQueueFirst;
                excessQueueFirst = i->acessQNext();
                i->setQueueNext( sentinelNode);
        }

        bool isEmptyStackQ() { return isEmptyExcessQ(); }
        bool isNotEmptyStackQ() { return isNotEmptyExcessQ(); }
        void rstStackQ() { rstExcessQueue(); }
        void stackQueuePush( NODE *i) {
                i->setQueueNext( excessQueueFirst);
                excessQueueFirst = i;
        }
        void stackQueuePop( NODE *i) {
                rmFromExcessQ( i);
        }
        // utils for buckets;
        void rstBucket( BUCKET *b) { b->setPositiveFirst( dNode); }
        bool isNotEmptyBucket( BUCKET *b) { return ( (b->positiveFirt()) != dNode); }
        void insertToBucket( NODE *i, BUCKET *b) {
                i->setBucketNext( b->positiveFirt() );
                b->positiveFirt()->setBucketPre( i);
                b->setPositiveFirst( i);
        }
        void getFromBucket( NODE *i, BUCKET *b) {
                i = b->positiveFirt();
                b->setPositiveFirst( i->acessBNext());
        }
        void rmFromBucket( NODE *i, BUCKET *b) {
                if ( i == b->positiveFirt() ) {
                        b->setPositiveFirst( i->acessBNext());
                } else {
                        i->acessBPre()->setBucketNext( i->acessBNext());
                        i->acessBNext()->setBucketPre( i->acessBPre());
                }
        }

        void updateCutOff() {
                if ( numberBadPriceIn + numberBadRelabel == 0) {
                        cutOffFactor = CUT_OFF_COEF2 * pow( (double)_n, CUT_OFF_POWER2 );
                        cutOffFactor = MAX ( cutOffFactor, CUT_OFF_MIN );
                        cutOff = cutOffFactor * epsilon;
                        cutOn = cutOff * CUT_OFF_GAP;
                } else {
                        cutOffFactor *= CUT_OFF_INCREASE;
                        cutOff = cutOffFactor * epsilon;
                        cutOn = cutOff * CUT_OFF_GAP;
                }
        }
        void exchange( ARC *a, ARC *b) {
                if ( a != b) {
                        ARC *sa = a->sister();
                        ARC *sb = b->sister();
                        long dCap;

                        tempUseArc.setResidualCapacity( a->acessResidualCapcity());
                        tempUseArc.setCost( a->cost());
                        tempUseArc.setHead( a->head());

                        a->setResidualCapacity( b->acessResidualCapcity());
                        a->setCost( b->cost());
                        a->setHead( b->head());

                        b->setResidualCapacity( tempUseArc.acessResidualCapcity());
                        b->setCost( tempUseArc.cost());
                        b->setHead( tempUseArc.head());

                        if ( a != sb) {
                                b->setSisterArc( sa);
                                a->setSisterArc( sb);
                                sa->setSisterArc( b);
                                sb->setSisterArc( a);
                        }

                        dCap = capArray[ a - arcArray];
                        capArray[ a - arcArray] = capArray[ b - arcArray];
                        capArray[ b - arcArray] = dCap;
                }
        }
};

class MinCostSolver
{
 public:
        typedef long long int excessType;
        typedef long long int priceType;

        class NODE;

        class ARC {
        public:
                long residualCapcity; // residual capacity;
                priceType _cost; // cost of arc;
                NODE *headNode; // head node;
                ARC *sisterArc; // opposite arc;
        public:
                ARC() {}
                ~ARC() {}

                void setResidualCapacity( long acessResidualCapcity) { residualCapcity = acessResidualCapcity; }
                void decreaseResidualCapacity( long delta) { residualCapcity -= delta; }
                void increaseResidualCapacity( long delta) { residualCapcity += delta; }
                void setCost( priceType cost) { _cost = cost; }
                void multiplyCost( priceType mult) { _cost *= mult; }
                void setHead( NODE *head) { headNode = head; }
                void setSisterArc( ARC *sister) { sisterArc = sister; }
                long acessResidualCapcity() { return residualCapcity; }
                priceType cost() { return _cost; }
                NODE *head() { return headNode; }
                ARC *sister() { return sisterArc; }
        };

        class NODE {
        public:
                excessType EXCESS;
                priceType _price;
                ARC *firstOutArc;
                ARC *currentOutArc;
                ARC *suspendArc;
                NODE *queueNext;
                NODE *bucketNext;
                NODE *bucketPre;
                long _rank;
                long _inp;
        public:
                NODE() {}
                ~NODE() {}

                void setExcess( excessType excess) { EXCESS = excess; }
                void decreaseExcess( long delta) { EXCESS -= delta; }
                void increaseExcess( long delta) { EXCESS += delta; }
                void setPrice( priceType price) { _price = price; }
                void decreasePrice( long delta) { _price -= delta; }
                void setFirst( ARC *first) { firstOutArc = first; }
                void setCurrent( ARC *current) { currentOutArc = current; }
                void increaseCurrent() { currentOutArc ++; }
                void setSuspended( ARC *suspended) { suspendArc = suspended; }
                void setQueueNext( NODE *q_next) { queueNext = q_next; }
                void setBucketNext( NODE *b_next) { bucketNext = b_next; }
                void setBucketPre( NODE *b_prev) { bucketPre = b_prev; }
                void setRank( long rank) { _rank = rank; }
                void setInp( long inp) { _inp = inp; }
                excessType acessExs() { return EXCESS; }
                priceType price() { return _price; }
                ARC *first() { return firstOutArc; }
                void decreaseFirst() { firstOutArc --; }
                void increaseFirst() { firstOutArc ++; }
                ARC *current() { return currentOutArc; }
                ARC *suspended() { return suspendArc; }
                NODE *acessQNext() { return queueNext; }
                NODE *acessBNext() { return bucketNext; }
                NODE *acessBPre() { return bucketPre; }
                long rank() { return _rank; }
                long inp() { return _inp; }
        };

        class BUCKET {
        private:
                // 1st node with positive excess or simply 1st node in the buket;
                NODE *firstInBucket;
        public:
        BUCKET( NODE *ptrFirst) : firstInBucket(ptrFirst) {}
                BUCKET() {}
                ~BUCKET() {}

                void setPositiveFirst( NODE *ptrFirst) { firstInBucket = ptrFirst; }
                NODE *positiveFirt() { return firstInBucket; }
        };

 public:
        long _n; // number of nodes
        long _m; // number of arcs

        long *capArray;
        NODE *nodeArray;
        NODE *sentinelNode;
        NODE *excessQueueFirst;
        NODE *excessQueueLast;
        ARC *arcArray;
        ARC *sentinelArc;

        BUCKET *bucketArray;
        BUCKET *lastBucket;
        long numberOfBucket;
        int timeForPriceIn;

        priceType epsilon;
        priceType _dn;
        priceType PricelowBound;
        priceType multiMaxCost;
        double scaleFactor;
        double cutOffFactor;
        double cutOn;
        double cutOff;
        excessType totalExcess;


        int priceFlag;

        int updateFlag;

        int maxCycleCancel;


        ARC tempUseArc;
        NODE tempUseNode;
        NODE *AddressOftmpNode;
        NODE *dNode;

        long nRel;
        long nRef;
        long numberOfNodeWithExcess;
        long numberPush;
        long numberRelabel;
        long numberDiscarge;
        long numberRefine;
        long numberUpdate;
        long numberScan;
        long numberPreScan;
        long numberPreScan1;
        long numberPreScan2;
        long numberBadPriceIn;
        long numberBadRelabel;
        long numberPRefine;

        bool numberZeroCycles;
        bool checkSolution;
        bool computeDuals;
        bool costRestart;
        bool printAnswer;
        long long int *nodeBalance;


        long nodeMin;
        long nodeMax;
    long *arcFirst;
        long *arcTail;
        long posCurrent;
        ARC *arcCurrent;
        ARC *arcNew;
        ARC *arcTmp;
        priceType maxCost; // maximum cost
        excessType _total_supply; // total supply
        excessType _total_demand; // total demand
        // pointers to the node structure
        NODE *iNode;
        NODE *jNode;

 public:
        MinCostSolver( long numNodes, long numArcs) {
                _n = numNodes;
                _m = numArcs;

                priceFlag = 0;
                updateFlag = 0;
                numberPush = 0;
                numberRelabel = 0;
                numberDiscarge = 0;
                numberRefine = 0;
                numberUpdate = 0;
                numberScan = 0;
                numberPreScan = 0;
                numberPreScan1 = 0;
                numberPreScan2 = 0;
                numberBadPriceIn = 0;
                numberBadRelabel = 0;
                numberPRefine = 0;
                numberZeroCycles = false;
                checkSolution = false;
                computeDuals = false;
                costRestart = false;
                printAnswer = true;
                // allocate arrays and prepare for "receiving" arcs;
                // will also reset posCurrent, etc.;
                allocateMemory();
        }
        ~MinCostSolver()
        {
            deleteMemory();
        }

        void faultHandeler( int cc);
        void allocateMemory();
        void deleteMemory();
        void setArc( long tail_node_id, long head_node_id,
                                  long low_bound, long up_bound, priceType cost);
        void setDemandOfNode( long id, long excess);
        void preProscessing();
        int initialization();
        void nodeScan( NODE *i);
        void priceUpdate();
        int relabel( NODE *i);
        int discharge( NODE *i);
        int priceIn();
        int refine();
        int priceRefine();
        void computePrices();
        void priceOut();
        int updateEpsilon();
        int checkFeasible();
        int checkCs();
        int checkEpsOpt();
        void initSolution();
        void csCostReinit();
        int cs2CostRestart( double *objective_cost);
        void printSolution();
        void printGraph();
        void finishup( double *objective_cost);
        int cs2( double *objective_cost);
        int runSolition();

        // shared utils;
        void increase_flow( NODE *i, NODE *j, ARC *a, long df) {
                i->decreaseExcess( df);
                j->increaseExcess( df);
                a->decreaseResidualCapacity( df);
                a->sister()->increaseResidualCapacity( df);
        }
        bool timeForUpdate() {
                return ( nRel > _n * UPDT_FREQ + numberOfNodeWithExcess * UPDT_FREQ_S);
        }
        // utils for excess queue;
        void rstExcessQueue() {
                for ( ; excessQueueFirst != NULL; excessQueueFirst = excessQueueLast ) {
                        excessQueueLast = excessQueueFirst->acessQNext();
                        excessQueueFirst->setQueueNext( sentinelNode);
                }
        }
//        如果它的next是_sentinel节点，就说明它不在excess 队列中
        bool notInExcessQueue( NODE *i) { return ( i->acessQNext() == sentinelNode); }
        bool isEmptyExcessQ() { return ( excessQueueFirst == NULL); }
        bool isNotEmptyExcessQ() { return ( excessQueueFirst != NULL); }
//        FIFO
        void insertToExcessQueue( NODE *i) {
                if ( isNotEmptyExcessQ() ) {
                        excessQueueLast->setQueueNext( i);
                } else {
                        excessQueueFirst = i;
                }
                i->setQueueNext( NULL);
                excessQueueLast = i;
        }
        void insertToFrontExcessQueue( NODE *i) {
                if ( isEmptyExcessQ() ) {
                        excessQueueLast = i;
                }
                i->setQueueNext( excessQueueFirst);
                excessQueueFirst = i;
        }
        void rmFromExcessQ( NODE *i) {
                i = excessQueueFirst;
                excessQueueFirst = i->acessQNext();
                i->setQueueNext( sentinelNode);
        }
        // utils for excess queue as a stack;
        bool isEmptyStackQ() { return isEmptyExcessQ(); }
        bool isNotEmptyStackQ() { return isNotEmptyExcessQ(); }
        void rstStackQ() { rstExcessQueue(); }
        void stackQueuePush( NODE *i) {
                i->setQueueNext( excessQueueFirst);
                excessQueueFirst = i;
        }
        void stackQueuePop( NODE *i) {
                rmFromExcessQ( i);
        }
        // utils for buckets;
        void rstBucket( BUCKET *b) { b->setPositiveFirst( dNode); }
        bool isNotEmptyBucket( BUCKET *b) { return ( (b->positiveFirt()) != dNode); }
        void insertToBucket( NODE *i, BUCKET *b) {
                i->setBucketNext( b->positiveFirt() );
                b->positiveFirt()->setBucketPre( i);
                b->setPositiveFirst( i);
        }
        void getFromBucket( NODE *i, BUCKET *b) {
                i = b->positiveFirt();
                b->setPositiveFirst( i->acessBNext());
        }
        void rmFromBucket( NODE *i, BUCKET *b) {
                if ( i == b->positiveFirt() ) {
                        b->setPositiveFirst( i->acessBNext());
                } else {
                        i->acessBPre()->setBucketNext( i->acessBNext());
                        i->acessBNext()->setBucketPre( i->acessBPre());
                }
        }
        // misc utils;
        void updateCutOff() {
                if ( numberBadPriceIn + numberBadRelabel == 0) {
                        cutOffFactor = CUT_OFF_COEF2 * pow( (double)_n, CUT_OFF_POWER2 );
                        cutOffFactor = MAX ( cutOffFactor, CUT_OFF_MIN );
                        cutOff = cutOffFactor * epsilon;
                        cutOn = cutOff * CUT_OFF_GAP;
                } else {
                        cutOffFactor *= CUT_OFF_INCREASE;
                        cutOff = cutOffFactor * epsilon;
                        cutOn = cutOff * CUT_OFF_GAP;
                }
        }
        void exchange( ARC *a, ARC *b) {
                if ( a != b) {
                        ARC *sa = a->sister();
                        ARC *sb = b->sister();
                        long dCap;

                        tempUseArc.setResidualCapacity( a->acessResidualCapcity());
                        tempUseArc.setCost( a->cost());
                        tempUseArc.setHead( a->head());

                        a->setResidualCapacity( b->acessResidualCapcity());
                        a->setCost( b->cost());
                        a->setHead( b->head());

                        b->setResidualCapacity( tempUseArc.acessResidualCapcity());
                        b->setCost( tempUseArc.cost());
                        b->setHead( tempUseArc.head());

                        if ( a != sb) {
                                b->setSisterArc( sa);
                                a->setSisterArc( sb);
                                sa->setSisterArc( b);
                                sb->setSisterArc( a);
                        }

                        dCap = capArray[ a - arcArray];
                        capArray[ a - arcArray] = capArray[ b - arcArray];
                        capArray[ b - arcArray] = dCap;
                }
        }
};
typedef struct
{
    int nodeNum;
    int demand;
}consumerDataType;

struct  averageExcessCostType
{
    int nodeID;
    double averageExcessCost;
    int nodeBaseCost;
    int serverID;
    int adjacentRank;
    int adjacentNode;
    long excessMax;
    int averageRouteCost;
    friend bool operator  <  (averageExcessCostType in1 ,averageExcessCostType in2)
    {
//        升序
        return in1.averageExcessCost < in2.averageExcessCost;
    }
    bool operator == (const int &in)
    {
        return nodeID == in;
    }
};

void deploy_server(char * graph[MAX_EDGE_NUM], int edge_num, char * filename);
void refineTheServer(graph &G, FindServer &findServerPro);
int findMinExcessCost(vector<averageExcessCostType> &vin, vector<DirectedEdge> &edgeVCec, FindServer &findServerPro);
int refineAccordingRank(vector<averageExcessCostType> &vin,FindServer &FindServer);
int findFuncForVector(vector<int> &vin,int val);
int erasecheck(int N,vector<int> &serverVec, graph &G);
int fitServer(int excess);
long calculateCost(vector<int> &serverNode);
#endif
