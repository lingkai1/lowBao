#include "deploy.h"
#include <stdio.h>
#include <string.h>

using namespace std;
double total_cost;
int currentNode;
int nodeNum;
int edgeNum;
int consumerNum;
int serverCostMax = 0;
int serverCostMid = 0;
int serverCostMin = 0x7FFF;
int serverCount;
vector<int> serverNode;
vector<serverDataType> serverCostVect; //  存储服务器的费用
vector<int> nodeBaseCost ;  //  存储节点的服务器架设费用
/////////////////////////////////////////////////////////////
//vector<double> averageExcessCostVec;
vector<averageExcessCostType>  averageExcessCostVecmoni;
vector<averageExcessCostType>  averageExcessCostVecmoniRaw;
vector<averageExcessCostType>  refineVect;
vector<outputEdgeType> usedARC;
int findInVect(vector<int> &vec, int val)
{
    vector<int>::iterator it;
    it = find(vec.begin(),vec.end(),val);
    if(it == vec.end())
        return 0;
    else
        return 1;
}
int findInVect(vector<ConsumerDatatype> &vec, int val)
{
    vector<ConsumerDatatype>::iterator it;
    for(it = vec.begin(); it < vec.end(); it++)
    {
        if((*it).nodeID == val)
        return 1;
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
//minCost
////////////////////////////////////////////////////////////////////////////////
#define INT_MAX 0x7FFF
#define WHITE 0
#define GREY  1
#define BLACK 2
#define OPEN( a ) ( a->acessResidualCapcity() > 0 )
#define CLOSED( a ) ( a->acessResidualCapcity() <= 0 )
#define REDUCED_COST( i, j, a ) ( i->price() + a->cost() - j->price() )
#define FEASIBLE( i, j, a ) ( i->price() + a->cost() < j->price() )
#define ADMISSIBLE( i, j, a ) ( OPEN( a ) && FEASIBLE( i, j, a ) )
#define SUSPENDED( i, a ) ( a < i->first() )


#define REMOVE_FROM_EXCESS_Q( i )				\
    {											\
        i           = excessQueueFirst;				\
        excessQueueFirst = i -> acessQNext();			\
        i ->setQueueNext( sentinelNode );		\
    }

#define STACKQ_POP( i )							\
    {											\
        i           = excessQueueFirst;				\
        excessQueueFirst = i -> acessQNext();			\
        i ->setQueueNext( sentinelNode );		\
    }

#define GET_FROM_BUCKET( i, b )					\
    {											\
        i    = ( b -> positiveFirt() );				\
        b ->setPositiveFirst( i -> acessBNext() );		\
    }
#define REMOVE_FROM_BUCKET( i, b )								\
    {															\
        if ( i == ( b -> positiveFirt() ) )							\
            b ->setPositiveFirst( i -> acessBNext() );					\
        else													\
            {													\
                ( i -> acessBPre() )->setBucketNext( i -> acessBNext() );	\
                ( i -> acessBNext() )->setBucketPre( i -> acessBPre() );	\
            }													\
    }

void MinCostSolver::faultHandeler( int cc)
{
    // abnormal finish
    printf ("\nError %d\n", cc );
    // 2 - problem is unfeasible
    // 5 - allocation fault
    // 6 - price overflow
    exit( cc);
}

void MinCostSolver::allocateMemory()
{
    // (1) allocate memory for 'nodes', 'arcs' and internal arrays;
//  后面需要设置一个  Sentinel Node  另一个是为了跳过 初始0点
    nodeArray = (NODE*) calloc ( _n+2,   sizeof(NODE) );
//    边数乘2   1是以为设置 边 Sentinel arc
    arcArray = (ARC*)  calloc ( 2*_m+1, sizeof(ARC) );
//    同样乘二 保存了当前边和sister 边的 capacity
    capArray = (long*) calloc ( 2*_m,   sizeof(long) );
//保存了sister 和当前边的tail node
    arcTail = (long*) calloc ( 2*_m,   sizeof(long) );
//    当前arc 所处的位置索引， 访问 capArray 和 _arc数组
    arcFirst = (long*) calloc ( _n+2,   sizeof(long) );
    // arcfirstOutArc [ 0 .. n+1 ] = 0 - initialized by calloc;
//设置所有的 excess 为0
    for ( NODE *in = nodeArray; in <= nodeArray + _n; in ++ ) {
        in->setExcess( 0);
    }
    if ( nodeArray == NULL || arcArray == NULL || arcFirst == NULL || arcTail == NULL) {
        printf("Error:  allocation error\n");
        exit( 1);
    }

    // (2) resets;
    posCurrent = 0;
    arcCurrent = arcArray; // set "current" pointer to the first arc *_arc 是初始位置
//    这两个先设置方向最大，一会改过来
    nodeMax = 0;
    nodeMin = _n;
    maxCost = 0;
    _total_supply = _total_demand = 0;
    // at this moment we are ready to add arcs and build the network,
    // by using setArc()...
}

void MinCostSolver::deleteMemory()
{
    if ( arcArray) free ( arcArray );
    if ( dNode) delete dNode;
    if ( capArray) free ( capArray );
    if ( bucketArray) free ( bucketArray );
    if ( checkSolution == true) free ( nodeBalance );
    if ( nodeArray) {
        nodeArray = nodeArray - nodeMin;
        free ( nodeArray );
    }
}

void MinCostSolver::setArc( long tail_node_id, long head_node_id,
                        long low_bound, long up_bound, // up_bound is basically capacity;
                        priceType cost)
{
    // DIMACS format:
    // c arc has <tail> <head> <capacity l.b.> <capacity u.b> <cost>

    if ( tail_node_id < 0 || tail_node_id > _n ||
         head_node_id < 0 || head_node_id > _n ) {
        printf("Error:  Arc with head or tail out of bounds\n");
        exit( 1);
    }
    if ( up_bound < 0 ) {
        up_bound = MAX_32;
        printf("Warning:  Infinite capacity\n");
    }
    if ( low_bound < 0 || low_bound > up_bound ) {
        printf("Error:  Wrong capacity bounds\n");
        exit( 1);
    }

    arcFirst[tail_node_id + 1] ++;
    arcFirst[head_node_id + 1] ++;
    iNode = nodeArray + tail_node_id;
    jNode = nodeArray + head_node_id;

    // store information about the arc   sister 的 tail 和当前tail 一起存储
    arcTail[posCurrent]   = tail_node_id;
    arcTail[posCurrent+1] = head_node_id;
//   只进行sethead 操作
    arcCurrent->setHead( jNode );
    arcCurrent->setResidualCapacity( up_bound - low_bound );
    capArray[posCurrent] = up_bound;
    arcCurrent->setCost( cost );
    arcCurrent->setSisterArc( arcCurrent + 1 );
//    设置 sister边（残存反向边）
    ( arcCurrent + 1 )->setHead( nodeArray + tail_node_id );
    ( arcCurrent + 1 )->setResidualCapacity( 0 );
    capArray[posCurrent+1] = 0;
//    sister边的cost和原cost 相反
    ( arcCurrent + 1 )->setCost( -cost );
//    两者互为sister
    ( arcCurrent + 1 )->setSisterArc( arcCurrent );

    iNode->decreaseExcess( low_bound );
    jNode->increaseExcess( low_bound );

    //
    if ( head_node_id < nodeMin ) nodeMin = head_node_id;
    if ( tail_node_id < nodeMin ) nodeMin = tail_node_id;
    if ( head_node_id > nodeMax ) nodeMax = head_node_id;
    if ( tail_node_id > nodeMax ) nodeMax = tail_node_id;

    if ( cost < 0 ) cost = -cost;
    if ( cost > maxCost && up_bound > 0 ) maxCost = cost;
    //  每次添加都添加当前边和sister的边
    arcCurrent += 2;
    posCurrent += 2;
}

void MinCostSolver::setDemandOfNode( long id, long excess)
{
    // set supply and demand of nodes; not used for transhipment nodes;
    if ( id < 0 || id > _n ) {
        printf("Error:  Unbalanced problem\n");
        exit( 1);
    }
//    可以设置负的excess
    (nodeArray + id)->setExcess( excess);
    if ( excess > 0) _total_supply += excess;
    if ( excess < 0) _total_demand -= excess;
}

void MinCostSolver::preProscessing()
{

    long i;
    long last, arc_num, arc_new_num;;
    long tail_node_id;
    NODE *head_p;
    ARC *arc_new, *arc_tmp;
    long up_bound;
    priceType cost; // arc cost;
    excessType cap_out; // sum of outgoing capacities
    excessType cap_in; // sum of incoming capacities

    if ( ABS( _total_supply - _total_demand ) > 0.5 ) {
        printf("Error:  Unbalanced problem\n");
        exit( 1);
    }

    // first arc from the first node
    ( nodeArray + nodeMin )->setFirst( arcArray );


//    node_min 一般是0
    for ( i = nodeMin + 1; i <= nodeMax + 1; i ++ ) {
//
        arcFirst[i] += arcFirst[i-1];
        ( nodeArray + i )->setFirst( arcArray + arcFirst[i] );
    }

    // scanning all the nodes except the last
    for ( i = nodeMin; i < nodeMax; i ++ )
    {
//        first - arc基地址？
        last = ( ( nodeArray + i + 1 )->first() ) - arcArray;

        for ( arc_num = arcFirst[i]; arc_num < last; arc_num ++ ) {
            tail_node_id = arcTail[arc_num];

            while ( tail_node_id != i ) {

                arc_new_num = arcFirst[tail_node_id];
                arcCurrent = arcArray + arc_num;
                arc_new = arcArray + arc_new_num;



                head_p = arc_new->head();
                arc_new->setHead( arcCurrent->head() );
                arcCurrent->setHead( head_p );

                up_bound          = capArray[arc_new_num];
                capArray[arc_new_num] = capArray[arc_num];
                capArray[arc_num]     = up_bound;

                up_bound = arc_new->acessResidualCapcity();
                arc_new->setResidualCapacity( arcCurrent->acessResidualCapcity() );
                arcCurrent->setResidualCapacity( up_bound) ;

                cost = arc_new->cost();
                arc_new->setCost( arcCurrent->cost() );
                arcCurrent->setCost( cost );

                if ( arc_new != arcCurrent->sister() ) {
                    arc_tmp = arc_new->sister();
                    arc_new->setSisterArc( arcCurrent->sister() );
                    arcCurrent->setSisterArc( arc_tmp );

                    arcCurrent->sister()->setSisterArc( arcCurrent );
                    arc_new->sister()->setSisterArc( arc_new );
                }

                arcTail[arc_num] = arcTail[arc_new_num];
                arcTail[arc_new_num] = tail_node_id;

                arcFirst[tail_node_id] ++ ;

                tail_node_id = arcTail[arc_num];
            }
        }

    }




    for ( NODE *ndp = nodeArray + nodeMin; ndp <= nodeArray + nodeMax; ndp ++ ) {
        cap_in  =   ( ndp->acessExs() );
        cap_out = - ( ndp->acessExs() );
        for ( arcCurrent = ndp->first(); arcCurrent != (ndp+1)->first();
              arcCurrent ++ ) {
            arc_num = arcCurrent - arcArray;
            if ( capArray[arc_num] > 0 ) cap_out += capArray[arc_num];
            if ( capArray[arc_num] == 0 )
                cap_in += capArray[ arcCurrent->sister() - arcArray ];
        }
    }
    if ( nodeMin < 0 || nodeMin > 1 ) {
        printf("Error:  Node ids must start from 0 or 1\n");
        exit( 1);
    }


    _n = nodeMax - nodeMin + 1;

    nodeArray = nodeArray + nodeMin;

    free ( arcFirst );
    free ( arcTail );
}

int MinCostSolver:: initialization()
{
    // initialization;
    // called after allocateMemory() and all nodes and arcs have been inputed;

    NODE *i; // current node
    ARC *a; // current arc
    ARC *a_stop;
    BUCKET *b; // current bucket
    long df;

    scaleFactor = (long) SCALE_DEFAULT;
//    设置末尾节点为 哨兵节点
    sentinelNode = nodeArray + _n;
    sentinelArc  = arcArray + _m;
//对所有节点进行初始化操作
    for ( i = nodeArray; i != sentinelNode; i ++ ) {
//        设置价格为0
        i->setPrice( 0);
//       suspend 节点什么意思
        i->setSuspended( i->first());
//        push 队列里放 放 哨兵
        i->setQueueNext( sentinelNode);
    }

    sentinelNode->setFirst( sentinelArc);
    sentinelNode->setSuspended( sentinelArc);

    // saturate negative arcs, e.g. in the circulation problem case
    for ( i = nodeArray; i != sentinelNode; i ++ ) {
        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {
            if ( a->cost() < 0) {
                if ( ( df = a->acessResidualCapcity()) > 0) {
                    increase_flow( i, a->head(), a, df);
                }
            }
        }
    }
//  节点数量+1
    _dn = _n + 1;
//    if ( numberZeroCycles == true) { // NO_ZERO_CYCLES
//        _dn = 2 * _dn;
//    }
//  遍历所有的arc
    for ( a = arcArray; a != sentinelArc; a ++ ) {
//  把所有边的cost 增大？
        a->multiplyCost( _dn);
    }


    if ((double) maxCost * (double) _dn > MAX_64) {
        printf("Warning:  Arc lengths too large, overflow possible\n");
    }
    multiMaxCost = maxCost * _dn;

    numberOfBucket = (long) (_dn * ceil(scaleFactor) + 2);

    bucketArray = (BUCKET*) calloc ( numberOfBucket, sizeof(BUCKET));
    if ( bucketArray == NULL )
        return ALLOCATION_FAULT;
//得到最后一个bucket的地址
    lastBucket = bucketArray + numberOfBucket;

    dNode = new NODE; // used as reference;
//init bucket
    for ( b = bucketArray; b != lastBucket; b ++ ) {
        rstBucket( b);
    }

    epsilon = multiMaxCost;
    if ( epsilon < 1) {
        epsilon = 1;
    }

    PricelowBound = -PRICEUP_UP_BOUND;

    cutOffFactor = CUT_OFF_COEF * pow( (double)_n, CUT_OFF_POWER);

    cutOffFactor = MAX( cutOffFactor, CUT_OFF_MIN);

    nRef = 0;

    priceFlag = 0;

    AddressOftmpNode = &tempUseNode;

    excessQueueFirst = NULL;

    return 0;
    //printGraph(); // debug;
}

void MinCostSolver::nodeScan( NODE *i)
{
    NODE *j; // opposite node
    ARC *a; // (i, j)
    ARC *a_stop; // first arc from the next node
    ARC *ra; // (j, i)
    BUCKET *b_old; // old bucket contained j
    BUCKET *b_new; // new bucket for j
    long i_rank;
    long j_rank; // ranks of nodes
    long j_new_rank;
    priceType rc; // reduced cost of (j, i)
    priceType dr; // rank difference

    numberScan ++;

    i_rank = i->rank();

    // scanning arcs;
    for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {

        ra = a->sister();

        if ( OPEN ( ra ) ) {
            j = a->head();
            j_rank = j->rank();

            if ( j_rank > i_rank ) {
                if ( ( rc = REDUCED_COST( j, i, ra ) ) < 0 ) {
                    j_new_rank = i_rank;
                } else {
                    dr = rc / epsilon;
                    j_new_rank = ( dr < numberOfBucket ) ? i_rank + (long)dr + 1 : numberOfBucket;
                }

                if ( j_rank > j_new_rank ) {
                    j->setRank( j_new_rank);
                    j->setCurrent( ra);

                    if ( j_rank < numberOfBucket ) {
                        b_old = bucketArray + j_rank;
                        REMOVE_FROM_BUCKET( j, b_old );
                    }

                    b_new = bucketArray + j_new_rank;
                    insertToBucket( j, b_new );
                }
            }
        }
    }

    i->decreasePrice( i_rank * epsilon);
    i->setRank( -1);
}

void MinCostSolver::priceUpdate()
{
    register NODE *i;
    excessType remain;
    // total excess of unscanned nodes with positive excess;
    BUCKET *b; // current bucket;
    priceType dp; // amount to be subtracted from prices;

    numberUpdate ++;

    for ( i = nodeArray; i != sentinelNode; i ++ ) {
        if ( i->acessExs() < 0 ) {
            insertToBucket( i, bucketArray );
            i->setRank( 0);
        } else {
            i->setRank( numberOfBucket);
        }
    }

    remain = totalExcess;
    if ( remain < 0.5 ) return;

    // scanning buckets, main loop;
    for ( b = bucketArray; b != lastBucket; b ++ ) {

        while ( isNotEmptyBucket( b) ) {

            GET_FROM_BUCKET( i, b );
            nodeScan( i );

            if ( i ->acessExs() > 0 ) {
                remain -= ( i->acessExs());
                if ( remain <= 0 ) break;
            }
        }
        if ( remain <= 0 ) break;
    }

    if ( remain > 0.5 ) updateFlag = 1;

    // finishup
    // changing prices for nodes which were not scanned during main loop;
    dp = ( b - bucketArray ) * epsilon;

    for ( i = nodeArray; i != sentinelNode; i ++ ) {

        if ( i->rank() >= 0 ) {
            if ( i->rank() < numberOfBucket ) {
                REMOVE_FROM_BUCKET( i, ( bucketArray + i->rank()) );
            }
            if ( i->price() > PricelowBound ) {
                i->decreasePrice( dp);
            }
        }
    }
}

int MinCostSolver::relabel( NODE *i)
{
    register ARC *a; // current arc from i
    register ARC *a_stop; // first arc from the next node
    register ARC *a_max; // arc which provides maximum price
    register priceType p_max; // current maximal price
    register priceType i_price; // price of node  i
    register priceType dp; // current arc partial residual cost

    p_max = PricelowBound;
    i_price = i->price();

    a_max = NULL;

//          1/2 arcs are scanned;
    for ( a = i->current() + 1, a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {
//     寻找pmax最大  max{p(w)-c(v,w)-epsilon}
        if ( OPEN(a) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
            if ( i_price < dp ) {
//                该边可以被push、
                i->setCurrent( a);
                return ( 1);
            }

            p_max = dp;
            a_max = a;
        }
    }

//          2/2 arcs are scanned;
    for ( a = i->first(), a_stop = i->current() + 1; a != a_stop; a ++ ) {
        if ( OPEN( a) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
            if ( i_price < dp ) {
                i->setCurrent( a);
                return ( 1);
            }
            p_max = dp;
            a_max = a;
        }
    }

    // finishup
    if ( p_max != PricelowBound ) {
        i->setPrice( p_max - epsilon);
        i->setCurrent( a_max);
    }
    else { // node can't be relabelled;
        if ( i->suspended() == i->first() ) {
            if ( i->acessExs() == 0 ) {
                i->setPrice( PricelowBound);
            } else {
                if ( nRef == 1 ) {
                    return UNFEASIBLE ;
                } else {
                    return PRICE_OFL ;
                }
            }
        } else { // node can't be relabelled because of suspended arcs;
            priceFlag = 1;
        }
    }

    numberRelabel ++;
    nRel ++;
    return ( 0);
}
//discharge操作，对节点进行push-relabel 直到excess 为0
int MinCostSolver::discharge( NODE *i)
{
    register ARC *a;// an arc from i
    register NODE *j; // head of a
    register long df; // amoumt of flow to be pushed through a
    excessType j_exc; // former excess of j

    numberDiscarge ++;

    a = i->current();
    j = a->head();

    if ( !ADMISSIBLE( i, j, a ) ) {
////////////////////////////////////////////////////////////////////////////////////
//        relabel( i );
        int temp = relabel( i );
        if(temp == UNFEASIBLE || temp == PRICE_OFL)
            return temp;
////////////////////////////////////////////////////////////////////////////////////
        a = i->current();
        j = a->head();
    }

    while ( 1 ) {

        j_exc = j->acessExs();
        if ( j_exc >= 0 ) {

            df = MIN( i->acessExs(), a->acessResidualCapcity() );
            if ( j_exc == 0) numberOfNodeWithExcess++;
            increase_flow( i, j, a, df ); // INCREASE_FLOW
            numberPush ++;

            if ( notInExcessQueue( j ) ) {
                insertToExcessQueue( j );
            }
        }
        else { // j_exc < 0;

            df = MIN( i->acessExs(), a->acessResidualCapcity() );
            increase_flow( i, j, a, df ); // INCREASE_FLOW
            numberPush ++;

            if ( j->acessExs() >= 0 ) {
                if ( j->acessExs() > 0 ) {
                    numberOfNodeWithExcess ++;
//////////////////////////////////////////////////////////////////////
//                    relabel( j );
                    int temp = relabel( j );
                    if(temp == UNFEASIBLE || temp == PRICE_OFL)
                        return temp;
//////////////////////////////////////////////////////////////////
                    insertToExcessQueue( j );
                }
                totalExcess += j_exc;
            }
            else {
                totalExcess -= df;
            }
        }

        if ( i->acessExs() <= 0) numberOfNodeWithExcess --;
        if ( i->acessExs() <= 0 || priceFlag ) break;
//////////////////////////////////////////////////////////////////
//        relabel( i );
        int temp = relabel( i );
        if(temp == UNFEASIBLE || temp == PRICE_OFL)
            return temp;
//////////////////////////////////////////////////////////////////
        a = i->current();
        j = a->head();
    }

    i->setCurrent( a);
    return 0;
}

int MinCostSolver::priceIn()
{
    NODE *i; // current node
    NODE *j;
    ARC *a; // current arc from i
    ARC *a_stop; // first arc from the next node
    ARC *b; // arc to be exchanged with suspended
    ARC *ra; // opposite to a
    ARC *rb; // opposite to b
    priceType rc; // reduced cost
    int n_in_bad; // number of priced_in arcs with negative reduced cost
    int bad_found; // if 1 we are at the second scan if 0 we are at the first scan
    excessType i_exc; // excess of i
    excessType df; // an amount to increase flow


    bad_found = 0;
    n_in_bad = 0;

 restart:

    for ( i = nodeArray; i != sentinelNode; i ++ ) {

        for ( a = i->first() - 1, a_stop = i->suspended() - 1; a != a_stop; a -- ) {

            rc = REDUCED_COST( i, a->head(), a );
            if ( ( rc < 0) && ( a->acessResidualCapcity() > 0) ) { // bad case;
                if ( bad_found == 0 ) {
                    bad_found = 1;
                    updateCutOff();
                    goto restart;
                }
                df = a->acessResidualCapcity();
                increase_flow( i, a->head(), a, df );

                ra = a->sister();
                j  = a->head();

                i->decreaseFirst();
                b = i->first();
                exchange( a, b );

                if ( SUSPENDED( j, ra ) ) {
                    j->decreaseFirst();
                    rb = j->first();
                    exchange( ra, rb );
                }

                n_in_bad ++;
            }
            else {
                if ( ( rc < cutOn ) && ( rc > -cutOn ) ) {
                    i->decreaseFirst();
                    b = i->first();
                    exchange( a, b );
                }
            }
        }
    }


    if ( n_in_bad != 0 ) {

        numberBadPriceIn ++;

        // recalculating excess queue;
        totalExcess = 0;
        numberOfNodeWithExcess = 0;
        rstExcessQueue();

        for ( i = nodeArray; i != sentinelNode; i ++ ) {
            i->setCurrent( i->first());
            i_exc = i->acessExs();
            if ( i_exc > 0 ) { // i is a source;
                totalExcess += i_exc;
                numberOfNodeWithExcess ++;
                insertToExcessQueue( i );
            }
        }

        insertToExcessQueue( AddressOftmpNode );
    }

    if ( timeForPriceIn == TIME_FOR_PRICE_IN2)
        timeForPriceIn = TIME_FOR_PRICE_IN3;
    if ( timeForPriceIn == TIME_FOR_PRICE_IN1)
        timeForPriceIn = TIME_FOR_PRICE_IN2;

    return ( n_in_bad);
}

int MinCostSolver::refine()
{
    NODE *i; // current node
    excessType i_exc; // excess of i
//    long np, nr, ns; // variables for additional print
    int pr_in_int; // current number of updates between price_in

//    np = numberPush;
//    nr = numberRelabel;
//    ns = numberScan;

    numberRefine ++;
//    和nRef 有什么区别。

    nRef ++;
//    第一次 nRef ==1 refine 次数为1
    nRel = 0;
    pr_in_int = 0;

    // initialize;
    totalExcess = 0;
    numberOfNodeWithExcess = 0;
    rstExcessQueue();

    timeForPriceIn = TIME_FOR_PRICE_IN1;

    for ( i = nodeArray; i != sentinelNode; i ++ ) {
        i->setCurrent( i->first());
        i_exc = i->acessExs();
        if ( i_exc > 0 ) { // i  is a source
            totalExcess += i_exc;
            numberOfNodeWithExcess++;
//            加入到要被push的excess 队列中
            insertToExcessQueue( i );
        }
    }

    if ( totalExcess <= 0 ) return 0;

    // (2) main loop

    while ( 1 ) {
//判断当前是否还有溢出的节点
        if ( isEmptyExcessQ() ) {
            if ( nRef > PRICE_OUT_START ) {
                pr_in_int = 0;
                priceIn();
            }

            if ( isEmptyExcessQ() ) break;
        }

        REMOVE_FROM_EXCESS_Q( i );

        // push all excess out of i
        if ( i->acessExs() > 0 ) {
//////////////////////////////////////////////////////////////////
//            discharge( i );
            int temp = discharge( i );
            if(temp == UNFEASIBLE || temp == PRICE_OFL)
                return temp;
//////////////////////////////////////////////////////////////////
            if ( timeForUpdate() || priceFlag ) {
                if ( i->acessExs() > 0 ) {
                    insertToExcessQueue( i );
                }

                if ( priceFlag && ( nRef > PRICE_OUT_START ) ) {
                    pr_in_int = 0;
                    priceIn();
                    priceFlag = 0;
                }

                priceUpdate();

                while ( updateFlag ) {
                    if ( nRef == 1 ) {
                        return UNFEASIBLE ;
                    } else {
                        updateFlag = 0;
                        updateCutOff();
                        numberBadRelabel ++;
                        pr_in_int = 0;
                        priceIn();
                        priceUpdate();
                    }
                }
                nRel = 0;

                if ( nRef > PRICE_OUT_START && (pr_in_int ++ > timeForPriceIn) ) {
                    pr_in_int = 0;
                    priceIn();
                }
            }
        }
    }

    return 0;
}

int MinCostSolver::priceRefine()
{
    NODE *i; // current node
    NODE *j; // opposite node
    NODE *ir; // nodes for passing over the negative cycle
    NODE *is;
    ARC *a; // arc (i,j)
    ARC *a_stop; // first arc from the next node
    ARC *ar;
    long bmax;            // number of farest nonempty bucket
    long i_rank;          // rank of node i
    long j_rank;         // rank of node j
    long j_new_rank;      // new rank of node j
    BUCKET *b;              // current bucket
    BUCKET *b_old;          // old and new buckets of current node
    BUCKET *b_new;
    priceType rc = 0; // reduced cost of a
    priceType dr; // ranks difference
    priceType dp;
    int cc;
    long df; // cycle capacity
    int nnc; // number of negative cycles cancelled during one iteration
    int snc; // total number of negative cycle cancelled

    numberPRefine ++;

    cc = 1;
    snc = 0;

    maxCycleCancel = ( nRef >= START_CYCLE_CANCEL) ? MAX_CYCLES_CANCELLED : 0;
    while ( 1 ) {

        nnc = 0;
        for ( i = nodeArray; i != sentinelNode; i ++ ) {
            i->setRank( 0);
            i->setInp( WHITE);
            i->setCurrent( i->first());
        }
        rstStackQ();

        for ( i = nodeArray; i != sentinelNode; i ++ ) {
            if ( i->inp() == BLACK ) continue;

            i->setBucketNext( NULL);

            // deapth first search
            while ( 1 ) {
                i->setInp( GREY);

                // scanning arcs from node i starting from current
                for ( a = i->current(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        if ( REDUCED_COST ( i, j, a ) < 0 ) {
                            if ( j->inp() == WHITE ) { // fresh node  - step forward
                                i->setCurrent( a);
                                j->setBucketNext( i);
                                i = j;
                                a = j->current();
                                a_stop = (j+1)->suspended();
                                break;
                            }

                            if ( j->inp() == GREY ) { // cycle detected
                                cc = 0;
                                nnc ++;
                                i->setCurrent( a);
                                is = ir = i;
                                df = MAX_32;

                                while ( 1 ) {
                                    ar = ir->current();
                                    if ( ar->acessResidualCapcity() <= df ) {
                                        df = ar->acessResidualCapcity();
                                        is = ir;
                                    }
                                    if ( ir == j ) break;
                                    ir = ir->acessBNext();
                                }

                                ir = i;

                                while ( 1 ) {
                                    ar = ir->current();
                                    increase_flow( ir, ar->head(), ar, df);
                                    if ( ir == j ) break;
                                    ir = ir->acessBNext();
                                }

                                if ( is != i ) {
                                    for ( ir = i; ir != is; ir = ir->acessBNext() ) {
                                        ir->setInp( WHITE);
                                    }
                                    i = is;
                                    a = is->current() + 1;
                                    a_stop = (is+1)->suspended();
                                    break;
                                }
                            }
                        }
                        // if j-color is BLACK - continue search from i
                    }
                } // all arcs from i are scanned

                if ( a == a_stop ) {
                    // step back
                    i->setInp( BLACK);
                    numberPreScan1 ++;
                    j = i->acessBNext();
                    stackQueuePush( i );
                    if ( j == NULL ) break;
                    i = j;
                    i->increaseCurrent();
                }

            } // end of deapth first search
        } // all nodes are scanned


        // () no negative cycle
        // computing longest paths with eps-precision

        snc += nnc;
        if ( snc < maxCycleCancel ) cc = 1;
        if ( cc == 0 ) break;
        bmax = 0;

        while ( isNotEmptyStackQ() ) {

            numberPreScan2 ++;
            STACKQ_POP( i );
            i_rank = i->rank();
            for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

                if ( OPEN( a ) ) {
                    j  = a->head();
                    rc = REDUCED_COST( i, j, a );

                    if ( rc < 0 ) { // admissible arc;
                        dr = (priceType) (( - rc - 0.5 ) / epsilon);
                        if (( j_rank = dr + i_rank ) < numberOfBucket ) {
                            if ( j_rank > j->rank() )
                                j->setRank( j_rank);
                        }
                    }
                }
            } // all arcs from i are scanned

            if ( i_rank > 0 ) {
                if ( i_rank > bmax ) bmax = i_rank;
                b = bucketArray + i_rank;
                insertToBucket( i, b );
            }
        } // end of while-cycle: all nodes are scanned - longest distancess are computed;


        if ( bmax == 0 ) // preflow is eps-optimal;
            { break; }


        for ( b = bucketArray + bmax; b != bucketArray; b -- ) {
            i_rank = b - bucketArray;
            dp = i_rank * epsilon;

            while ( isNotEmptyBucket( b) ) {
                GET_FROM_BUCKET( i, b );
                numberPreScan ++;

                for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        j_rank = j->rank();
                        if ( j_rank < i_rank ) {
                            rc = REDUCED_COST( i, j, a );
                            if ( rc < 0 ) {
                                j_new_rank = i_rank;
                            } else {
                                dr = rc / epsilon;
                                j_new_rank = ( dr < numberOfBucket ) ? i_rank - ( (long)dr + 1 ) : 0;
                            }
                            if ( j_rank < j_new_rank ) {
                                if ( cc == 1 ) {
                                    j->setRank( j_new_rank);
                                    if ( j_rank > 0 ) {
                                        b_old = bucketArray + j_rank;
                                        REMOVE_FROM_BUCKET( j, b_old );
                                    }
                                    b_new = bucketArray + j_new_rank;
                                    insertToBucket( j, b_new );
                                }
                                else {
                                    df = a->acessResidualCapcity();
                                    increase_flow( i, j, a, df );
                                }
                            }
                        }
                    } // end if opened arc
                } // all arcs are scanned

                i->decreasePrice( dp);

            } // end of while-cycle: the bucket is scanned
        } // end of for-cycle: all buckets are scanned

        if ( cc == 0 ) break;

    } // end of main loop
    if ( cc == 0 ) {
        for ( i = nodeArray; i != sentinelNode; i ++) {
            for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                if ( REDUCED_COST( i, a->head(), a ) < - epsilon ) {
                    if ( ( df = a->acessResidualCapcity() ) > 0 ) {
                        increase_flow( i, a->head(), a, df );
                    }
                }
            }
        }
    }

    return ( cc );
}

void MinCostSolver::computePrices()
{
    NODE *i; // current node
    NODE *j; // opposite node
    ARC *a; // arc (i,j)
    ARC *a_stop; // first arc from the next node
    long bmax; // number of farest nonempty bucket
    long i_rank; // rank of node i
    long j_rank; // rank of node j
    long j_new_rank; // new rank of node j
    BUCKET *b; // current bucket
    BUCKET *b_old; // old and new buckets of current node
    BUCKET *b_new;
    priceType rc; // reduced cost of a
    priceType dr; // ranks difference
    priceType dp;
    int cc; // return code: 1 - flow is epsilon optimal 0 - refine is needed

    numberPRefine ++;
    cc = 1;

    // (1) main loop
    // while negative cycle is found or eps-optimal solution is constructed
    while ( 1 ) {

        for ( i = nodeArray; i != sentinelNode; i ++) {
            i->setRank( 0);
            i->setInp( WHITE);
            i->setCurrent( i->first());
        }
        rstStackQ();

        for ( i = nodeArray; i != sentinelNode; i ++ ) {
            if ( i->inp() == BLACK ) continue;

            i->setBucketNext( NULL);
            // depth first search
            while ( 1 ) {
                i->setInp( GREY);

                // scanning arcs from node i
                for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        if ( REDUCED_COST( i, j, a ) < 0 ) {
                            if ( j->inp() == WHITE ) { // fresh node  - step forward
                                i->setCurrent( a);
                                j->setBucketNext( i);
                                i = j;
                                a = j->current();
                                a_stop = (j+1)->suspended();
                                break;
                            }

                            if ( j->inp() == GREY ) { // cycle detected; should not happen
                                cc = 0;
                            }
                        }
                        // if j-color is BLACK - continue search from i
                    }
                } // all arcs from i are scanned

                if ( a == a_stop ) {
                    // step back
                    i->setInp( BLACK);
                    numberPreScan1 ++;
                    j = i->acessBNext();
                    stackQueuePush( i );
                    if ( j == NULL ) break;
                    i = j;
                    i->increaseCurrent();
                }

            } // end of deapth first search
        } // all nodes are scanned
        if ( cc == 0 ) break;
        bmax = 0;

        while ( isNotEmptyStackQ() ) {
            numberPreScan2 ++;
            STACKQ_POP( i );
            i_rank = i->rank();
            for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                if ( OPEN( a ) ) {
                    j  = a->head();
                    rc = REDUCED_COST( i, j, a );


                    if ( rc < 0 ) {// admissible arc
                        dr = - rc;
                        if (( j_rank = dr + i_rank ) < numberOfBucket ) {
                            if ( j_rank > j->rank() )
                                j->setRank( j_rank);
                        }
                    }
                }
            } // all arcs from i are scanned

            if ( i_rank > 0 ) {
                if ( i_rank > bmax ) bmax = i_rank;
                b = bucketArray + i_rank;
                insertToBucket( i, b );
            }
        } // end of while-cycle: all nodes are scanned - longest distancess are computed;

        if ( bmax == 0 )
            { break; }

        for ( b = bucketArray + bmax; b != bucketArray; b -- ) {
            i_rank = b - bucketArray;
            dp = i_rank;

            while ( isNotEmptyBucket( b) ) {
                GET_FROM_BUCKET( i, b );
                numberPreScan ++;

                for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        j_rank = j->rank();
                        if ( j_rank < i_rank ) {
                            rc = REDUCED_COST( i, j, a );

                            if ( rc < 0 ) {
                                j_new_rank = i_rank;
                            } else {
                                dr = rc;
                                j_new_rank = ( dr < numberOfBucket ) ? i_rank - ( (long)dr + 1 ) : 0;
                            }
                            if ( j_rank < j_new_rank ) {
                                if ( cc == 1 ) {
                                    j->setRank( j_new_rank);
                                    if ( j_rank > 0 ) {
                                        b_old = bucketArray + j_rank;
                                        REMOVE_FROM_BUCKET( j, b_old );
                                    }
                                    b_new = bucketArray + j_new_rank;
                                    insertToBucket( j, b_new );
                                }
                            }
                        }
                    } // end if opened arc
                } // all arcs are scanned

                i->decreasePrice( dp);

            } // end of while-cycle: the bucket is scanned
        } // end of for-cycle: all buckets are scanned

        if ( cc == 0 ) break;

    } // end of main loop
}

void MinCostSolver::priceOut()
{
    NODE *i; // current node
    ARC *a; // current arc from i
    ARC *a_stop; // first arc from the next node
    ARC *b; // arc to be exchanged with suspended
    double numbelCutOff; // -cut_off
    double rc; // reduced cost

    numbelCutOff = - cutOff;

    for ( i = nodeArray; i != sentinelNode; i ++) {
        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

            rc = REDUCED_COST( i, a->head(), a );
            if ( ( rc > cutOff && CLOSED(a->sister()) ) ||
                 ( rc < numbelCutOff && CLOSED(a) ) ) { // suspend the arc

                b = i->first();
                i->increaseFirst();
                exchange( a, b );
            }
        }
    }
}

int MinCostSolver::updateEpsilon()
{
    // decrease epsilon after epsilon-optimal flow is constructed;
    if ( epsilon <= 1 ) return ( 1 );

    epsilon = (priceType) (ceil ( (double) epsilon / scaleFactor ));
    cutOff = cutOffFactor * epsilon;
    cutOn = cutOff * CUT_OFF_GAP;

    return ( 0 );
}

int MinCostSolver::checkFeasible()
{
    if ( checkSolution == false)
        return ( 0);

    NODE *i;
    ARC *a, *a_stop;
    long fa;
    int ans = 1;

    for ( i = nodeArray; i != sentinelNode; i ++) {
        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
            if ( capArray[ ARC_NUMBER(a) ] > 0) {
                fa = capArray[ ARC_NUMBER(a) ] - a->acessResidualCapcity();
                if ( fa < 0) {
                    ans = 0;
                    break;
                }
                nodeBalance[ i - nodeArray ] -= fa;
                nodeBalance[ a->head() - nodeArray ] += fa;
            }
        }
    }

    for ( i = nodeArray; i != sentinelNode; i ++) {
        if ( nodeBalance[ i - nodeArray ] != 0) {
            ans = 0;
            break;
        }
    }

    return ( ans);
}

int MinCostSolver::checkCs()
{
    // check complimentary slackness;
    NODE *i;
    ARC *a, *a_stop;

    for ( i = nodeArray; i != sentinelNode; i ++) {
        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

            if ( OPEN(a) && (REDUCED_COST(i, a->head(), a) < 0) ) {
                return ( 0);
            }
        }
    }
    return(1);
}

int MinCostSolver::checkEpsOpt()
{
    NODE *i;
    ARC *a, *a_stop;

    for ( i = nodeArray; i != sentinelNode; i ++) {
        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

            if ( OPEN(a) && (REDUCED_COST(i, a->head(), a) < - epsilon) ) {
                return ( 0);
            }
        }
    }
    return(1);
}

void MinCostSolver::initSolution()
{
    ARC *a; // current arc (i,j)
    NODE *i; // tail of a
    NODE *j; // head of a
    long df; // residual capacity

    for ( a = arcArray; a != sentinelArc; a ++ ) {
        if ( a->acessResidualCapcity() > 0 && a->cost() < 0 ) {
            df = a->acessResidualCapcity();
            i  = a->sister()->head();
            j  = a->head();
            increase_flow( i, j, a, df );
        }
    }
}

void MinCostSolver::csCostReinit()
{
    if ( costRestart == false)
        return;

    NODE *i; // current node
    ARC *a;          // current arc
    ARC *a_stop;
    BUCKET *b; // current bucket
    priceType rc, minc, sum;


    for ( b = bucketArray; b != lastBucket; b ++) {
        rstBucket( b);
    }

    rc = 0;
    for ( i = nodeArray; i != sentinelNode; i ++) {
        rc = MIN(rc, i->price());
        i->setFirst( i->suspended());
        i->setCurrent( i->first());
        i->setQueueNext( sentinelNode);
    }

    // make prices nonnegative and multiply
    for ( i = nodeArray; i != sentinelNode; i ++) {
        i->setPrice( (i->price() - rc) * _dn);
    }

    // multiply arc costs
    for (a = arcArray; a != sentinelArc; a ++) {
        a->multiplyCost( _dn);
    }

    sum = 0;
    for ( i = nodeArray; i != sentinelNode; i ++) {
        minc = 0;
        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
            if ( (OPEN(a) && ((rc = REDUCED_COST(i, a->head(), a)) < 0)) )
                minc = MAX( epsilon, -rc);
        }
        sum += minc;
    }

    epsilon = ceil(sum / _dn);

    cutOffFactor = CUT_OFF_COEF * pow((double)_n, CUT_OFF_POWER);

    cutOffFactor = MAX( cutOffFactor, CUT_OFF_MIN);

    nRef = 0;

    numberRefine = numberDiscarge = numberPush = numberRelabel = 0;
    numberUpdate = numberScan = numberPRefine = numberPreScan = numberPreScan1 =
        numberBadPriceIn = numberBadRelabel = 0;

    priceFlag = 0;

    excessQueueFirst = NULL;
}

int MinCostSolver::cs2CostRestart( double *objective_cost)
{
    // restart after a cost update;
    if ( costRestart == false)
        return 0;

    int cc; // for storing return code;

    printf("c \nc ******************************\n");
    printf("c Restarting after a cost update\n");
    printf("c ******************************\nc\n");

    csCostReinit();

    printf ("c Init. epsilon = %6.0f\n", epsilon);
    cc = updateEpsilon();

    if (cc != 0) {
        printf("c Old solution is optimal\n");
    }
    else {
        do { // scaling loop
            while ( 1 ) {
                if ( ! priceRefine() )
                    break;

                if ( nRef >= PRICE_OUT_START ) {
                    if ( priceIn() )
                        break;
                }
                if ((cc = updateEpsilon ()))
                    break;
            }
            if (cc) break;
//////////////////////////////////////////////
//            refine();
            int temp = refine();
            if(temp == UNFEASIBLE || temp == PRICE_OFL)
                return temp;
//////////////////////////////////////////
            if ( nRef >= PRICE_OUT_START ) {
                priceOut();
            }
            if ( updateEpsilon() )
                break;
        } while ( cc == 0 );
    }

    finishup( objective_cost );
     return 0;
}

void MinCostSolver::printSolution()
{


    NODE *i;
    ARC *a;
    long ni;

    for ( i = nodeArray; i < nodeArray + _n; i ++ ) {
        ni = NODE_NUMBER( i );
        for ( a = i->suspended(); a != (i+1)->suspended(); a ++) {
            if ( capArray[ ARC_NUMBER (a) ]  > 0 ) { // 不要sister边。
                if(capArray[ ARC_NUMBER(a) ] - a->acessResidualCapcity() > 0)
                    // 把具有流量的边push到 usedARC容器里。
                    usedARC.push_back({NODE_NUMBER(i), NODE_NUMBER(a->head()),capArray[ ARC_NUMBER(a) ] - a->acessResidualCapcity(),a->cost()});

            }
        }
    }
}

void MinCostSolver::printGraph()
{
    NODE *i;
    ARC *a;
    long ni, na;
    printf ("\nGraph: %d\n", _n);
    for ( i = nodeArray; i < nodeArray + _n; i ++ ) {
        ni = NODE_NUMBER( i );
        printf("\nNode %d", ni);
        for ( a = i->suspended(); a != (i+1)->suspended(); a ++) {
            na = ARC_NUMBER( a );
            printf("\n {%d} %d -> %d  cap: %d  cost: %d", na,
                ni, NODE_NUMBER(a->head()), capArray[ARC_NUMBER(a)], a->cost());
        }
    }
}

void MinCostSolver::finishup( double *objective_cost)
{
    ARC *a; // current arc
    long na; // corresponding position in capacity array
    double obj_internal = 0; // objective
    priceType cs; // actual arc cost
    long flow; // flow through an arc
    NODE *i;

    // (1) NO_ZERO_CYCLES?
    if ( numberZeroCycles == true) {
        for ( a = arcArray; a != sentinelArc; a ++ ) {
            if ( a->cost() == 1) {
                assert( a->sister()->cost() == -1);
                a->setCost( 0);
                a->sister()->setCost( 0);
            }
        }
    }

    // (2)
    for ( a = arcArray, na = 0; a != sentinelArc ; a ++, na ++ ) {
        cs = a->cost() / _dn;
        if ( capArray[na]  > 0 && (flow = capArray[na] - a->acessResidualCapcity()) != 0 )
            obj_internal += (double) cs * (double) flow;
        a->setCost( cs);
    }

    for ( i = nodeArray; i != sentinelNode; i ++) {
        i->setPrice( (i->price() / _dn));
    }

    // (3) COMP_DUALS?
    if ( computeDuals == true) {
        computePrices();
    }

    *objective_cost = obj_internal;
}

int MinCostSolver::cs2( double *objective_cost)
{
    // the main calling function;
    int cc = 0; // for storing return code;


    // (1) update epsilon first;
    updateEpsilon();


    // (2) scaling loop;
    do {
//////////////////////////////////////////////
//            refine();
             int temp = refine();
             if(temp == UNFEASIBLE || temp == PRICE_OFL)
                return temp;
//////////////////////////////////////////

        if ( nRef >= PRICE_OUT_START )
            priceOut();

        if ( updateEpsilon() )
            break;

        while (1) {
            if ( ! priceRefine() )
                break;

            if ( nRef >= PRICE_OUT_START ) {
                if ( priceIn() ) break;
                if ( (cc = updateEpsilon()) ) break;
            }
        }
    } while ( cc == 0 );


    // (3) finishup;
    finishup( objective_cost );
    return 0;
}
int MinCostSolver::runSolition()
{
    double objective_cost;


    // (4) ordering, etc.;
    preProscessing();


    // () CHECK_SOLUTION?
    if ( checkSolution == true) {
        nodeBalance = (long long int *) calloc (_n+1, sizeof(long long int));
        for ( NODE *i = nodeArray; i < nodeArray + _n; i ++ ) {
            nodeBalance[i - nodeArray] = i->acessExs();
        }
    }


    // (5) initializations;
    _m = 2 * _m;
    if(initialization() == ALLOCATION_FAULT)
        return ALLOCATION_FAULT; // works already with 2*m;

    int temp = cs2( &objective_cost );
    if(temp == UNFEASIBLE || temp == PRICE_OFL)
       return temp;

    total_cost = objective_cost;

    if ( checkSolution == true ) {
        printf("c checking feasibility...\n");
        if ( checkFeasible() )
            printf("c ...OK\n");
        else
            printf("c ERROR: solution infeasible\n");
        printf("c computing prices and checking CS...\n");
        computePrices();
        if ( checkCs() )
            printf("c ...OK\n");
        else
            printf("ERROR: CS violation\n");
    }

    // () PRINT_ANS?
    if ( printAnswer == true ) {
        printSolution();
    }

    // () cleanup;
//    deleteMemory();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
//find Server
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FindServer::faultHandeler( int cc)
{
    // abnormal finish
    printf ("\nError %d\n", cc );
    exit( cc);
}

void FindServer::allocateMemory()
{
    // (1) allocate memory for 'nodes', 'arcs' and internal arrays;
//  后面需要设置一个  Sentinel Node  另一个是为了跳过 初始0点
    nodeArray = (NODE*) calloc ( _n+2,   sizeof(NODE) );
//    边数乘2   1是以为设置 边 Sentinel arc
    arcArray = (ARC*)  calloc ( 2*_m+1, sizeof(ARC) );
//    同样乘二 保存了当前边和sister 边的 capacity
    capArray = (long*) calloc ( 2*_m,   sizeof(long) );
//保存了sister 和当前边的tail node
    arcTail = (long*) calloc ( 2*_m,   sizeof(long) );
//    当前arc 所处的位置索引， 访问 capArray 和 _arc数组
    arcFirst = (long*) calloc ( _n+2,   sizeof(long) );
    // arcfirstOutArc [ 0 .. n+1 ] = 0 - initialized by calloc;
//设置所有的 excess 为0
    for ( NODE *in = nodeArray; in <= nodeArray + _n; in ++ ) {
        in->setExcess( 0);
    }
    if ( nodeArray == NULL || arcArray == NULL || arcFirst == NULL || arcTail == NULL) {
        printf("Error:  Memory allocation problem inside CS2\n");
        exit( 1);
    }

    // (2) resets;
    posCurrent = 0;
    arcCurrent = arcArray; // set "current" pointer to the first arc *_arc 是初始位置
//    这两个先设置方向最大，一会改过来
    nodeMax = 0;
    nodeMin = _n;
    maxCost = 0;
    _total_supply = _total_demand = 0;
    // at this moment we are ready to add arcs and build the network,
    // by using setArc()...
}

void FindServer::deleteMemory()
{
    if ( arcArray) free ( arcArray );
    if ( dNode) delete dNode;
    if ( capArray) free ( capArray );
    if ( bucketArray) free ( bucketArray );
    if ( checkSolution == true) free ( nodeBalance );
    if ( nodeArray) {
        nodeArray = nodeArray - nodeMin;
        free ( nodeArray );
    }
}

void FindServer::setArc( long tail_node_id, long head_node_id,
                        long low_bound, long up_bound, // up_bound is basically capacity;
                        priceType cost)
{
    // DIMACS format:
    // c arc has <tail> <head> <capacity l.b.> <capacity u.b> <cost>

    if ( tail_node_id < 0 || tail_node_id > _n ||
         head_node_id < 0 || head_node_id > _n ) {
        printf("Error:  Arc with head or tail out of bounds inside CS2\n");
        exit( 1);
    }
    if ( up_bound < 0 ) {
        up_bound = MAX_32;
        printf("Warning:  Infinite capacity replaced by BIGGEST_FLOW\n");
    }
    if ( low_bound < 0 || low_bound > up_bound ) {
        printf("Error:  Wrong capacity bounds inside CS2\n");
        exit( 1);
    }

    arcFirst[tail_node_id + 1] ++;
    arcFirst[head_node_id + 1] ++;
    iNode = nodeArray + tail_node_id;
    jNode = nodeArray + head_node_id;

    // store information about the arc   sister 的 tail 和当前tail 一起存储
    arcTail[posCurrent]   = tail_node_id;
    arcTail[posCurrent+1] = head_node_id;
//
    arcCurrent->setHead( jNode );
    arcCurrent->setResidualCapacity( up_bound - low_bound );
    capArray[posCurrent] = up_bound;
    arcCurrent->setCost( cost );
    arcCurrent->setSisterArc( arcCurrent + 1 );
//    设置head 和tail的int 型
    arcCurrent->headInt = head_node_id;
    arcCurrent->tailInt = tail_node_id;
//    设置 sister边（残存反向边）
    ( arcCurrent + 1 )->setHead( nodeArray + tail_node_id );
    ( arcCurrent + 1 )->setResidualCapacity( 0 );
    capArray[posCurrent+1] = 0;
//    sister边的cost和原cost 相反
    ( arcCurrent + 1 )->setCost( -cost );
//    两者互为sister
    ( arcCurrent + 1 )->setSisterArc( arcCurrent );

    iNode->decreaseExcess( low_bound );
    jNode->increaseExcess( low_bound );

    // searching for minimum and maximum node
    if ( head_node_id < nodeMin ) nodeMin = head_node_id;
    if ( tail_node_id < nodeMin ) nodeMin = tail_node_id;
    if ( head_node_id > nodeMax ) nodeMax = head_node_id;
    if ( tail_node_id > nodeMax ) nodeMax = tail_node_id;

    if ( cost < 0 ) cost = -cost;
    if ( (cost > maxCost && up_bound > 0) && (head_node_id != nodeNum+1))
        maxCost = cost;

    // prepare for next arc to be added;  每次添加都添加当前边和sister的边
    arcCurrent += 2;
    posCurrent += 2;
}

void FindServer::setDemandOfNode( long id, long excess)
{
    // set supply and demand of nodes; not used for transhipment nodes;
    if ( id < 0 || id > _n ) {
        printf("Error:  Unbalanced problem inside CS2\n");
        exit( 1);
    }
//    可以设置负的excess
    (nodeArray + id)->setExcess( excess);
    if ( excess > 0) _total_supply += excess;
    if ( excess < 0) _total_demand -= excess;
}

void FindServer::preProscessing()
{
    // called after the arcs were just added and before runSolition();
    // ordering arcs - linear time algorithm;
    long i;
    long last, arc_num, arc_new_num;;
    long tail_node_id;
    NODE *head_p;
    ARC *arc_new, *arc_tmp;
    long up_bound;
    priceType cost; // arc cost;
    excessType cap_out; // sum of outgoing capacities
    excessType cap_in; // sum of incoming capacities

    if ( ABS( _total_supply - _total_demand ) > 0.5 ) {
        printf("Error:  Unbalanced problem inside CS2\n");
        exit( 1);
    }

    // first arc from the first node 0
    ( nodeArray + nodeMin )->setFirst( arcArray );



    for ( i = nodeMin + 1; i <= nodeMax + 1; i ++ ) {
//        ？ 为啥要加起来啊！
        arcFirst[i] += arcFirst[i-1];
        ( nodeArray + i )->setFirst( arcArray + arcFirst[i] );
    }

    // scanning all the nodes except the last
    for ( i = nodeMin; i < nodeMax; i ++ )
    {
//        first - arc基地址？
        last = ( ( nodeArray + i + 1 )->first() ) - arcArray;

        for ( arc_num = arcFirst[i]; arc_num < last; arc_num ++ ) {
            tail_node_id = arcTail[arc_num];

            while ( tail_node_id != i ) {
                arc_new_num = arcFirst[tail_node_id];
                arcCurrent = arcArray + arc_num;
                arc_new = arcArray + arc_new_num;

                // arc_current must be cited in the position arc_new
                // swapping these arcs:

                head_p = arc_new->head();
                arc_new->setHead( arcCurrent->head() );
                arcCurrent->setHead( head_p );

                up_bound          = capArray[arc_new_num];
                capArray[arc_new_num] = capArray[arc_num];
                capArray[arc_num]     = up_bound;

                up_bound = arc_new->acessResidualCapcity();
                arc_new->setResidualCapacity( arcCurrent->acessResidualCapcity() );
                arcCurrent->setResidualCapacity( up_bound) ;

                cost = arc_new->cost();
                arc_new->setCost( arcCurrent->cost() );
                arcCurrent->setCost( cost );

                if ( arc_new != arcCurrent->sister() ) {
                    arc_tmp = arc_new->sister();
                    arc_new->setSisterArc( arcCurrent->sister() );
                    arcCurrent->setSisterArc( arc_tmp );

                    arcCurrent->sister()->setSisterArc( arcCurrent );
                    arc_new->sister()->setSisterArc( arc_new );
                }

                arcTail[arc_num] = arcTail[arc_new_num];
                arcTail[arc_new_num] = tail_node_id;

                // we increase arcfirstOutArc[tail_node_id]
                arcFirst[tail_node_id] ++ ;

                tail_node_id = arcTail[arc_num];
            }
        }
        // all arcs outgoing from  i  are in place
    }
    // arcs are ordered by now!


    // testing network for possible excess overflow
    for ( NODE *ndp = nodeArray + nodeMin; ndp <= nodeArray + nodeMax; ndp ++ ) {
        cap_in  =   ( ndp->acessExs() );
        cap_out = - ( ndp->acessExs() );
        for ( arcCurrent = ndp->first(); arcCurrent != (ndp+1)->first();
              arcCurrent ++ ) {
            arc_num = arcCurrent - arcArray;
            if ( capArray[arc_num] > 0 ) cap_out += capArray[arc_num];
            if ( capArray[arc_num] == 0 )
                cap_in += capArray[ arcCurrent->sister() - arcArray ];
        }
    }
    if ( nodeMin < 0 || nodeMin > 1 ) {
        printf("Error:  Node ids must start from 0 or 1 inside CS2\n");
        exit( 1);
    }

    // adjustments due to nodes' ids being between nodeMin - nodeMax;
    _n = nodeMax - nodeMin + 1;
//    重定义指针
    nodeArray = nodeArray + nodeMin;

    // () free internal memory, not needed anymore inside CS2;
    free ( arcFirst );
    free ( arcTail );
}

int FindServer:: initialization()
{
    // initialization;
    // called after allocateMemory() and all nodes and arcs have been inputed;

    NODE *i; // current node
    ARC *a; // current arc
    ARC *a_stop;
    BUCKET *b; // current bucket
    long df;

    scaleFactor = (long) SCALE_DEFAULT;
//    设置末尾节点为 哨兵节点
    sentinelNode = nodeArray + _n;
    sentinelArc  = arcArray + _m;
//对所有节点进行初始化操作
    for ( i = nodeArray; i != sentinelNode; i ++ ) {
//        设置价格为0
        i->setPrice( 0);
//        设置 到超级汇点的流量，averagetotal
        i->flewToSuperSink = 0;
        //i->averageExcessCost = MAX_64;
//       suspend 只有在这里赋值了，所以它永远等于first'
        i->setSuspended( i->first());
//        push 队列里放 放 哨兵
        i->setQueueNext( sentinelNode);

    }

    sentinelNode->setFirst( sentinelArc);
    sentinelNode->setSuspended( sentinelArc);

    // saturate negative arcs, e.g. in the circulation problem case
    int countiii = 0;
    for ( i = nodeArray; i != sentinelNode; i ++ ) {
        countiii++;
//        if(countiii>=1199)
//            cout << "wait";
        i->totalCost = 0;
        i->averageExcessCost = MAX_32;
        i->serverId = 0x7FFF;
//给消费节点进行各个 处理变量的初始化
        if(i->acessExs()>0)
        {
           i->serverId = fitServer(i->acessExs());
           if(i->serverId > 0)
           {
            int serverCost = serverCostVect[i->serverId].cost;
//               有满足要求的服务器
            i->averageExcessCost = (serverCost*(_n+1) + nodeBaseCost[NODE_NUMBER(i)]*(_n+1))/i->acessExs();
            averageExcessCostVecmoniRaw.push_back({NODE_NUMBER(i),i->averageExcessCost,nodeBaseCost[NODE_NUMBER(i)],i->serverId,0,0,i->acessExs(),0});
            i->totalCostMax = serverCost * (_n+1);
            i->excessMax = i->acessExs();
           }
         else
         {
//             没有满足要求的服务器
                averageExcessCostVecmoniRaw.push_back({NODE_NUMBER(i),MAX_32,nodeBaseCost[NODE_NUMBER(i)],i->serverId,0,0,1,0});
                i->excessMax = 1;
                i->totalCostMax = MAX_32;
         }
        }
        else
        {
            averageExcessCostVecmoniRaw.push_back({NODE_NUMBER(i),MAX_32,nodeBaseCost[NODE_NUMBER(i)],i->serverId,0,0,1,0});
            i->excessMax = 1;
            i->totalCostMax = MAX_32;
        }
//        FOR_ACTIVEarcArray_a_FROM_i
        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {
            if ( a->cost() < 0) {
//                初始化的时候所有新建的sister负边都没有residual
                if ( ( df = a->acessResidualCapcity()) > 0) {
                    increase_flow( i, a->head(), a, df);
                }
            }
        }
    }
//  节点数量+1
    _dn = _n + 1;
//    if ( numberZeroCycles == true) { // NO_ZERO_CYCLES
//        _dn = 2 * _dn;
//    }
//  遍历所有的arc
    for ( a = arcArray; a != sentinelArc; a ++ ) {
//  这个乘是干啥用的，把所有边的cost 增大
        a->multiplyCost( _dn);
    }



    if ((double) maxCost * (double) _dn > MAX_64) {
        printf("Warning:  Arc lengths too large, overflow possible\n");
    }
    multiMaxCost = maxCost * _dn;
    numberOfBucket = (long) (_dn * ceil(scaleFactor) + 2);

    bucketArray = (BUCKET*) calloc ( numberOfBucket, sizeof(BUCKET));
    if ( bucketArray == NULL )
        return ALLOCATION_FAULT;
//得到最后一个bucket的地址
    lastBucket = bucketArray + numberOfBucket;

    dNode = new NODE; // used as reference;
//init bucket
    for ( b = bucketArray; b != lastBucket; b ++ ) {
        rstBucket( b);
    }

    epsilon = multiMaxCost;
    if ( epsilon < 1) {
        epsilon = 1;
    }

    PricelowBound = -PRICEUP_UP_BOUND;

    cutOffFactor = CUT_OFF_COEF * pow( (double)_n, CUT_OFF_POWER);

    cutOffFactor = MAX( cutOffFactor, CUT_OFF_MIN);

    nRef = 0;

    priceFlag = 0;

    AddressOftmpNode = &tempUseNode;

    excessQueueFirst = NULL;

    return 0;
    //printGraph(); // debug;
}

void FindServer::nodeScan( NODE *i)
{
    NODE *j; // opposite node
    ARC *a; // (i, j)
    ARC *a_stop; // first arc from the next node
    ARC *ra; // (j, i)
    BUCKET *b_old; // old bucket contained j
    BUCKET *b_new; // new bucket for j
    long i_rank;
    long j_rank; // ranks of nodes
    long j_new_rank;
    priceType rc; // reduced cost of (j, i)
    priceType dr; // rank difference

    numberScan ++;

    i_rank = i->rank();

    // scanning arcs;
    for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {

        ra = a->sister();

        if ( OPEN ( ra ) ) {
            j = a->head();
            j_rank = j->rank();

            if ( j_rank > i_rank ) {
                if ( ( rc = REDUCED_COST( j, i, ra ) ) < 0 ) {
                    j_new_rank = i_rank;
                } else {
                    dr = rc / epsilon;
                    j_new_rank = ( dr < numberOfBucket ) ? i_rank + (long)dr + 1 : numberOfBucket;
                }

                if ( j_rank > j_new_rank ) {
                    j->setRank( j_new_rank);
                    j->setCurrent( ra);

                    if ( j_rank < numberOfBucket ) {
                        b_old = bucketArray + j_rank;
                        REMOVE_FROM_BUCKET( j, b_old );
                    }

                    b_new = bucketArray + j_new_rank;
                    insertToBucket( j, b_new );
                }
            }
        }
    }

    i->decreasePrice( i_rank * epsilon);
    i->setRank( -1);
}

void FindServer::priceUpdate()
{
    register NODE *i;
    excessType remain;
    // total excess of unscanned nodes with positive excess;
    BUCKET *b; // current bucket;
    priceType dp; // amount to be subtracted from prices;

    numberUpdate ++;
// 把excess 大于零和小于零的分开
    for ( i = nodeArray; i != sentinelNode; i ++ ) {
        if ( i->acessExs() < 0 ) {
            insertToBucket( i, bucketArray );
            i->setRank( 0); // 在buket
        } else {
            i->setRank( numberOfBucket);// 不在bucket
        }
    }

    remain = totalExcess;
    // 剩下的excess还有多少
    if ( remain < 0.5 ) return;

    // scanning buckets, main loop;
    for ( b = bucketArray; b != lastBucket; b ++ ) {

        while ( isNotEmptyBucket( b) ) {

            GET_FROM_BUCKET( i, b );
            nodeScan( i );

            if ( i ->acessExs() > 0 ) {
                remain -= ( i->acessExs());
                if ( remain <= 0 ) break;
            }
        }
        if ( remain <= 0 ) break;
    }

    if ( remain > 0.5 ) updateFlag = 1;

    // finishup
    // changing prices for nodes which were not scanned during main loop;
    dp = ( b - bucketArray ) * epsilon;

    for ( i = nodeArray; i != sentinelNode; i ++ ) {

        if ( i->rank() >= 0 ) {
            if ( i->rank() < numberOfBucket ) {
                REMOVE_FROM_BUCKET( i, ( bucketArray + i->rank()) );
            }
            if ( i->price() > PricelowBound ) {
                i->decreasePrice( dp);
            }
        }
    }
}
// 在relabel里决定是否去当前节点是否去超级汇点。
int FindServer::relabel( NODE *i)
{
     ARC *a; // current arc from i
     ARC *a_stop; // first arc from the next node
     ARC *a_max; // arc which provides maximum price
     priceType p_max; // current maximal price
     priceType i_price; // price of node  i
     priceType dp; // current arc partial residual cost
    int flowToSuperNode = 0;// 用来判断relabel 之后是不是去了超级汇点
    p_max = PricelowBound;
    i_price = i->price();

    a_max = NULL;
//current to last arc are scaned
//          1/2 arcs are scanned;
    for ( a = i->current() + 1, a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {
//     寻找pmax最大  max{p(w)-c(v,w)-epsilon}
//        判断是否是超级汇点的边,如果是，则改变这条边的cost为 servercost/excess


        if((a->head() - nodeArray) == nodeNum+1)
            if(findInVect(serverNode,a->head() - nodeArray))
                a->_cost = ceil((double)serverCostMax * (_n+1)/(i->acessExs()+i->flewToSuperSink));
            else
                a->_cost = ceil((double)serverCostMax * (_n+1)/i->acessExs());

//找 p(w)-c(v,w) 最大
        if ( OPEN(a) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
//            再判断最大的是不是去 超级汇点的边
            if( (a->head() - nodeArray) == nodeNum+1)
                flowToSuperNode = 1;
            else
                flowToSuperNode = 0;
            if ( i_price < dp ) {
//                该边可以被push 出去
//                如果去了超级汇点，则服务器数量加1。
                if(flowToSuperNode == 1)
                {
                    serverCount ++;
                    vector<int>::iterator it;
                    it = find(serverNode.begin(),serverNode.end(),a->head() - nodeArray);
                    if(it == serverNode.end())
                    serverNode.push_back(i- nodeArray);
                    i->flewToSuperSink += i->acessExs();
                }
                i->setCurrent(a);
                return ( 1);
            }

            p_max = dp;  // 最大可以设置的price
            a_max = a;
        }
    }

//         first to current arc are scanned;
    for ( a = i->first(), a_stop = i->current() + 1; a != a_stop; a ++ ) {
 //     寻找pmax最大  max{p(w)-c(v,w)-epsilon}
  //        判断是否是超级汇点的边,如果是，则改变这条边的cost为 servercost/excess
        if((a->head() - nodeArray) == nodeNum+1)
            a->_cost = ceil((double)serverCostMax* (_n+1)/i->acessExs());

        if ( OPEN( a) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
 //            再判断最大的是不是去超级汇点的边
 //            再判断最大的是不是去 超级汇点的边
             if( (a->head() - nodeArray) == nodeNum+1)
                flowToSuperNode = 1;
             else
                flowToSuperNode = 0;
            if ( i_price < dp ) {
//                该边可以被push 出去
//                如果去了超级汇点，则服务器数量加1。
               if(flowToSuperNode == 1)
               {
                   serverCount ++;
                   vector<int>::iterator it;
                   it = find(serverNode.begin(),serverNode.end(),a->head() - nodeArray);
                   if(it == serverNode.end())
                   serverNode.push_back(i - nodeArray);
                   i->flewToSuperSink += i->acessExs();
               }
                i->setCurrent(a);
                return ( 1);
            }
            p_max = dp;
            a_max = a;
        }
    }

    // finishup
    if ( p_max != PricelowBound ) {
//        更新price： relabel
        i->setPrice( p_max - epsilon);
        i->setCurrent( a_max);
    }
    else { // node can't be relabelled;  已经到了最低点了
        if ( i->suspended() == i->first() ) {
            if ( i->acessExs() == 0 ) {
                i->setPrice( PricelowBound);
            } else {
                if ( nRef == 1 ) {
                    return UNFEASIBLE ;
                } else {
                    return PRICE_OFL ;
                }
            }
        } else { // node can't be relabelled because of suspended arcs;
            priceFlag = 1;
        }
    }
    if(flowToSuperNode == 1)
    {
        serverCount ++;
        vector<int>::iterator it;
        it = find(serverNode.begin(),serverNode.end(),a->head() - nodeArray);
        if(it == serverNode.end())
        serverNode.push_back(i - nodeArray);
        i->flewToSuperSink += i->acessExs();
    }
    numberRelabel ++;
    nRel ++;
    return ( 0);
}
//discharge操作，对节点进行push-relabel 直到excess 为0
int FindServer::discharge( NODE *i)
{
    register ARC *a;// an arc from i
    register NODE *j; // head of a
    register long df; // amoumt of flow to be pushed through a
    excessType j_exc; // former excess of j

    numberDiscarge ++;

    a = i->current();
    j = a->head();
    if(i-nodeArray == 43)
        j_exc = 0;
//    在对节点进行操作的时候，直接对
    if((a->head() - nodeArray) == nodeNum+1)
        a->_cost = ceil((double)serverCostMax* (_n+1)/i->acessExs());
//  当前的要被push的边，是不是admissible 如果是，则直接push，不是，则进行relabel操作，relabel操作 包括选择admissible的边或者relabel
    if ( !ADMISSIBLE( i, j, a ) ) {
////////////////////////////////////////////////////////////////////////////////////
//        relabel( i );
        int temp = relabel( i );
        if(temp == UNFEASIBLE || temp == PRICE_OFL)
            return temp;
////////////////////////////////////////////////////////////////////////////////////
//     重新定位要进行push操作的current边
        a = i->current();
        j = a->head();
    }
    else if(j-nodeArray == nodeNum+1)//  admissble 并且这条边时通向超级汇点的
    {
//        更新cost
        a->_cost = ceil((double)serverCostMax * (_n+1)/(i->acessExs()+i->flewToSuperSink));
//        推掉所有的excess
        increase_flow( i, j, a, i->acessExs() );
        numberOfNodeWithExcess --;
        numberPush ++;
        total_cost -= i->acessExs();
        return 0;
    }
    while ( 1 ) {

        j_exc = j->acessExs();
//        目的节点j的excess大于零
        if ( j_exc >= 0 ) {

            df = MIN( i->acessExs(), a->acessResidualCapcity() );
            if ( j_exc == 0) numberOfNodeWithExcess++;
            increase_flow( i, j, a, df ); // INCREASE_FLOW
            numberPush ++;
//           初始化的时候把Node的qnext都设置为 sentinel_node
            if ( notInExcessQueue( j ) ) {
                insertToExcessQueue( j );
            }
        }
//         如果目的节点excess小于0 则要把目的节点 加到push队列中
        else { // j_exc < 0;

            df = MIN( i->acessExs(), a->acessResidualCapcity() );
            increase_flow( i, j, a, df ); // INCREASE_FLOW
            numberPush ++;

            if ( j->acessExs() >= 0 ) {
//                Push完之后 输出点的excess 超了
                if ( j->acessExs() > 0 ) {
                    numberOfNodeWithExcess ++;
//////////////////////////////////////////////////////////////////////
//                    relabel( j );  第二次relabel 对j
//                  为啥对这个点也要进行relabel？，，  overhead？？？
                    int temp = relabel( j );
                    if(temp == UNFEASIBLE || temp == PRICE_OFL)
                        return temp;
//////////////////////////////////////////把j也加入到push队列中~
                    insertToExcessQueue( j );
                }
                totalExcess += j_exc;
            }
            else {
                totalExcess -= df;
            }
        }
//当前再溢出
        if ( i->acessExs() <= 0) numberOfNodeWithExcess --;
        if ( i->acessExs() <= 0 || priceFlag ) break;
//////////////////////////////////////////////////////////////////
//        relabel( i );
        //如果节点i依然 excess 则接着进行relabel操作选择合适的边进行push
        int temp = relabel( i );
        if(temp == UNFEASIBLE || temp == PRICE_OFL)
            return temp;
//////////////////////////////////////////////////////////////////
//        上面relabel了一下，重新更新了current 的值，取这个current的值
        a = i->current();
        j = a->head();
    }

    i->setCurrent( a);
    return 0;
}

int FindServer::priceIn()
{
    NODE *i; // current node
    NODE *j;
    ARC *a; // current arc from i
    ARC *a_stop; // first arc from the next node
    ARC *b; // arc to be exchanged with suspended
    ARC *ra; // opposite to a
    ARC *rb; // opposite to b
    priceType rc; // reduced cost
    int n_in_bad; // number of priced_in arcs with negative reduced cost
    int bad_found; // if 1 we are at the second scan if 0 we are at the first scan
    excessType i_exc; // excess of i
    excessType df; // an amount to increase flow


    bad_found = 0;
    n_in_bad = 0;

 restart:

    for ( i = nodeArray; i != sentinelNode; i ++ ) {

        for ( a = i->first() - 1, a_stop = i->suspended() - 1; a != a_stop; a -- ) {

            rc = REDUCED_COST( i, a->head(), a );
            if ( ( rc < 0) && ( a->acessResidualCapcity() > 0) ) { // bad case;
                if ( bad_found == 0 ) {
                    bad_found = 1;
                    updateCutOff();
                    goto restart;
                }
                df = a->acessResidualCapcity();
                increase_flow( i, a->head(), a, df );

                ra = a->sister();
                j  = a->head();

                i->decreaseFirst();
                b = i->first();
                exchange( a, b );

                if ( SUSPENDED( j, ra ) ) {
                    j->decreaseFirst();
                    rb = j->first();
                    exchange( ra, rb );
                }

                n_in_bad ++;
            }
            else {
                if ( ( rc < cutOn ) && ( rc > -cutOn ) ) {
                    i->decreaseFirst();
                    b = i->first();
                    exchange( a, b );
                }
            }
        }
    }


    if ( n_in_bad != 0 ) {

        numberBadPriceIn ++;

        // recalculating excess queue;
        totalExcess = 0;
        numberOfNodeWithExcess = 0;
        rstExcessQueue();

        for ( i = nodeArray; i != sentinelNode; i ++ ) {
            i->setCurrent( i->first());
            i_exc = i->acessExs();
            if ( i_exc > 0 ) { // i is a source;
                totalExcess += i_exc;
                numberOfNodeWithExcess ++;
                insertToExcessQueue( i );
            }
        }

        insertToExcessQueue( AddressOftmpNode );
    }

    if ( timeForPriceIn == TIME_FOR_PRICE_IN2)
        timeForPriceIn = TIME_FOR_PRICE_IN3;
    if ( timeForPriceIn == TIME_FOR_PRICE_IN1)
        timeForPriceIn = TIME_FOR_PRICE_IN2;

    return ( n_in_bad);
}

int FindServer::refine()
{
    NODE *i; // current node
    excessType i_exc; // excess of i
//    long np, nr, ns; // variables for additional print
    int pr_in_int; // current number of updates between price_in

//    np = numberPush;
//    nr = numberRelabel;
//    ns = numberScan;

    numberRefine ++;
//    和nRef 有什么区别。

    nRef ++;
//    第一次 nRef ==1 refine 次数为1
    nRel = 0;
    pr_in_int = 0;

    // initialize;
    totalExcess = 0;
    numberOfNodeWithExcess = 0;
    rstExcessQueue();

    timeForPriceIn = TIME_FOR_PRICE_IN1;

    for ( i = nodeArray; i != sentinelNode; i ++ ) {
        i->setCurrent( i->first());
        i_exc = i->acessExs();
        if ( i_exc > 0 ) { // i  is a source
            totalExcess += i_exc;
            numberOfNodeWithExcess++;
//            加入到要被push的excess 队列中
            insertToExcessQueue( i );
        }
    }

    if ( totalExcess <= 0 ) return 0;

    // (2) main loop

    while ( 1 ) {
//判断当前是否还有溢出的节点
        if ( isEmptyExcessQ() ) {
            if ( nRef > PRICE_OUT_START ) {
                pr_in_int = 0;
                priceIn();
            }

            if ( isEmptyExcessQ() ) break;
        }

        REMOVE_FROM_EXCESS_Q( i );

        // push all excess out of i
        if ( i->acessExs() > 0 ) {
//////////////////////////////////////////////////////////////////
//            discharge( i );
            int temp = discharge( i );
            if(temp == UNFEASIBLE || temp == PRICE_OFL)
                return temp;
//////////////////////////////////////////////////////////////////
            if ( timeForUpdate() || priceFlag ) {
                if ( i->acessExs() > 0 ) {
                    insertToExcessQueue( i );
                }

                if ( priceFlag && ( nRef > PRICE_OUT_START ) ) {
                    pr_in_int = 0;
                    priceIn();
                    priceFlag = 0;
                }

                priceUpdate();

                while ( updateFlag ) {
                    if ( nRef == 1 ) {
                        return UNFEASIBLE ;
                    } else {
                        updateFlag = 0;
                        updateCutOff();
                        numberBadRelabel ++;
                        pr_in_int = 0;
                        priceIn();
                        priceUpdate();
                    }
                }
                nRel = 0;

                if ( nRef > PRICE_OUT_START && (pr_in_int ++ > timeForPriceIn) ) {
                    pr_in_int = 0;
                    priceIn();
                }
            }
        }
    }

    return 0;
}

int FindServer::priceRefine()
{
    NODE *i; // current node
    NODE *j; // opposite node
    NODE *ir; // nodes for passing over the negative cycle
    NODE *is;
    ARC *a; // arc (i,j)
    ARC *a_stop; // first arc from the next node
    ARC *ar;
    long bmax;            // number of farest nonempty bucket
    long i_rank;          // rank of node i
    long j_rank;         // rank of node j
    long j_new_rank;      // new rank of node j
    BUCKET *b;              // current bucket
    BUCKET *b_old;          // old and new buckets of current node
    BUCKET *b_new;
    priceType rc = 0; // reduced cost of a
    priceType dr; // ranks difference
    priceType dp;
    int cc;
    // return code: 1 - flow is epsilon optimal
    // 0 - refine is needed
    long df; // cycle capacity
    int nnc; // number of negative cycles cancelled during one iteration
    int snc; // total number of negative cycle cancelled

    numberPRefine ++;

    cc = 1;
    snc = 0;

    maxCycleCancel = ( nRef >= START_CYCLE_CANCEL) ? MAX_CYCLES_CANCELLED : 0;


    // (1) main loop
    // while negative cycle is found or eps-optimal solution is constructed
    while ( 1 ) {

        nnc = 0;
        for ( i = nodeArray; i != sentinelNode; i ++ ) {
            i->setRank( 0);
            i->setInp( WHITE);
            i->setCurrent( i->first());
        }
        rstStackQ();

        for ( i = nodeArray; i != sentinelNode; i ++ ) {
            if ( i->inp() == BLACK ) continue;

            i->setBucketNext( NULL);

            // deapth first search
            while ( 1 ) {
                i->setInp( GREY);

                // scanning arcs from node i starting from current
                for ( a = i->current(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        if ( REDUCED_COST ( i, j, a ) < 0 ) {
                            if ( j->inp() == WHITE ) { // fresh node  - step forward
                                i->setCurrent( a);
                                j->setBucketNext( i);
                                i = j;
                                a = j->current();
                                a_stop = (j+1)->suspended();
                                break;
                            }

                            if ( j->inp() == GREY ) { // cycle detected
                                cc = 0;
                                nnc ++;
                                i->setCurrent( a);
                                is = ir = i;
                                df = MAX_32;

                                while ( 1 ) {
                                    ar = ir->current();
                                    if ( ar->acessResidualCapcity() <= df ) {
                                        df = ar->acessResidualCapcity();
                                        is = ir;
                                    }
                                    if ( ir == j ) break;
                                    ir = ir->acessBNext();
                                }

                                ir = i;

                                while ( 1 ) {
                                    ar = ir->current();
                                    increase_flow( ir, ar->head(), ar, df);
                                    if ( ir == j ) break;
                                    ir = ir->acessBNext();
                                }

                                if ( is != i ) {
                                    for ( ir = i; ir != is; ir = ir->acessBNext() ) {
                                        ir->setInp( WHITE);
                                    }
                                    i = is;
                                    a = is->current() + 1;
                                    a_stop = (is+1)->suspended();
                                    break;
                                }
                            }
                        }
                        // if j-color is BLACK - continue search from i
                    }
                } // all arcs from i are scanned

                if ( a == a_stop ) {
                    // step back
                    i->setInp( BLACK);
                    numberPreScan1 ++;
                    j = i->acessBNext();
                    stackQueuePush( i );
                    if ( j == NULL ) break;
                    i = j;
                    i->increaseCurrent();
                }

            } // end of deapth first search
        } // all nodes are scanned


        // () no negative cycle
        // computing longest paths with eps-precision

        snc += nnc;
        if ( snc < maxCycleCancel ) cc = 1;
        if ( cc == 0 ) break;
        bmax = 0;

        while ( isNotEmptyStackQ() ) {

            numberPreScan2 ++;
            STACKQ_POP( i );
            i_rank = i->rank();
            for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

                if ( OPEN( a ) ) {
                    j  = a->head();
                    rc = REDUCED_COST( i, j, a );

                    if ( rc < 0 ) { // admissible arc;
                        dr = (priceType) (( - rc - 0.5 ) / epsilon);
                        if (( j_rank = dr + i_rank ) < numberOfBucket ) {
                            if ( j_rank > j->rank() )
                                j->setRank( j_rank);
                        }
                    }
                }
            } // all arcs from i are scanned

            if ( i_rank > 0 ) {
                if ( i_rank > bmax ) bmax = i_rank;
                b = bucketArray + i_rank;
                insertToBucket( i, b );
            }
        } // end of while-cycle: all nodes are scanned - longest distancess are computed;


        if ( bmax == 0 ) // preflow is eps-optimal;
            { break; }


        for ( b = bucketArray + bmax; b != bucketArray; b -- ) {
            i_rank = b - bucketArray;
            dp = i_rank * epsilon;

            while ( isNotEmptyBucket( b) ) {
                GET_FROM_BUCKET( i, b );
                numberPreScan ++;

                for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        j_rank = j->rank();
                        if ( j_rank < i_rank ) {
                            rc = REDUCED_COST( i, j, a );
                            if ( rc < 0 ) {
                                j_new_rank = i_rank;
                            } else {
                                dr = rc / epsilon;
                                j_new_rank = ( dr < numberOfBucket ) ? i_rank - ( (long)dr + 1 ) : 0;
                            }
                            if ( j_rank < j_new_rank ) {
                                if ( cc == 1 ) {
                                    j->setRank( j_new_rank);
                                    if ( j_rank > 0 ) {
                                        b_old = bucketArray + j_rank;
                                        REMOVE_FROM_BUCKET( j, b_old );
                                    }
                                    b_new = bucketArray + j_new_rank;
                                    insertToBucket( j, b_new );
                                }
                                else {
                                    df = a->acessResidualCapcity();
                                    increase_flow( i, j, a, df );
                                }
                            }
                        }
                    } // end if opened arc
                } // all arcs are scanned

                i->decreasePrice( dp);

            } // end of while-cycle: the bucket is scanned
        } // end of for-cycle: all buckets are scanned

        if ( cc == 0 ) break;

    } // end of main loop



    // (2) finish
    // if refine needed - saturate non-epsilon-optimal arcs;

    if ( cc == 0 ) {
        for ( i = nodeArray; i != sentinelNode; i ++) {
            for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                if ( REDUCED_COST( i, a->head(), a ) < - epsilon ) {
                    if ( ( df = a->acessResidualCapcity() ) > 0 ) {
                        increase_flow( i, a->head(), a, df );
                    }
                }
            }
        }
    }

    return ( cc );
}

void FindServer::computePrices()
{
    NODE *i; // current node
    NODE *j; // opposite node
    ARC *a; // arc (i,j)
    ARC *a_stop; // first arc from the next node
    long bmax; // number of farest nonempty bucket
    long i_rank; // rank of node i
    long j_rank; // rank of node j
    long j_new_rank; // new rank of node j
    BUCKET *b; // current bucket
    BUCKET *b_old; // old and new buckets of current node
    BUCKET *b_new;
    priceType rc; // reduced cost of a
    priceType dr; // ranks difference
    priceType dp;
    int cc; // return code: 1 - flow is epsilon optimal 0 - refine is needed

    numberPRefine ++;
    cc = 1;

    // (1) main loop
    // while negative cycle is found or eps-optimal solution is constructed
    while ( 1 ) {

        for ( i = nodeArray; i != sentinelNode; i ++) {
            i->setRank( 0);
            i->setInp( WHITE);
            i->setCurrent( i->first());
        }
        rstStackQ();

        for ( i = nodeArray; i != sentinelNode; i ++ ) {
            if ( i->inp() == BLACK ) continue;

            i->setBucketNext( NULL);
            // depth first search
            while ( 1 ) {
                i->setInp( GREY);

                // scanning arcs from node i
                for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        if ( REDUCED_COST( i, j, a ) < 0 ) {
                            if ( j->inp() == WHITE ) { // fresh node  - step forward
                                i->setCurrent( a);
                                j->setBucketNext( i);
                                i = j;
                                a = j->current();
                                a_stop = (j+1)->suspended();
                                break;
                            }

                            if ( j->inp() == GREY ) { // cycle detected; should not happen
                                cc = 0;
                            }
                        }
                        // if j-color is BLACK - continue search from i
                    }
                } // all arcs from i are scanned

                if ( a == a_stop ) {
                    // step back
                    i->setInp( BLACK);
                    numberPreScan1 ++;
                    j = i->acessBNext();
                    stackQueuePush( i );
                    if ( j == NULL ) break;
                    i = j;
                    i->increaseCurrent();
                }

            } // end of deapth first search
        } // all nodes are scanned


        // no negative cycle
        // computing longest paths

        if ( cc == 0 ) break;
        bmax = 0;

        while ( isNotEmptyStackQ() ) {
            numberPreScan2 ++;
            STACKQ_POP( i );
            i_rank = i->rank();
            for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                if ( OPEN( a ) ) {
                    j  = a->head();
                    rc = REDUCED_COST( i, j, a );


                    if ( rc < 0 ) {// admissible arc
                        dr = - rc;
                        if (( j_rank = dr + i_rank ) < numberOfBucket ) {
                            if ( j_rank > j->rank() )
                                j->setRank( j_rank);
                        }
                    }
                }
            } // all arcs from i are scanned

            if ( i_rank > 0 ) {
                if ( i_rank > bmax ) bmax = i_rank;
                b = bucketArray + i_rank;
                insertToBucket( i, b );
            }
        } // end of while-cycle: all nodes are scanned - longest distancess are computed;

        if ( bmax == 0 )
            { break; }

        for ( b = bucketArray + bmax; b != bucketArray; b -- ) {
            i_rank = b - bucketArray;
            dp = i_rank;

            while ( isNotEmptyBucket( b) ) {
                getFromBucket( i, b );
                numberPreScan ++;

                for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        j_rank = j->rank();
                        if ( j_rank < i_rank ) {
                            rc = REDUCED_COST( i, j, a );

                            if ( rc < 0 ) {
                                j_new_rank = i_rank;
                            } else {
                                dr = rc;
                                j_new_rank = ( dr < numberOfBucket ) ? i_rank - ( (long)dr + 1 ) : 0;
                            }
                            if ( j_rank < j_new_rank ) {
                                if ( cc == 1 ) {
                                    j->setRank( j_new_rank);
                                    if ( j_rank > 0 ) {
                                        b_old = bucketArray + j_rank;
                                        rmFromBucket( j, b_old );
                                    }
                                    b_new = bucketArray + j_new_rank;
                                    insertToBucket( j, b_new );
                                }
                            }
                        }
                    } // end if opened arc
                } // all arcs are scanned

                i->decreasePrice( dp);

            } // end of while-cycle: the bucket is scanned
        } // end of for-cycle: all buckets are scanned

        if ( cc == 0 ) break;

    } // end of main loop
}

void FindServer::priceOut()
{
    NODE *i; // current node
    ARC *a; // current arc from i
    ARC *a_stop; // first arc from the next node
    ARC *b; // arc to be exchanged with suspended
    double numbelCutOff; // -cut_off
    double rc; // reduced cost

    numbelCutOff = - cutOff;

    for ( i = nodeArray; i != sentinelNode; i ++) {
        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

            rc = REDUCED_COST( i, a->head(), a );
            if ( ( rc > cutOff && CLOSED(a->sister()) ) ||
                 ( rc < numbelCutOff && CLOSED(a) ) ) { // suspend the arc

                b = i->first();
                i->increaseFirst();
                exchange( a, b );
            }
        }
    }
}

int FindServer::updateEpsilon()
{
    // decrease epsilon after epsilon-optimal flow is constructed;
    if ( epsilon <= 1 ) return ( 1 );

    epsilon = (priceType) (ceil ( (double) epsilon / scaleFactor ));
    cutOff = cutOffFactor * epsilon;
    cutOn = cutOff * CUT_OFF_GAP;

    return ( 0 );
}

int FindServer::checkFeasible()
{
    if ( checkSolution == false)
        return ( 0);

    NODE *i;
    ARC *a, *a_stop;
    long fa;
    int ans = 1;

    for ( i = nodeArray; i != sentinelNode; i ++) {
        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
            if ( capArray[ ARC_NUMBER(a) ] > 0) {
                fa = capArray[ ARC_NUMBER(a) ] - a->acessResidualCapcity();
                if ( fa < 0) {
                    ans = 0;
                    break;
                }
                nodeBalance[ i - nodeArray ] -= fa;
                nodeBalance[ a->head() - nodeArray ] += fa;
            }
        }
    }

    for ( i = nodeArray; i != sentinelNode; i ++) {
        if ( nodeBalance[ i - nodeArray ] != 0) {
            ans = 0;
            break;
        }
    }

    return ( ans);
}

int FindServer::checkCs()
{
    // check complimentary slackness;
    NODE *i;
    ARC *a, *a_stop;

    for ( i = nodeArray; i != sentinelNode; i ++) {
        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

            if ( OPEN(a) && (REDUCED_COST(i, a->head(), a) < 0) ) {
                return ( 0);
            }
        }
    }
    return(1);
}

int FindServer::checkEpsOpt()
{
    NODE *i;
    ARC *a, *a_stop;

    for ( i = nodeArray; i != sentinelNode; i ++) {
        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

            if ( OPEN(a) && (REDUCED_COST(i, a->head(), a) < - epsilon) ) {
                return ( 0);
            }
        }
    }
    return(1);
}

void FindServer::initSolution()
{
    ARC *a; // current arc (i,j)
    NODE *i; // tail of a
    NODE *j; // head of a
    long df; // residual capacity

    for ( a = arcArray; a != sentinelArc; a ++ ) {
        if ( a->acessResidualCapcity() > 0 && a->cost() < 0 ) {
            df = a->acessResidualCapcity();
            i  = a->sister()->head();
            j  = a->head();
            increase_flow( i, j, a, df );
        }
    }
}

void FindServer::csCostReinit()
{
    if ( costRestart == false)
        return;

    NODE *i; // current node
    ARC *a;          // current arc
    ARC *a_stop;
    BUCKET *b; // current bucket
    priceType rc, minc, sum;


    for ( b = bucketArray; b != lastBucket; b ++) {
        rstBucket( b);
    }

    rc = 0;
    for ( i = nodeArray; i != sentinelNode; i ++) {
        rc = MIN(rc, i->price());
        i->setFirst( i->suspended());
        i->setCurrent( i->first());
        i->setQueueNext( sentinelNode);
    }

    // make prices nonnegative and multiply
    for ( i = nodeArray; i != sentinelNode; i ++) {
        i->setPrice( (i->price() - rc) * _dn);
    }

    // multiply arc costs
    for (a = arcArray; a != sentinelArc; a ++) {
        a->multiplyCost( _dn);
    }

    sum = 0;
    for ( i = nodeArray; i != sentinelNode; i ++) {
        minc = 0;
        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
            if ( (OPEN(a) && ((rc = REDUCED_COST(i, a->head(), a)) < 0)) )
                minc = MAX( epsilon, -rc);
        }
        sum += minc;
    }

    epsilon = ceil(sum / _dn);

    cutOffFactor = CUT_OFF_COEF * pow((double)_n, CUT_OFF_POWER);

    cutOffFactor = MAX( cutOffFactor, CUT_OFF_MIN);

    nRef = 0;

    numberRefine = numberDiscarge = numberPush = numberRelabel = 0;
    numberUpdate = numberScan = numberPRefine = numberPreScan = numberPreScan1 =
        numberBadPriceIn = numberBadRelabel = 0;

    priceFlag = 0;

    excessQueueFirst = NULL;
}

int FindServer::cs2CostRestart( double *objective_cost)
{
    // restart after a cost update;
    if ( costRestart == false)
        return 0;

    int cc; // for storing return code;

    printf("c \nc ******************************\n");
    printf("c Restarting after a cost update\n");
    printf("c ******************************\nc\n");

    csCostReinit();

    printf ("c Init. epsilon = %6.0f\n", epsilon);
    cc = updateEpsilon();

    if (cc != 0) {
        printf("c Old solution is optimal\n");
    }
    else {
        do { // scaling loop
            while ( 1 ) {
                if ( ! priceRefine() )
                    break;

                if ( nRef >= PRICE_OUT_START ) {
                    if ( priceIn() )
                        break;
                }
                if ((cc = updateEpsilon ()))
                    break;
            }
            if (cc) break;
//////////////////////////////////////////////
//            refine();
            int temp = refine();
            if(temp == UNFEASIBLE || temp == PRICE_OFL)
                return temp;
//////////////////////////////////////////
            if ( nRef >= PRICE_OUT_START ) {
                priceOut();
            }
            if ( updateEpsilon() )
                break;
        } while ( cc == 0 );
    }

    finishup( objective_cost );
     return 0;
}

void FindServer::printSolution()
{
    if ( printAnswer == false)
        return;

    NODE *i;
    ARC *a;
    long ni;
    priceType cost;
//    printf ("c\ns %.0l\n", cost );

    for ( i = nodeArray; i < nodeArray + _n; i ++ ) {
        ni = NODE_NUMBER( i );
        for ( a = i->suspended(); a != (i+1)->suspended(); a ++) {
            if ( capArray[ ARC_NUMBER (a) ]  > 0 ) {
//                printf("f %7ld %7ld %10ld\n",
//                       ni, NODE_NUMBER(a->head()), capArray[ ARC_NUMBER(a) ] - a->acessResidualCapcity());
            }
        }
    }

    // COMP_DUALS?
    if ( computeDuals == true) { // find minimum price;
        cost = MAX_32;
        for ( i = nodeArray; i != sentinelNode; i ++) {
            cost = MIN(cost, i->price());
        }
        for ( i = nodeArray; i != sentinelNode; i ++) {
//            printf("p %7ld %7.2lld\n", NODE_NUMBER(i), i->price() - cost);
        }
    }

//    printf("c\n");
}

void FindServer::printGraph()
{
    NODE *i;
    ARC *a;
    long ni, na;
//    printf ("\nGraph: %d\n", _n);
    for ( i = nodeArray; i < nodeArray + _n; i ++ ) {
        ni = NODE_NUMBER( i );
//        printf("\nNode %d", ni);
        for ( a = i->suspended(); a != (i+1)->suspended(); a ++) {
            na = ARC_NUMBER( a );
            printf("\n {%d} %d -> %d  cap: %d  cost: %d", na,
                ni, NODE_NUMBER(a->head()), capArray[ARC_NUMBER(a)], a->cost());
        }
    }
}

void FindServer::finishup( double *objective_cost)
{
    ARC *a; // current arc
    long na; // corresponding position in capacity array
    double obj_internal = 0; // objective
    priceType cs; // actual arc cost
    long flow; // flow through an arc
    NODE *i;

    // (1) NO_ZERO_CYCLES?
    if ( numberZeroCycles == true) {
        for ( a = arcArray; a != sentinelArc; a ++ ) {
            if ( a->cost() == 1) {
                assert( a->sister()->cost() == -1);
                a->setCost( 0);
                a->sister()->setCost( 0);
            }
        }
    }

    // (2)
    for ( a = arcArray, na = 0; a != sentinelArc ; a ++, na ++ ) {
        cs = a->cost() / _dn;
        if ( capArray[na]  > 0 && (flow = capArray[na] - a->acessResidualCapcity()) != 0 )
            obj_internal += (double) cs * (double) flow;
        a->setCost( cs);
    }

    for ( i = nodeArray; i != sentinelNode; i ++) {
        i->setPrice( (i->price() / _dn));
    }

    // (3) COMP_DUALS?
    if ( computeDuals == true) {
        computePrices();
    }

    *objective_cost = obj_internal;
}

int FindServer::cs2( double *objective_cost)
{
    // the main calling function;
    int cc = 0; // for storing return code;


    // (1) update epsilon first;
    updateEpsilon();


    // (2) scaling loop;
    do {
//////////////////////////////////////////////
//            refine();
             int temp = refine();
             if(temp == UNFEASIBLE || temp == PRICE_OFL)
                return temp;
//////////////////////////////////////////

        if ( nRef >= PRICE_OUT_START )
            priceOut();

        if ( updateEpsilon() )
            break;


    } while ( cc == 0 );


    // (3) finishup;
    finishup( objective_cost );
    return 0;
}

int FindServer::runSolition()
{


    double objective_cost;


    // (4) ordering, etc.;
    preProscessing();


    // () CHECK_SOLUTION?
    if ( checkSolution == true) {
        nodeBalance = (long long int *) calloc (_n+1, sizeof(long long int));
        for ( NODE *i = nodeArray; i < nodeArray + _n; i ++ ) {
            nodeBalance[i - nodeArray] = i->acessExs();
        }
    }


    // (5) initializations;
    _m = 2 * _m;
    if(initialization() == ALLOCATION_FAULT)
        return ALLOCATION_FAULT; // works already with 2*m;

    int temp = cs2( &objective_cost );
    if(temp == UNFEASIBLE || temp == PRICE_OFL)
       return temp;
    ////////////////////////////////////////////////
    total_cost = objective_cost;

    if ( printAnswer == true ) {
        printSolution();
    }

    // () cleanup;
    deleteMemory();
    return 0;
}

void FindServer::increase_flow(FindServer::NODE *i, FindServer::NODE *j, FindServer::ARC *a, long df)
{
    //从出点获取averageRouteCost
    i->averageRouteCost = (double)i->totalCost/i->acessExs();
    averageExcessCostVecmoniRaw[NODE_NUMBER(i)].averageRouteCost = i->averageRouteCost;
    i->decreaseExcess( df);
    j->increaseExcess( df);
    a->decreaseResidualCapacity( df);
    a->sister()->increaseResidualCapacity( df);
    // 这条路径上产生（单位路径费用+前面的单位路径费用）*流大小
    double flowCost = (i->averageRouteCost+a->cost())*df;
//    入点更新：（费用都是由入点来付） 入点j，出点i。
    j->totalCost += flowCost;  // 路径总费用更新
    //单位流总成本更新
    //选择合适的服务器
    j->serverId = fitServer(j->acessExs());
    int serverCost = INT_MAX;
//    可以找到可用的服务器
    if(j->serverId > 0)
    {
        serverCost = serverCostVect[j->serverId].cost;
//        更新averageExcessCost
        double averageCostTemp = (double)(j->totalCost + serverCost*(_n+1) + nodeBaseCost[NODE_NUMBER(j)]*(_n+1)) / j->acessExs();
        if(averageCostTemp < j->averageExcessCost)
         {
            // 更新averageExcessCost
         j->averageExcessCost = averageCostTemp;
            //              主要更新3个量 ,因为nodeID不改变 所以就不push了
         averageExcessCostVecmoniRaw[NODE_NUMBER(j)].averageExcessCost = j->averageExcessCost;
         averageExcessCostVecmoniRaw[NODE_NUMBER(j)].nodeBaseCost = nodeBaseCost[NODE_NUMBER(j)];
         averageExcessCostVecmoniRaw[NODE_NUMBER(j)].serverID = j->serverId;
         averageExcessCostVecmoniRaw[NODE_NUMBER(j)].excessMax = j->acessExs();
         j->totalCostMax = j->totalCost;
         j->excessMax = j->acessExs();
        }
        if(j->acessExs() > 0)// 只有超级汇点的excess 会小于0
        j->averageRouteCost = j->totalCost / j->acessExs();// 入点需要更新averageRouteCost 出点不需要
        averageExcessCostVecmoniRaw[NODE_NUMBER(j)].averageRouteCost = j->averageRouteCost;
    }
    else
    {
//       无法找到可用的服务器
//        只更新了路径总费用totalCost，不更新averageExcessCost 直接下一步
    }
//    else
//       j->averageRouteCost = -1;


// 出点更新（费用是由入点来付，出点只需要把excess 删掉流行了）
   i->totalCost = i->acessExs() *  i->averageRouteCost;
   // 如果i的excess 也大于零
   if(i->acessExs() > 0)
   {
       i->serverId = fitServer(i->acessExs());
       if(i->serverId > 0)
       {
         serverCost = serverCostVect[j->serverId].cost;
         double averageCostTemp = (i->totalCost + serverCost*(_n+1) + nodeBaseCost[NODE_NUMBER(i)]*(_n+1)) / i->acessExs();
         if(averageCostTemp < i->averageExcessCost)
         {
           i->averageExcessCost = averageCostTemp;
//           只更新这三个量
           averageExcessCostVecmoniRaw[NODE_NUMBER(i)].averageExcessCost = i->averageExcessCost;
           averageExcessCostVecmoniRaw[NODE_NUMBER(i)].nodeBaseCost = nodeBaseCost[NODE_NUMBER(i)];
           averageExcessCostVecmoniRaw[NODE_NUMBER(i)].serverID = i->serverId;
            averageExcessCostVecmoniRaw[NODE_NUMBER(j)].excessMax = i->acessExs();
           i->totalCostMax = i->totalCost;
           i->excessMax = i->acessExs();
         }
       }
      else
       {
//          没有找到也不更新
       }

   }
}

int start;
int nodeEnd;
int capacity;
int cost;

int consumerIndex;
int consumerNode;
int consumerDemand;

int totalDemand = 0;
int countEdge = 0;
int countConsumer = 0;



vector<ConsumerDatatype> consumerVector;


vector<int> nodeServerDelete;
int countServerDeleted = 0;
double costPre;
double costNow;
long long costFinal =  MAX_32;
vector<int>::iterator it;
MinCostSolver *mcfProblem;
int viliate(graph &G, vector<int> serverVec)
{
//   服务器数量（排除不具有服务器资格的节点）
//    先清空使用边数组
    usedARC.clear();
    mcfProblem = new MinCostSolver( nodeNum+1, 2*edgeNum+serverVec.size());
    int j = 1;
    for (vector<DirectedEdge> v:G.adj)
    {

        for(DirectedEdge e:v)
        {
            mcfProblem->setArc( e.startNode, e.endNode, 0, e.capacity, e.cost);
        }
    }
//    添加消费节点
    for(ConsumerDatatype cons:consumerVector)
    {
        mcfProblem->setDemandOfNode(cons.nodeID,cons.demand);
    }
    // 添加服务器节点
    for(int i : serverVec)
    {
        if(averageExcessCostVecmoniRaw[i].serverID < 0)
            return -2;
        int cap = serverCostVect[averageExcessCostVecmoniRaw[i].serverID].capacity;

        mcfProblem->setArc(i,nodeNum,0,cap,0);
    }
   // 添加超级汇点
    mcfProblem->setDemandOfNode( nodeNum, -totalDemand);
    if(mcfProblem->runSolition() == -2)
    {
         delete mcfProblem;
        return -2;
    }
    else
    {
         delete mcfProblem;
        return total_cost + calculateCost(serverVec);
    }

}
//You need to complete the function
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
//init graph
    DirectedEdge DirectedEdgeTemp;

     char numTemp[5] = {' ',' ',' ',' '};// 存1w 还是1K？
//     当前行读了多少
     int count = 0;
//     numTemp写的位置
     int j = 0;
	 int line = 0;// 当前读的文件行数
    // Output demo
//1.      取节点数和边数
     for( int i = 0; topo[0][i] != '\n'; i++)
     {
         if( topo[0][i] == ' ')
         {
             switch(count)
             {
                case 0: nodeNum = atoi(numTemp); break;
                case 1: edgeNum = atoi(numTemp); break;
             }
             count ++;
             j = 0;
             numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' ';
         }
         else
         {
             numTemp[j] = topo[0][i];
             j++;
         }
     }

     consumerNum = atoi(numTemp);
     //         重置临时变量
     count = 0;j=0;
     numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' ';
     graph graphPtrData(nodeNum);
     FindServer findServerProblem( nodeNum+1, 2*edgeNum+nodeNum);
////////////////////////////////////////
//2.      读服务器费用和输出能力 第二行开始
	line = 2;
	while(topo[line][0] != '\n')
	{
	 // 	读一行
	 int capacityTemp;
	 int costTemp;
	 for( int i = 0; topo[line][i] != '\n'; i++)
     {
         if( topo[line][i] == ' ')
         {
             switch(count)
             {
                case 0:  break;// 服务器序号这个量不管了
                case 1: capacityTemp = atoi(numTemp); break; // 读服务器cap
             }
             count ++;
             j = 0;
             numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' ';
         }
         else
         {
             numTemp[j] = topo[line][i];
             j++;
         }
     }
	 costTemp = atoi(numTemp); // 读服务器费用。
//     serverCostVect.push_back({capacityTemp，costTemp});
     serverCostVect.push_back({capacityTemp,costTemp});
     //    重置临时变量
	 count = 0;j=0;
     numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' ';
	 line ++;
	}
    line ++;//碰到首元素是空格了，切换到下一行
////////////////////////////////////////
//3.      读节点架设成本
	while(topo[line][0] != '\n')
	{
	 // 	读一行
	 for( int i = 0; topo[line][i] != '\n'; i++)
     {
             numTemp[j] = topo[line][i];
             j++;
             if(topo[line][i] == ' ')
			 {
				 j = 0;// 前一个就不要了 从头开始读
				 numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' '; 
			 }		 
     }
	 nodeBaseCost.push_back(atoi(numTemp)); // 把架设成本放到nodeBaseCost 中
     //    重置临时变量
	 j = 0; 
     numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' ';
	 line ++;
	}
	// 读服务器
	line ++; //碰到首元素是空格了，切换到下一行    
////////////////////////////////////////
//4.       读边
	int edgelineStart = line;
     for(; line < edgeNum +edgelineStart; line++)
     {
         for(int i = 0; topo[line][i] != '\n';i++)
         {
             if( topo[line][i] == ' ')
             {
                 switch(count)
                 {
                    case 0: DirectedEdgeTemp.startNode = atoi(numTemp); break;
                    case 1: DirectedEdgeTemp.endNode = atoi(numTemp); break;
                    case 2: DirectedEdgeTemp.capacity = atoi(numTemp); break;
//                    case 3: cost = atoi(numTemp); break;
                 }
                 count ++;
                 j = 0;
                 numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' ';
             }
             else
             {
                 numTemp[j] = topo[line][i];
                 j++;
             }
         }
         DirectedEdgeTemp.cost = atoi(numTemp);
//         重置临时变量
         count = 0;j=0;
         numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' ';
//         添加边

         findServerProblem.setArc( DirectedEdgeTemp.startNode, DirectedEdgeTemp.endNode, 0, DirectedEdgeTemp.capacity, DirectedEdgeTemp.cost);
         findServerProblem.setArc( DirectedEdgeTemp.endNode, DirectedEdgeTemp.startNode, 0, DirectedEdgeTemp.capacity, DirectedEdgeTemp.cost);
         countEdge++;
          graphPtrData.addEdge(DirectedEdgeTemp);
         {
              int temp = DirectedEdgeTemp.startNode;
          DirectedEdgeTemp.startNode = DirectedEdgeTemp.endNode;
          DirectedEdgeTemp.endNode = temp;
          }
          graphPtrData.addEdge(DirectedEdgeTemp);
    }
////////////////////////////////////////
//5.       读消费节点
 for(line = edgeNum + edgelineStart + 1; line < line_num; line++)
     {
         for(int i = 0; topo[line][i] != '\n';i++)
         {
             if( topo[line][i] == ' ')
             {
                 switch(count)
                 {
                    case 0: consumerIndex = atoi(numTemp); break;
                    case 1: consumerNode = atoi(numTemp); break;
//                    case 2: capacity = atoi(numTemp); break;
//                    case 3: cost = atoi(numTemp); break;
                 }
                 count ++;
                 j = 0;
                 numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' ';
             }
             else
             {
                 numTemp[j] = topo[line][i];
                 j++;
             }
         }
         consumerDemand = atoi(numTemp);
//         重置临时变量
         count = 0;j=0;
         numTemp[0] = ' ';numTemp[1] = ' ';numTemp[2] = ' ';numTemp[3] = ' ';numTemp[4] = ' ';
//         添加消费节点的需求
         findServerProblem.setDemandOfNode(consumerNode, consumerDemand);
         totalDemand += consumerDemand;
         consumerVector.push_back({consumerNode,consumerDemand});
//          consumerVector.push_back(consumerNode);
         countConsumer++;
     }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//      汇点   0号不要
     for(int i = 0;i<=nodeNum-1;i++)
     {
        findServerProblem.setArc( i, nodeNum, 0, totalDemand+1, MAX_64-1);
     }

    for(serverDataType d: serverCostVect)
    {
        if(d.cost > serverCostMax)
            serverCostMax = d.cost;
        serverCostMid += d.cost;
    }
    serverCostMid = serverCostMid/serverCostVect.size();
    findServerProblem.setDemandOfNode(nodeNum, -totalDemand);
    findServerProblem.runSolition();
    averageExcessCostVecmoniRaw.erase(averageExcessCostVecmoniRaw.end()-1);
    averageExcessCostVecmoni = averageExcessCostVecmoniRaw;
    sort(averageExcessCostVecmoni.begin(), averageExcessCostVecmoni.end());
    refineTheServer(graphPtrData,findServerProblem);
    long  costResult ;
    stack<int> unCheckedRank1Node;
    vector<int> RankCheckedVector;
    vector<int> serverNode;
    vector<int> juncNode;
    int iterNum = 0;
    int countCostUp = 0;
    for(averageExcessCostType i:refineVect)
    {
        serverNode.push_back(i.nodeID);
        costResult = viliate(graphPtrData,serverNode);
        if( costResult == UNFEASIBLE) // 不可达的情况
        {
            if(i.adjacentRank >= 1)
                //checkRank1 = 1; //有一个 rank为1的点push进去了 ，注意
                unCheckedRank1Node.push(i.nodeID);
                continue;
       }
        else if(costResult < costFinal) // cost 接着降低
        {
            //问题开始有解了，开始尝试去删除前面被push进去的 rank不为0的节点
            costFinal = costResult;
//            if(nodeNum < 800)
            if(i.adjacentRank >= 1)
                //checkRank1 = 1; //有一个 rank为1的点push进去了 ，注意
            unCheckedRank1Node.push(i.nodeID);
            continue;
       }
        else                            // cost 开始升高
        {
  //            先把点扔了
            serverNode.pop_back(); // 增加了费用 pop掉
            if(countCostUp > 0.4 * consumerNum)
                break;
            countCostUp ++;
//           再检查前面的点，之后会变的无序
           while(unCheckedRank1Node.size())
           {

               int N = unCheckedRank1Node.top();
               unCheckedRank1Node.pop();
               costResult = erasecheck(N, serverNode, graphPtrData);

               if(costResult > costFinal || costResult == -2)// 变大了或者变的不可达了
                {
                   RankCheckedVector.push_back(N);  // 这个点是有用的点
                   serverNode.push_back(N);         // 重新压回去
                }
               else
               {
                   costFinal = costResult;
                   juncNode.push_back(N);// 这个已经被扔掉的点感觉是不是被扔的太快了 越在前面的 约容易被
//                   这个点没用，扔了
                   countCostUp = 0;
               }
           }

           int minNode = -2;
//          检查完rankstack的点之后，进行置换操作////////////////////////////////////////////////////////////////////////////////////////////////////////
//           serverNode.push_back(i.nodeID);
//           costResult = viliate(graphPtrData,serverNode);
//           if(costResult < costFinal && costResult !=-2)
//           {
//               costFinal = costResult;
//               countCostUp = 0;
//               continue;
//           }
//           serverNode.pop_back();
/////////////////////////////////////////////////////////////////////////////////////////////////////
           minNode = -2;
//            加进来的点如果是消费者节点，则要和前面已经在的节点进行互换
           if(findInVect(consumerVector,i.nodeID))
           {
//                和前面所有的点进行删除插入比较
               vector<int> serverTemp = serverNode;
               for(int j=serverTemp.size()-1; j>0 ;j--)
               {
                   int temp = serverTemp[j]; //暂存
                   serverTemp[j] = i.nodeID; //换成是这个消费节点
                   costResult = viliate(graphPtrData,serverTemp);
                   if (costResult < costFinal && costResult != UNFEASIBLE)// 说明这个节点更好 那就置换掉
                   {
                       minNode = j;
                       costFinal = costResult; // 更新costFinal
                   }
                   serverTemp[j] = temp;// 变回去
               }
//               迭代完了，如果minNode ！=-2 说明曾找到过更小的，则把对应的位换掉。
               if(minNode != -2)
               {
//                  serverNode[minNode] = i.nodeID;//  换成这个节点。
                    countCostUp = 0;
                  //////////////////////
//                  int needCheckNode = serverNode[minNode]; // 保存这个要被可能被删掉的点
                  serverNode[minNode] = i.nodeID;//  换成这个节点。
                   minNode = -2; // 重置
                     vector<int>::iterator iter;
//                  找最小的删除， 从末尾开始找
                     int NodeToDelete;
                 do
                 {
                  NodeToDelete = -1;
                  iter = serverNode.end()-1;
                  for(int j=serverNode.size()-1; j>0 ;j--)
                  {
                    int temp = (*iter);
                   serverNode.erase(iter);
                   iter --;
                   costResult = viliate(graphPtrData,serverNode);
                   if(costResult < costFinal && costResult!=-2)// 可以删除 找最小
                    {
                       costFinal = costResult;
                       NodeToDelete = temp;  //  把对应的Node 号记下来。
                    }
                   serverNode.push_back(temp);

                 }
                 if(NodeToDelete != -1) // 找到能删除最小的点了
                  {
                      countCostUp = 0;
                      for(iter=serverNode.begin();iter < serverNode.end();iter++)
                      {
                          if((*iter) == NodeToDelete) // 找到并删掉
                              serverNode.erase(iter);
                      }
                  }
                }while(NodeToDelete != -1); // 一直找最小的删除。
               } // end 删除点
           }// end 加进来的是消费点
         }// end cost升高

        iterNum ++;
        cout << costFinal <<" iter:"<<iterNum<<endl;
    }// 循环结束
    //////////////////////////////////////////
    //////////////////////////////////////////
    // 前面的垃圾回收工作:  把前面 因为 rank=1 扔掉的点再回收一次
    // 1 加点
    if(nodeNum<700 && nodeNum > 200)
    for(int i : juncNode)
    {
        serverNode.push_back(i);
        costResult = viliate(graphPtrData,serverNode);
        if (costResult < costFinal && costResult != UNFEASIBLE)// 说明这个节点更好 那就置换掉
        {

            costFinal = costResult; // 更新costFinal
         }
        else
          serverNode.pop_back();
    }


    cout << "junc recycle" <<endl;
    int minNode = -2;
    int countIter = 0;
    if(nodeNum > 200)
    for (int i = juncNode.size()-1; i > MAX(juncNode.size()- 0.3 * consumerNum,0) ;i--)
    {
//        一次循环找最小
        int nodeNumber = juncNode[i];
     for(int j=serverNode.size()-1; j>0 ;j--)
     {
         if(countIter > 10000)
             break;
         countIter ++;
       int temp = serverNode[j]; //暂存
       serverNode[j] = nodeNumber; //换成是这个节点
       costResult = viliate(graphPtrData,serverNode);
       if (costResult < costFinal && costResult != UNFEASIBLE)// 说明这个节点更好 那就置换掉
         {
           minNode = j;
           costFinal = costResult; // 更新costFinal
         }
         serverNode[j] = temp;// 变回去
     }
// 迭代完了，如果minNode ！=-2 说明曾找到过可以替换的点，则把对应的位换掉。
    if(minNode != -2)
    {

        serverNode[minNode] = nodeNumber;
         minNode = -2;
    }
   }
    viliate(graphPtrData,serverNode);
    vector<int> pathArcPos;
    string fileWrite;
    int pathNum = 0;
    char a[60];
    int nodeId = 0;
   //FILE * fopen(const char * path,const char * mode);
    for(ConsumerDatatype cons: consumerVector)
    {
//        一个一个找消费节点的路径
        vector<outputEdgeType>::iterator it;
        int minFlow;  // 初始化容量最小的边 为初始边

      while(1)
      {
          pathArcPos.clear();
          if(cons.demand == 0)
              break;
          it = find(usedARC.begin(),usedARC.end(),cons.nodeID);//还能找到 说明还有消费节点还有具有residual的边 前面有重载==、
          if(it == usedARC.end())// 找不着了 说明全走完了
              break;
          minFlow = MIN((*it).flow,cons.demand);  // 初始化容量最小的边 为初始边的residual
          pathArcPos.push_back(it - usedARC.begin());// 把这条边压入路径容器
//        找容量最小的那条边的值，做为该路径的值。
         while((*it).endNode!=nodeNum)// 直到到服务器节点 一直找
         {
            // 找下一条边，通往服务器
            it = find(usedARC.begin(),usedARC.end(),(*it).endNode);
            // 把这条边的索引压入路径容器
            pathArcPos.push_back(it - usedARC.begin());
            // 如果新找到的边容量更小，则更新容量
            if((*it).flow < minFlow)
                minFlow = (*it).flow;

//            路径的终点
         }

//        找完单条路径之后。更新usedArc； 并输出路径
         for(int i = pathArcPos.size()-2; i>=0;i--)
            {
                usedARC[pathArcPos[i]].flow -= minFlow;
//                cout<<usedARC[pathArcPos[i]].endNode<<"--"<<'('<<usedARC[pathArcPos[i]].cost<<')'
//                   <<"->";
//                if(i != )
                sprintf(a,"%d ", usedARC[pathArcPos[i]].endNode);
                fileWrite += a;
//                sprintf()
            }
                 cons.demand -= minFlow;
//                cout<<usedARC[pathArcPos[0]].startNode<<" |||  flow:"<<minFlow<<endl;
                int serverID = averageExcessCostVecmoniRaw[usedARC[pathArcPos[0]].startNode].serverID;
                sprintf(a,"%d %d %d %d\n", usedARC[pathArcPos[0]].startNode,nodeId,minFlow,serverID);
                fileWrite += a;
                pathNum++;
        }
        nodeId ++;
    }

//    char * topo_file = (char *)"17\n\n0 8 0 20\n21 8 0 20\n9 11 1 13\n21 22 2 20\n23 22 2 8\n1 3 3 11\n24 3 3 17\n27 3 3 26\n24 3 3 10\n18 17 4 11\n1 19 5 26\n1 16 6 15\n15 13 7 13\n4 5 8 18\n2 25 9 15\n0 7 10 10\n23 24 11 23";
    char numPath[20];
    sprintf(numPath,"%d\n\n",pathNum);
    string s;
    s += numPath;
    s += fileWrite;
    //write_result(numPath, filename);
    write_result(s.c_str(), filename);
    cout << costFinal <<endl;
}

void refineTheServer(graph &G,FindServer &findServerPro)
{
   int refineClusterCount = 0;

   for(int i = 0;i < averageExcessCostVecmoni.size(); i++)
//       遍历这个待选点的所有边
   {
       if(refineClusterCount >= consumerNum*2)
           break;
       int rankSet = findMinExcessCost(refineVect,G.adj[averageExcessCostVecmoni[i].nodeID],findServerPro);
//       小于1 就计数
//       if(rankSet <= 2)
//       {
          averageExcessCostVecmoni[i].adjacentRank = rankSet;
          refineVect.push_back(averageExcessCostVecmoni[i]);
          refineClusterCount ++;
//       }
   }
}

int findMinExcessCost(vector<averageExcessCostType> &vin, vector<DirectedEdge> &edgeVCec, FindServer &findServerPro)
{
    vector<averageExcessCostType>::iterator it;
//    vector<averageExcessCostType> findResult;
    int rank0Count = 0;
    int rankMax = 0;
    int findFlag = 0;
    int rankTemp;
//    遍历所有的边
    for (DirectedEdge e:edgeVCec)
    {
//        查找这条出度的边有没有
        it = find(vin.begin(),vin.end(),e.endNode);
        if(it != vin.end())
        {
//            找着了至少一个
            findFlag = 1;

//          if(e.)
//        下面是判断 这个节点是不是从上一个节点流过来的
            // 首先判断这个节点的excess 是不是 比较小
          if(averageExcessCostVecmoniRaw[e.startNode].excessMax <= 1.5 * e.capacity)
         {
          int maxFlow = MIN(averageExcessCostVecmoniRaw[e.endNode].excessMax,e.capacity);
          long long routeCost = maxFlow * (e.cost  * (nodeNum+3) + averageExcessCostVecmoniRaw[e.endNode].averageRouteCost);

          double averageExcessCost = (routeCost + (nodeBaseCost[e.startNode]+serverCostVect[averageExcessCostVecmoniRaw[e.startNode].serverID].cost) * (nodeNum+3))/maxFlow;
////         说明该点是由上一点汇入的- -，rank 应该急剧升高
         if(averageExcessCost - 1 <= averageExcessCostVecmoniRaw[e.startNode].averageExcessCost)
            { rankTemp = (*it).adjacentRank+10;return rankTemp;}
          }
//          如果不是 那样，则进行正常处理 直接返回1级相邻
          else if((*it).adjacentRank == 0)
          {
                rank0Count ++;
                return rank0Count;
          }
          else if((*it).adjacentRank > rankMax)
          {
            rankMax = (*it).adjacentRank;
          }
        }
    }
//    没有找到
    if(findFlag == 0)
        return 0;

    if(rank0Count > 1)
        return 0;
    else
        return rankMax+1;
}

int findFuncForVector(vector<int> &vin,int val)
{
    vector<int>::iterator it;
    it = find(vin.begin(),vin.end(),val);
    if(it == vin.end())
        return 0;
    else
        return 1;
}
int erasecheck(int N,vector<int> &serverVec, graph &G)
{
//    while(st.size())
//    {
//        int N = st.top();
//        st.pop();
        for(vector<int>::iterator iter = serverVec.begin(); iter!=serverVec.end();)
        {
            if((*iter) == N)
                iter = serverVec.erase(iter);// 删除对应的节点
            else
                iter ++;
        }
       return viliate(G,serverVec);
//    }
}



















//如果没有找到返回-1
int fitServer(int excess)
{
//    遍历所有
    int minIndex = 0x7FFF;
    int minCost = 0x7FFF;
    for(int i = 0;i < serverCostVect.size();i++)
    {

      if((serverCostVect[i].capacity >= excess) && (serverCostVect[i].cost < minCost))
        {
          minIndex = i;
          minCost = serverCostVect[i].cost;
        }
    }
    if( minIndex == 0x7FFF)
        return -1;
    else
        return minIndex;
}

long calculateCost(vector<int> &serverNode)
{
    long cost = 0;
    for(int i: serverNode)
    {
        cost += averageExcessCostVecmoniRaw[i].nodeBaseCost + serverCostVect[averageExcessCostVecmoniRaw[i].serverID].cost;
    }
    return cost;
}
