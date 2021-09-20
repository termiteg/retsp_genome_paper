#=== configuration ===
CHAIN_RESULTS_DIR=lastz_results

TARGET_DIR=Rspe
QUERY_DIR=Mnat

AXTCHAIN_LINEAR_GAP=medium
AXTCHAIN_MIN_SCORE=5000

CHAINNET_MIN_SPACE=100
CHAINNET_MIN_SCORE=10000

#====


TARGET=$TARGET_DIR
QUERY=$QUERY_DIR

## step 1
## netting

## step 2

chainMergeSort $CHAIN_RESULTS_DIR/*chain > all.chain
chainPreNet all.chain $TARGET.sizes $QUERY.sizes all.pre.chain

## step 3

chainNet all.pre.chain -minSpace=$CHAINNET_MIN_SPACE minScore=$CHAINNET_MIN_SCORE $TARGET.sizes $QUERY.sizes target.net query.net

## step 4

netSyntenic target.net target.syn.net
netSyntenic query.net query.syn.net

