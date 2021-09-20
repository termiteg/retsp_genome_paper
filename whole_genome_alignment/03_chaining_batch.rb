#=== conf ===
lastz_results_dir = "lastz_results"

conf = <<'EOS'
TARGET_DIR=Rspe
QUERY_DIR=Mnat

AXTCHAIN_LINEAR_GAP=medium
AXTCHAIN_MIN_SCORE=5000

CHAINNET_MIN_SPACE=100
CHAINNET_MIN_SCORE=10000

EOS
#===

## step 1

Dir["#{lastz_results_dir}/*.psl"].each do |psl|
  m = /run_lastz_(.+)\.sh.o.psl/.match(psl)
  p  qname = m[1]
  scriptf = "run_chain_#{qname}.sh"
  File.open(scriptf, "w") do |o|
    o.puts "#=== conf ==="
    o.puts conf
    o.puts "#==="
    out = "#{psl}.chain"
    cmd = "axtChain #{psl} $TARGET_DIR $QUERY_DIR stdout -linearGap=$AXTCHAIN_LINEAR_GAP -minScore=$AXTCHAIN_MIN_SCORE -psl |  chainAntiRepeat  $TARGET_DIR $QUERY_DIR stdin #{out}; "
    o.puts cmd
  end
end
  
