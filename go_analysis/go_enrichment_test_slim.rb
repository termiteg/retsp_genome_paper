#=== conf ===
go_categorizer_out = "GOCat_RspeOGS1.0.renamed.out"
target_f = ARGV[0]
#"tmp.target.txt"
bg_f = ARGV[1]
"tmp.bg.txt"
max = (ARGV[2] || 1000000).to_i
min = (ARGV[3] || 5).to_i
goslimf = (ARGV[4] || "goslim_generic.obo")

outbase = File.basename(target_f, ".txt") + ".gomat_slim"

#===

## read go_slim data

slimgos = {}
File.open(goslimf).each do |l|
  if m = /^id: (GO:\d+)/.match(l)
    go = m[1]
    slimgos[go] = true
  end
end

#slimgos

## read go_categorizer_out

go2genes = {}
go_info = {}
go_curr = nil
File.open(go_categorizer_out).each do |l|
  if m = /^>(GO:\d+)/.match(l)
    go_curr = m[1]
    raise if go2genes.has_key?(go_curr)
    go2genes[go_curr] = []
    a = l.chomp.split(/\t/)
    go_info[go_curr] = {
      "id" => go_curr,
      "desc" => a[1],
      "category" => a[2]}
  else
    go2genes[go_curr] += l.chomp.split
  end

end

go2genes

gene2gos = {}
go2genes.each do |go, genes|
  genes.each do |g|
    unless gene2gos.has_key?(g)
      gene2gos[g] = []
    end
    gene2gos[g] << go
  end
end

#p gene2gos

##


genes_target = File.open(target_f).readlines.map{|x| x.chomp.strip}
#p genes_target.length
gos_target = []
genes_target.each do |gene|
  unless gene2gos.has_key?(gene)
    gene2gos[gene] = []
  end
  gene2gos[gene].each do |go|
    gos_target << go
  end
end
#p gos_target.size

genes_bg = File.open(bg_f).readlines.map{|x| x.chomp.strip}
#p genes_bg.length
gos_bg = []
genes_bg.each do |gene|
  unless gene2gos.has_key?(gene)
    gene2gos[gene] = []
  end
  gene2gos[gene].each do |go|
    gos_bg << go
  end
end
# p gos_bg.size

gos_target_count = Hash.new(0)
gos_target.each do |go|
  gos_target_count[go] += 1
end

#p gos_target_count

gos_bg_count = Hash.new(0)
gos_bg.each do |go|
  gos_bg_count[go] += 1
end

#p gos_bg_count

%w{BP CC MF}.each do |cat|

  outf = "#{outbase}.#{cat}.txt"
  o = File.open(outf, "w") 

  o.puts "# category: #{cat}"
  o.puts "# target_file: #{target_f}"
  o.puts "# background_file: #{bg_f}"
  o.puts "# num_genes_target: #{genes_target.size}"
  o.puts "# num_genes_background: #{genes_bg.size}"
  o.puts "# GO terms filter"
  o.puts "#    max: #{max}"
  o.puts "#    min" #{min}"
  o.puts "#    GO slims: yes"
  o.puts "# GO slim terms: #{slimgos.keys.size}"
  o.puts "# GO slim definition: #{goslimf}"
  o.puts "#"
  o.puts "#gene\ttarget\tbg"

  gos = gos_bg_count.keys.select{|go| go_info[go]["category"] == cat}


  gos.sort.each do |go|
    if gos_bg_count[go] <= max && gos_bg_count[go] >= min && slimgos[go]
      o.puts [go, gos_target_count[go], gos_bg_count[go]].join("\t")
    else
    end
  end

  o.close
  matrix_file = outf
  output_file = File.basename(matrix_file, ".txt") + ".fisher.Rout"

  rcode = <<"EOS"
d <- read.delim("#{matrix_file}", comment.char="#", head=F, row.names=1)
colnames(d) <- c("target", "bg")
pval <- apply(d, 1, 
  function(x){(fisher.test(matrix(c(x[1], #{genes_target.size}-x[1], x[2], #{genes_bg.size}-x[2]), 2), alternative="greater" ) )$p.value})
qval <- p.adjust(pval)
d <- cbind(d, pvalue=pval, FDR=qval)
write.table(d, "#{output_file}", sep="\t", quote=F)
EOS
  rcode_f = File.basename(matrix_file, ".Rout") + ".R"
  File.open(rcode_f, "w"){|o|
    o.puts rcode
  }
  system "Rscript --slave --vanilla #{rcode_f}"

  rout = output_file
  outputf = File.basename(rout, ".Rout") + ".txt"
  o = File.open(outputf, "w")

  ## evaluate output
  File.open(rout).each do |l|
    unless /^GO/.match(l)
      o.puts l
      next
    end
    a = l.chomp.split(/\t/)
    if a[4].to_f < 0.01
      sig = "**"
    elsif a[4].to_f < 0.05
      sig = "*"
    else
      sig = nil
    end
    out = [a, sig, go_info[a[0]]["desc"]].flatten
    o.puts out.join("\t")
  end

  o.close

end


