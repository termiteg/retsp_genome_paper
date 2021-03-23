require 'stringio'

$query_fastaf = "Rspe_OGS1.0.pep.fa"
$odb8_gene_file ="ODB8_EukOGs_genes_Arthropoda-6656.txt"
$odb_fastaf = "ODB8.Arthropoda.fa"
$candidate_table = "Rspe_OGS1.0.pep.fa_ALL.fasta.vs.ODB8.Arthropoda.fa.blast.top.OGID.wHitsOnly.txt"

$ogid_to_genes = Hash.new

File.open($odb8_gene_file).each do |l|
  a = l.chomp.split(/\t/, -1)
  ogid = a[1]
  pid = a[3]
  unless $ogid_to_genes.has_key?(ogid)
    $ogid_to_genes[ogid] = []
  end
  $ogid_to_genes[ogid] << pid
end


def evaluate_membership(query_id, query_ogid, tophit_id, tophit_bitscore)

  tophit_bitscore = tophit_bitscore.to_f

  cmd = "blastdbcmd -db #{$odb_fastaf} -entry #{$ogid_to_genes[query_ogid].join(',')}"

  res1 = nil
  IO.popen(cmd){|io|
    res1 = io.read
  }

  res2 = nil
  cmd2 = "blastdbcmd -db #{$query_fastaf} -entry #{query_id}"
  IO.popen(cmd2){|io|
    res2 = io.read
  }

  groupfasta = res1 + "\n" + res2

  groupfasta_file = "tmp.#{query_ogid}-#{query_id}.groupfasta.pep.fa"
  File.open(groupfasta_file, "w"){|o|
    o.puts groupfasta
  }

  ## build blastdb of target OG group
  cmd3= "makeblastdb -in #{groupfasta_file} -dbtype prot -logfile makeblastdb.log"
  system cmd3

  ## retrieve tophit seq
  cmd5 = "blastdbcmd -db #{$odb_fastaf} -entry #{tophit_id}"

  tophitseq = nil
  IO.popen(cmd5){|io|
    tophitseq = io.read
  }

  cmd4 = "blastp -db #{groupfasta_file} -num_alignments 10000 -num_descriptions 10000 -seg yes -soft_masking yes -outfmt 6 "
  
  reverse_blast_output = nil
  IO.popen(cmd4, "r+") do |io|
    io.puts tophitseq
    io.close_write
    reverse_blast_output =  io.read
  end

  cmd = "rm #{groupfasta_file}.p* "
  system cmd

  ## parse reverse_blast_output

  STDERR.puts reverse_blast_output

  hitscores = []
  StringIO.open(reverse_blast_output).each do |l|
    a = l.chomp.split(/\t/)
    hitid = a[1].sub(/^lcl\|/, '')
    bitscore = a[11].strip.to_f
  
    prev = hitscores[-1]
    if prev == nil
      hitscores << [hitid, bitscore]
    elsif prev[0] == hitid && prev[1] > bitscore
      # do nothing
    else
      hitscores << [hitid, bitscore]
    end
  end

  query_record = hitscores.find{|h| h[0] == query_id}
  if query_record
    query_score = query_record[1]
  else
    ## nohit
    STDERR.puts "#{query_id} is not found in reverse blast result."
    result = [query_id, query_ogid, tophit_id, $ogid_to_genes[query_ogid].size, nil, nil, nil, nil, false]
    return result
  end

  query_at = hitscores.index{|x| x[0] == query_id}
  scores_wo_query = hitscores.map{|h| h[1]}
  scores_wo_query.delete_at(query_at)
  minscore = scores_wo_query.sort.first

  passed = false
  if query_score >= minscore
    passed = true
  end

  meanscore = hitscores.map{|h| h[1]}.sort.inject(0){|i, j| j += i} / $ogid_to_genes[query_ogid].size.to_f
  medianscore = hitscores.map{|h| h[1]}.sort[$ogid_to_genes[query_ogid].size / 2]

  result = [query_id, query_ogid, tophit_id, $ogid_to_genes[query_ogid].size, query_score, minscore, meanscore, medianscore, passed]
  
  return result
end

## process all entries 

File.open($candidate_table).each do |l|
  a = l.chomp.split(/\t/)
  query_id = a[0]
  query_ogid = a[7]
  bitscore = a[6].to_i
  tophit = a[2]
  result = evaluate_membership(query_id, query_ogid, tophit, bitscore)
  puts result.join("\t")
end
