file = ARGV[0]
 
File.open(file).each do |line|
  next unless line =~ /exon/
  f = line.slice(/((.+?\t){8})/, 1)
  att = line.chomp.sub(/((.+?\t){8})/, "")
  gid = att.slice(/ID=(RS\d+)/, 1)
  enum = att.slice(/\.exon(\d+)/, 1)
  print f
  print "gene_id \"#{gid}\"; exon_number \"#{enum}\"; gene_name \"#{gid}\";\n"
end
