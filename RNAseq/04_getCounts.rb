file = ARGV[0]
 
File.open(file).each do |line|
  next if line =~ /^#/
  x = line.sub(/^(.+?\t){6}/, "")
  if line =~ /^Geneid/
    s = x.gsub(/(\.\.\/thout\/)|(\/accepted_hits.bam)/, "")
    puts s
  else
    print line.slice(/^(.+?)\t/, 1) + "\t"
    puts x
  end
end
