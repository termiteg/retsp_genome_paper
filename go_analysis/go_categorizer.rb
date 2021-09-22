    #######
    ###  go_categorizer.rb groups genes by GO category.
    ###
    ###  The input file needs to have tab-delimitted two columns,
    ###  the first for "gene ID" and the second for "GO ID".
    ###
    ###  go_categorizer.rb needs the file "go_daily-termdb.rdf-xml.gz",
    ###  which is available from the GO web site
    ###  http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz
    ###  This GO file must be placed in the same directory as this ruby script.
    ######
     
    require 'optparse'
    require 'zlib'
     
    opt = OptionParser.new
    inputfile = String.new
    infocont = nil
     
    opt.on('-i FILE', '--input FILE', 'GO ID list [required]') {|i| inputfile = i}
    opt.on('-c', '--infocont', 'Calculate information content (take much longer) [default: FALSE]') {|c| infocont = "T"}
    opt.on('-h', '--help', 'Show this message') {puts opt ; exit}
    opt.parse!(ARGV)
     
    scriptdir = File.dirname(__FILE__)
    gofilepath = File.expand_path("go_daily-termdb.rdf-xml.gz", scriptdir)
     
    f = Zlib::GzipReader.open(gofilepath)
     
    ### linking higher- and lower-level GO IDs
    acc2name = {}
    a_inc_b = Hash.new { |hash, key| hash[key] = [] }
     
    f.each("</go:term>") do |e|
     
       next if /<go:definition>OBSOLETE/ =~ e
     
       goacc = e.sub(/.*<go:accession>(.+?)<\/go:accession>.*/m, '\1')
       goname = e.sub(/.*<go:name>(.+?)<\/go:name>.*/m, '\1')
       acc2name[goacc] = goname
     
       next if goacc == "all"
     
       higher = e.scan(/(?<=<go:is_a rdf:resource="http:\/\/www.geneontology.org\/go#).+(?=")/)
       higher << e.scan(/(?<=<go:part_of rdf:resource="http:\/\/www.geneontology.org\/go#).+(?=")/)
     
       higher_f = higher.flatten
       higher_f.each do |a|
          next if goacc == a
          a_inc_b[a] << goacc
       end
     
    end
     
    f.close
     
     
    ### making a hash of {"a GO ID" => "all its descendant GO IDs"}
    ### counting a number of descendants
     
    allpath = Hash.new { |hash, key| hash[key] = [] }
    a_inc_b.each_key do |key|
       allpath[key] = []
    end
     
    allpath.each_key do |key|
       i = 1 ; x = 0
       allpath[key][x] = a_inc_b[key]
       while i > 0
          array = []
          allpath[key][x].each do |down|
             t = a_inc_b[down]
             next if t.size == 0
             array << t
          end
          unless array.size == 0
             allpath[key] << array.flatten
          end
          i = array.size
          x += 1
       end
    end
     
     
    ###  print GO IDs and their descendant GO IDs in a new file,
    ###  which is used for GO Partition analysis (go_partition.rb).
    fname = inputfile.sub(/.+\/(.+)\..*/, '\1')
    path_file = File.open("go_daily-termdb.rdf-xml.path", "w")
     
    allpath.each do |k, v|
       des = v.flatten.uniq
       d = ""
       d << des.shift(5).join(" ") + "\n" while des.size > 0
       path_file.print("<goterm>\n<goid>#{k}</goid>\n<descendants>\n")
       path_file.print(d)
       path_file.print("</descendants>\n</goterm>\n\n")
    end
     
     
    ### making a hash of {"GO ID" => "seq IDs"}
    go2seq = Hash.new { |hash, key| hash[key] = [] }
    allseq = []
    File.open(inputfile).each do |line|
       a = line.chomp.split("\t")
       go2seq[a[1]] << a[0]
       allseq << a[0]
    end
     
     
    ### making a table of GO ID, GO name, Num Genes, Infomation Content and Gene ID
     
    a_inc_b["all"].each do |top|
     
       allpath[top].flatten.uniq.each do |e| ## [[level_2], [level_3], ..., [level_n]].flatten
     
          seqid = go2seq[e]
     
          allpath[e].flatten.each { |d| seqid << go2seq[d] }
          nr = seqid.flatten.uniq
          n = nr.size
     
          unless n == 0
     
             if acc2name[top] == "biological_process"
                cate = "BP"
             elsif acc2name[top] == "cellular_component"
                cate = "CC"
             elsif acc2name[top] == "molecular_function"
                cate = "MF"
             end
     
             print ">" + e + "\t" + acc2name[e] + "\t" + cate + "\t" + n.to_s
     
             if infocont
                include Math
                pVn = n.to_f / allseq.uniq.size.to_f
                ic = - log2(pVn)
                print "\t" + ic.to_s + "\n"
             else
                print "\n"
             end
     
             print nr.shift(5).join(" ") + "\n" while nr.size > 0
     
          end
     
       end
     
    end

