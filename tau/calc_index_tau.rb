
# Calculate tissue specificity index, tau (Yanai et al. 2005)
def calc_tissue_specificity_index(values, offset)
  # values : Array
  # offset : add this number to each value
  values = values.map{|v| v + offset}
  max = values.sort[-1]
  raise "max should be > 0" if max == 0
  numerator = values.map{|v| 1 - v / max}.inject(0){|sum, x| sum += x}
  denominator = (values.size - 1)
  tau = numerator / denominator
end

if __FILE__ == $0
  inputf = ARGV[0]
  offset = (ARGV[1] || 1).to_f
  pretty_out = true

  ## print header
  puts "# Calculate tissue specificity index, tau, according to Yanai's method (Yanai et al. 2005)"
  puts "# source: #{inputf}"
  puts "# offset: #{offset}"
  puts "# date:   #{Time.now}"
  puts "# script: #{$0} (#{File.expand_path($0)})"
  puts "#"

  ## main part
  File.open(inputf).each do |l|
      a = l.chomp.split(/\t/)
    if /\#/.match(l)
      #    puts l

      header = []
      header << "Name" << "tau" << a[1..-1]
        puts "#" + header.join("\t")
      next
    end
    
    id = a.shift
    values = a.map{|b| b.to_f}
    values = values.map{|v| Math.log2(v + 1)} # rpkm + 1 to
    tau = calc_tissue_specificity_index(values, offset)
    if pretty_out
      out_values = a.flatten.map{|x| sprintf("%.2f", x)}.join("\t")
    else
      out_values = a.flatten
    end
    puts [id, tau, out_values].join("\t")

  end

end

