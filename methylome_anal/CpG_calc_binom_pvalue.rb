#=== conf ===
FALSE_POS = 0.006
#===

source = ARGV[0]

#source = "tmp.CHG.txt"
#compile_mpileup_good_C_G_sites_ver3_BOTH.rb.out"

#chromosomepositionstrandcontextcontext3%methylatedsum_methylatedsum_unmethylatedcoveragedata_0data_1data_2data_3data_4data_5data_6data_7
#scaffold_0 72 + CHG CTG 0.0 0 21 21 0,3 0,3 0,0 0,1 0,4 0,6 0,3 0,1
#scaffold_0 74 - CHG CAG NaN 0 0  0  0,0 0,0 0,0 0,0 0,0 0,0 0,0 0,0

#require 'rinruby'



class Pbinom

  def initialize(p)
    @cash = {}
    @prob = p
  end

  attr_reader :prop

  def calc(x, y)
    # @cash
    if @cash.has_key?([x, y])
#      STDERR.puts "Already cashed #{x} #{y}"
      return @cash.fetch([x, y])
    elsif y == 0
      return "NA"
    else
      cmd = "Rscript --vanilla --slave -e 'binom.test(#{x}, #{y}, #{@prob}, alternative=\"greater\")$p.value'"
      rout = nil
      IO.popen(cmd){|io| rout =io.read}
      res = rout.sub(/\[1\]/, "").strip.to_f
      @cash[[x, y]] = res
      return res
    end
  end

end

dbinom = Pbinom.new(FALSE_POS)

File.open(source).each do |l|
  if /^\#/.match(l)
    puts l
    next
  end
  a = l.chomp.split(/\t/)
  chr = a[0]
  pos = a[1].to_i
  strand = a[2]
  context  = a[3]
  context3 = a[4]
  mC_perc  = a[5]
  mC = a[6].to_i
  nonmC = a[7].to_i
  cov = a[8].to_i
  pvalue = dbinom.calc(mC, cov)
  puts [a, pvalue].flatten.join("\t")
end
