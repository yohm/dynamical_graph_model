require 'pp'

# take logarithmic binning
file = ARGV[0]

counts = Hash.new() {|hash,key| hash[key] = [] }
File.open(file).each do |line|
  dat = line.split.map(&:to_i)
  key = dat[0]
  next if key == 0
  bin_value = 2**((Math.log(key)/Math.log(2.0)).to_i)
  counts[bin_value] = dat[1..-1].zip( counts[bin_value] ).map {|x,y| x + y.to_i }
end

counts.sort_by {|key,val| key}.each do |key,val|
  mapped = val.map {|v| v.to_f / key }
  $stdout.puts "#{key} #{mapped.join(' ')}"
end

