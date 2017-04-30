require "bundler/gem_tasks"

task :doc do
  dir = "ext/numo/ffte"
  src = %w[ffte.c]
  path = src.map{|s| File.join(dir,s)}
  sh "cd #{dir}; ruby extconf.rb; make #{src.join(' ')}"
  sh "yard doc -o yard -r README.md #{path.join(' ')}"
end

task :cleandoc do
  rm_rf %w[yard .yardoc]
end
