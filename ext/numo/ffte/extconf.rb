require 'rbconfig.rb'
require 'numo/narray'

FFTE="ffte-6.0"

if !File.exist?(FFTE)
  system "tar xzf #{FFTE}.tgz"
end

require 'mkmf'

#$CFLAGS="-g -O0"
#$INCFLAGS = "-I../ext -I../ext/types #$INCFLAGS"

$LOAD_PATH.each do |x|
  if File.exist? File.join(x,'numo/numo/narray.h')
    $INCFLAGS = "-I#{x}/numo " + $INCFLAGS
    break
  end
end

$objs = (%w[
kernel
factor
fft235
zfft1d
zfft2d
zfft3d
zdfft2d
zdfft3d
dzfft2d
dzfft3d
].map{|x| File.join(FFTE,x)}+
%w[
ffte
]).map{|x| x+".o"}

fflags = ""
# GNU FORTRAN v4
if have_library("gfortran")
  $defs.push "-fPIC -DGNU_FORTRAN"
  fc = "gfortran"
  if false # have_library("gomp")
    fflags += "-fopenmp"
  end
# GNU FORTRAN v3
elsif have_library("g77")
  $defs.push "-fPIC -DGNU_FORTRAN"
  fc = "g77"
elsif have_library('f2c')
  $defs.push "-DF2C"
else
  puts "Fortran compiler not found"
  exit(1)
end

if !have_header('numo/narray.h')
  print <<EOL
  Header numo/narray.h was not found. Give pathname as follows:
  % ruby extconf.rb --with-narray-include=narray_h_dir
EOL
  exit(1)
end

create_makefile('numo/ffte')

if $makefile_created
  puts "appending extra install tasks to Makefile"
  File.open("Makefile","a") do |w|
    w.print <<EOL

F77 = #{fc}
FFLAGS = -O3 -fPIC -Iffte-6.0 #{fflags} -fomit-frame-pointer
EOL
  end
end
